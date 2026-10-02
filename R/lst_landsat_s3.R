#' Search the USGS Landsat STAC API
#'
#' Searches the official USGS Landsat Collection 2 STAC catalog for scenes
#' intersecting an area of interest and date range. Returns public HTTPS
#' asset URLs (served via landsatlook.usgs.gov) that do not require AWS
#' credentials — this is generally the easiest way to get individual
#' Landsat bands without setting up requester-pays S3 access.
#'
#' \strong{Implementation note}: on the native USGS STAC server, surface
#' reflectance bands (SR_B4/B5/B6, i.e. red/NIR/SWIR) and the thermal band
#' (ST_B10) live in \emph{two separate collections} —
#' \code{"landsat-c2l2-sr"} and \code{"landsat-c2l2-st"} respectively — even
#' though they come from the same acquisition. A search against only one of
#' them will never return the other's bands. This function queries both
#' collections and merges the results by scene (matching on the shared
#' acquisition identifier, stripping the collection-specific \code{_SR}/
#' \code{_ST} suffix), so every row has both thermal and spectral asset URLs
#' populated whenever available.
#'
#' @param aoi An \code{sf} polygon (EPSG:4326) defining the search area.
#' @param start_date,end_date Character dates, \code{"YYYY-MM-DD"}.
#' @param max_cloud_cover Maximum scene cloud cover percentage. Default 20.
#' @param stac_url STAC search endpoint. Default the USGS LandsatLook STAC server.
#' @param limit Maximum number of results \emph{per collection}. Default 10.
#' @return A data.frame with one row per scene: \code{id} (base scene
#'   identifier, collection suffix stripped), \code{datetime},
#'   \code{cloud_cover}, \code{st_b10_href} (thermal band asset URL),
#'   \code{qa_pixel_href} (QA band asset URL), \code{sr_b4_href} (red,
#'   surface reflectance), \code{sr_b5_href} (NIR), \code{sr_b6_href}
#'   (SWIR-1) — the latter three enable NDVI/NDBI covariates via
#'   \code{\link{compute_ndvi}}/\code{\link{compute_ndbi}}. Ordered by
#'   increasing cloud cover. A scene may have \code{NA} in some href
#'   columns if it was only found in one of the two collections (e.g. an ST
#'   product not yet available for a very recent acquisition).
#' @export
landsat_stac_search <- function(aoi, start_date, end_date,
                                 max_cloud_cover = 20,
                                 stac_url = "https://landsatlook.usgs.gov/stac-server/search",
                                 limit = 10) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("Package 'httr2' is required for STAC search. Install it with ",
         "install.packages('httr2').", call. = FALSE)
  }
  bbox <- as.numeric(sf::st_bbox(sf::st_transform(aoi, 4326)))

  sr_df <- .stac_search_one_collection(stac_url, bbox, start_date, end_date,
                                        "landsat-c2l2-sr", max_cloud_cover, limit)
  st_df <- .stac_search_one_collection(stac_url, bbox, start_date, end_date,
                                        "landsat-c2l2-st", max_cloud_cover, limit)

  .merge_landsat_scenes(sr_df, st_df)
}

#' @keywords internal
.stac_search_one_collection <- function(stac_url, bbox, start_date, end_date,
                                         collection, max_cloud_cover, limit) {
  body <- list(
    bbox = bbox,
    datetime = paste0(start_date, "T00:00:00Z/", end_date, "T23:59:59Z"),
    collections = list(collection),
    query = list(`eo:cloud_cover` = list(lt = max_cloud_cover)),
    limit = limit
  )

  req <- httr2::request(stac_url) |> httr2::req_body_json(body)
  resp <- httr2::req_perform(req)
  result <- httr2::resp_body_json(resp)

  empty <- data.frame(base_id = character(0), datetime = character(0),
                       cloud_cover = numeric(0), st_b10_href = character(0),
                       qa_pixel_href = character(0), sr_b4_href = character(0),
                       sr_b5_href = character(0), sr_b6_href = character(0),
                       stringsAsFactors = FALSE)
  if (length(result$features) == 0) return(empty)

  do.call(rbind, lapply(result$features, function(f) {
    assets <- f$assets
    data.frame(
      base_id = .strip_scene_suffix(f$id),
      datetime = f$properties$datetime,
      cloud_cover = f$properties$`eo:cloud_cover` %||% NA_real_,
      st_b10_href = .find_asset_href(assets, c("ST_B10", "st_b10", "lwir11")) %||% NA_character_,
      qa_pixel_href = .find_asset_href(assets, c("QA_PIXEL", "qa_pixel")) %||% NA_character_,
      sr_b4_href = .find_asset_href(assets, c("SR_B4", "sr_b4", "red")) %||% NA_character_,     # rosso
      sr_b5_href = .find_asset_href(assets, c("SR_B5", "sr_b5", "nir08")) %||% NA_character_,   # NIR
      sr_b6_href = .find_asset_href(assets, c("SR_B6", "sr_b6", "swir16")) %||% NA_character_,  # SWIR-1
      stringsAsFactors = FALSE
    )
  }))
}

#' @keywords internal
.strip_scene_suffix <- function(id) {
  sub("_(SR|ST)$", "", id)
}

#' Merge SR- and ST-collection scene rows into one row per acquisition
#'
#' @param sr_df,st_df Data frames from \code{.stac_search_one_collection()}
#'   for the SR and ST collections respectively (same \code{base_id} column
#'   used as the join key).
#' @return Merged data.frame, one row per unique \code{base_id}, ordered by
#'   \code{cloud_cover} ascending. Column \code{id} added as an alias of
#'   \code{base_id} for backward compatibility.
#' @keywords internal
.merge_landsat_scenes <- function(sr_df, st_df) {
  href_cols <- c("st_b10_href", "qa_pixel_href", "sr_b4_href", "sr_b5_href", "sr_b6_href")

  all_ids <- union(sr_df$base_id, st_df$base_id)
  if (length(all_ids) == 0) {
    out <- data.frame(id = character(0), datetime = character(0), cloud_cover = numeric(0))
    for (col in href_cols) out[[col]] <- character(0)
    return(out)
  }

  rows <- lapply(all_ids, function(bid) {
    sr_row <- sr_df[sr_df$base_id == bid, , drop = FALSE]
    st_row <- st_df[st_df$base_id == bid, , drop = FALSE]

    pick <- function(col) {
      sr_val <- if (nrow(sr_row) > 0) sr_row[[col]][1] else NA
      st_val <- if (nrow(st_row) > 0) st_row[[col]][1] else NA
      if (!is.na(sr_val)) sr_val else st_val
    }

    datetime <- if (nrow(sr_row) > 0) sr_row$datetime[1] else st_row$datetime[1]
    cloud_cover <- if (nrow(sr_row) > 0) sr_row$cloud_cover[1] else st_row$cloud_cover[1]

    row <- data.frame(id = bid, datetime = datetime, cloud_cover = cloud_cover,
                       stringsAsFactors = FALSE)
    for (col in href_cols) row[[col]] <- pick(col)
    row
  })

  merged <- do.call(rbind, rows)
  merged[order(merged$cloud_cover), ]
}

#' @keywords internal
.find_asset_href <- function(assets, candidate_names) {
  nm <- names(assets)
  hit <- nm[nm %in% candidate_names]
  if (length(hit) == 0) return(NULL)
  assets[[hit[1]]]$href
}

#' Download a Landsat asset from an href returned by \code{\link{landsat_stac_search}}
#'
#' \strong{Importante (aggiornamento 2025)}: il viewer LandsatLook pubblico
#' e' stato dismesso (13 marzo 2025) e \code{landsatlook.usgs.gov/data/...}
#' ora reindirizza al login ERS per le richieste non autenticate — quindi
#' questa funzione, da sola, \emph{non} scarica piu' con successo i file,
#' anche se l'URL e' corretto. Per un download che funziona senza
#' approvazione manuale dell'account, usa
#' \code{\link{get_landsat_lst_pc}} (Planetary Computer) o
#' \code{\link{get_sentinel3_lst}} (Copernicus Data Space). Questa
#' funzione resta utile per altri asset pubblici (es. thumbnail, alcuni
#' file di metadati) che potrebbero non richiedere login.
#'
#' @param url Asset URL, e.g. from \code{\link{landsat_stac_search}}'s
#'   \code{st_b10_href} / \code{qa_pixel_href} columns.
#' @param dest Destination file path. Default a temp file preserving the
#'   original extension.
#' @return \code{dest}, invisibly.
#' @export
landsat_download_asset <- function(url, dest = tempfile(fileext = paste0(".", tools::file_ext(url)))) {
  if (requireNamespace("httr2", quietly = TRUE)) {
    req <- httr2::request(url)
    httr2::req_perform(req, path = dest)
  } else {
    utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  }
  invisible(dest)
}

#' Download a Landsat Collection 2 object directly from the requester-pays S3 bucket
#'
#' Use this only if you need direct S3 access (e.g. bulk downloads at scale)
#' and have an AWS account set up for requester-pays billing. For most
#' single-scene use cases, \code{\link{landsat_download_asset}} with a URL
#' from \code{\link{landsat_stac_search}} is simpler and free of AWS setup.
#'
#' @param key S3 object key within the \code{usgs-landsat} bucket.
#' @param dest Destination file path.
#' @param bucket S3 bucket name. Default \code{"usgs-landsat"}.
#' @param region AWS region. Default \code{"us-west-2"}.
#' @return \code{dest}, invisibly.
#' @export
landsat_s3_download <- function(key, dest, bucket = "usgs-landsat", region = "us-west-2") {
  if (!requireNamespace("aws.s3", quietly = TRUE)) {
    stop("Package 'aws.s3' is required for direct S3 access. Install it with ",
         "install.packages('aws.s3'). Requires AWS credentials configured ",
         "(AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY) with requester-pays billing enabled.",
         call. = FALSE)
  }
  aws.s3::save_object(object = key, bucket = bucket, file = dest,
                       region = region, request_pays = TRUE)
  invisible(dest)
}

#' Decode the Landsat Collection 2 QA_PIXEL band
#'
#' Decodes the 16-bit QA_PIXEL band accompanying Landsat Collection 2
#' Level-2 products into a binary keep/discard mask.
#'
#' @param qa A \code{SpatRaster} of raw QA_PIXEL values.
#' @param exclude Character vector of flag bits to exclude (mask out).
#'   One or more of \code{"fill"}, \code{"dilated_cloud"}, \code{"cirrus"},
#'   \code{"cloud"}, \code{"cloud_shadow"}, \code{"snow"}. Default excludes
#'   fill, cloud, cloud shadow, cirrus, and dilated cloud (i.e. keeps clear
#'   land/water pixels, snow included).
#' @return A binary \code{SpatRaster} mask (1 = keep, 0 = discard).
#' @export
decode_landsat_qa_pixel <- function(qa, exclude = c("fill", "dilated_cloud", "cirrus",
                                                      "cloud", "cloud_shadow")) {
  bit_pos <- c(fill = 0, dilated_cloud = 1, cirrus = 2, cloud = 3,
               cloud_shadow = 4, snow = 5, clear = 6, water = 7)
  exclude_bits <- bit_pos[exclude]
  if (anyNA(exclude_bits)) stop("Unknown flag in `exclude`.", call. = FALSE)

  mask <- terra::app(qa, function(x) {
    bad <- Reduce(`|`, lapply(exclude_bits, function(b) bitwAnd(as.integer(x), bitwShiftL(1L, b)) != 0))
    as.numeric(!bad)
  })
  mask
}

#' Read and process a Landsat Collection 2 Level-2 ST_B10 (LST) file
#'
#' Applies the official Collection 2 scale/offset to convert raw ST_B10
#' digital numbers to degrees Celsius, and optionally masks using the
#' QA_PIXEL band.
#'
#' @param st_b10_path Path to the ST_B10 GeoTIFF.
#' @param qa_path Optional path to the matching QA_PIXEL GeoTIFF.
#' @param exclude See \code{\link{decode_landsat_qa_pixel}}.
#' @return A single-layer \code{SpatRaster} in degrees Celsius.
#' @export
read_landsat_lst <- function(st_b10_path, qa_path = NULL,
                              exclude = c("fill", "dilated_cloud", "cirrus",
                                          "cloud", "cloud_shadow")) {
  r <- terra::rast(st_b10_path)
  r <- r * 0.00341802 + 149.0   # -> Kelvin
  r <- r - 273.15                # -> Celsius

  if (!is.null(qa_path)) {
    qa <- terra::rast(qa_path)
    mask <- decode_landsat_qa_pixel(qa, exclude = exclude)
    r <- terra::mask(r, mask, maskvalues = 0)
  }
  names(r) <- "LST_celsius"
  r
}

#' End-to-end Landsat 8/9 LST acquisition via the public USGS STAC catalog
#'
#' Searches for the least-cloudy scene in the date range, downloads its
#' ST_B10 and QA_PIXEL assets over public HTTPS (no AWS credentials
#' required), and returns a processed LST raster in degrees Celsius.
#'
#' @inheritParams landsat_stac_search
#' @param outdir Local directory to store downloaded files. Default a temp dir.
#' @param qc_filter Whether to mask low-quality pixels using QA_PIXEL.
#'   Default \code{TRUE}.
#' @return A single-layer \code{SpatRaster} in degrees Celsius.
#' @export
get_landsat_lst <- function(aoi, start_date, end_date, max_cloud_cover = 20,
                             outdir = tempfile("landsat_"), qc_filter = TRUE) {
  scenes <- landsat_stac_search(aoi, start_date, end_date, max_cloud_cover = max_cloud_cover)
  if (nrow(scenes) == 0) {
    stop("No Landsat scenes found for the given AOI/date range/cloud threshold.", call. = FALSE)
  }
  best <- scenes[1, ]
  if (is.na(best$st_b10_href)) {
    stop("Selected scene has no ST_B10 asset available.", call. = FALSE)
  }
  message(sprintf("Using scene %s (%s, %.1f%% cloud cover)", best$id, best$datetime, best$cloud_cover))

  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  st_dest <- file.path(outdir, paste0(best$id, "_ST_B10.tif"))
  landsat_download_asset(best$st_b10_href, st_dest)

  qa_dest <- NULL
  if (qc_filter && !is.na(best$qa_pixel_href)) {
    qa_dest <- file.path(outdir, paste0(best$id, "_QA_PIXEL.tif"))
    landsat_download_asset(best$qa_pixel_href, qa_dest)
  }

  read_landsat_lst(st_dest, qa_path = qa_dest)
}
