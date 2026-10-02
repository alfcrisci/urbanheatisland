#' Authenticate with Sentinel Hub (OAuth2 client credentials)
#'
#' Sentinel Hub uses OAuth2 client-credentials authentication. Create an
#' OAuth client in the Sentinel Hub Dashboard
#' (\url{https://apps.sentinel-hub.com/dashboard/}) to obtain a client ID
#' and secret.
#'
#' @param client_id Sentinel Hub OAuth client ID. If \code{NULL}, read from
#'   the \code{SH_CLIENT_ID} environment variable.
#' @param client_secret Sentinel Hub OAuth client secret. If \code{NULL},
#'   read from the \code{SH_CLIENT_SECRET} environment variable.
#' @return A list with \code{token} (character) and \code{expires_at} (POSIXct).
#' @export
shub_authenticate <- function(client_id = NULL, client_secret = NULL) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("Package 'httr2' is required for Sentinel Hub functions. Install it with ",
         "install.packages('httr2').", call. = FALSE)
  }
  client_id <- client_id %||% Sys.getenv("SH_CLIENT_ID")
  client_secret <- client_secret %||% Sys.getenv("SH_CLIENT_SECRET")

  if (identical(client_id, "") || identical(client_secret, "")) {
    stop(
      "Sentinel Hub credentials not supplied. Pass client_id/client_secret ",
      "or set SH_CLIENT_ID / SH_CLIENT_SECRET environment variables.",
      call. = FALSE
    )
  }

  req <- httr2::request("https://services.sentinel-hub.com/oauth/token") |>
    httr2::req_body_form(
      grant_type = "client_credentials",
      client_id = client_id,
      client_secret = client_secret
    )

  resp <- httr2::req_perform(req)
  body <- httr2::resp_body_json(resp)

  list(
    token = body$access_token,
    expires_at = Sys.time() + as.numeric(body$expires_in)
  )
}

#' Search the Sentinel Hub Catalog API for Landsat 8/9 scenes
#'
#' @param token Bearer token from \code{\link{shub_authenticate}}.
#' @param aoi An \code{sf} polygon (EPSG:4326) defining the search area.
#' @param start_date,end_date Character dates, \code{"YYYY-MM-DD"}.
#' @param collection Sentinel Hub collection ID. Default
#'   \code{"landsat-ot-l2"} (Landsat 8/9 OLI-TIRS Collection 2 Level-2).
#'   Use \code{"landsat-ot-l1"} for Level-1.
#' @param max_cloud_cover Maximum scene cloud cover percentage. Default 20.
#' @param limit Maximum number of results. Default 10.
#' @return A data.frame with one row per scene: \code{id}, \code{datetime},
#'   \code{cloud_cover}, ordered by increasing cloud cover.
#' @export
shub_catalog_search <- function(token, aoi, start_date, end_date,
                                 collection = "landsat-ot-l2",
                                 max_cloud_cover = 20, limit = 10) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("Package 'httr2' is required for Sentinel Hub functions.", call. = FALSE)
  }
  bbox <- as.numeric(sf::st_bbox(sf::st_transform(aoi, 4326)))

  body <- list(
    bbox = bbox,
    datetime = paste0(start_date, "T00:00:00Z/", end_date, "T23:59:59Z"),
    collections = list(collection),
    limit = limit,
    filter = sprintf("eo:cloud_cover < %s", max_cloud_cover)
  )

  req <- httr2::request("https://services.sentinel-hub.com/api/v1/catalog/1.0.0/search") |>
    httr2::req_auth_bearer_token(token) |>
    httr2::req_body_json(body)

  resp <- httr2::req_perform(req)
  result <- httr2::resp_body_json(resp)

  if (length(result$features) == 0) {
    return(data.frame(id = character(0), datetime = character(0), cloud_cover = numeric(0)))
  }

  df <- do.call(rbind, lapply(result$features, function(f) {
    data.frame(
      id = f$id,
      datetime = f$properties$datetime,
      cloud_cover = f$properties$`eo:cloud_cover` %||% NA_real_
    )
  }))
  df[order(df$cloud_cover), ]
}

#' Default evalscript for Landsat thermal (ST_B10) retrieval
#'
#' Returns a Sentinel Hub evalscript that outputs the Landsat Collection 2
#' Level-2 surface temperature band (ST_B10), already converted to degrees
#' Celsius using the official scale (mult 0.00341802, add 149.0, then
#' Kelvin-to-Celsius).
#'
#' @return Character evalscript (Sentinel Hub evalscript v3).
#' @export
shub_evalscript_lst <- function() {
  "//VERSION=3
function setup() {
  return {
    input: [{ bands: [\"ST_B10\"] }],
    output: { bands: 1, sampleType: \"FLOAT32\" }
  };
}
function evaluatePixel(sample) {
  // ST_B10 raw DN -> Kelvin -> Celsius (Landsat Collection 2 Level-2 scale factors)
  let kelvin = sample.ST_B10 * 0.00341802 + 149.0;
  return [kelvin - 273.15];
}"
}

#' Request a Landsat LST raster via the Sentinel Hub Process API
#'
#' @param token Bearer token from \code{\link{shub_authenticate}}.
#' @param aoi An \code{sf} polygon (EPSG:4326) defining the request area.
#' @param datetime A single scene datetime string (e.g. from
#'   \code{\link{shub_catalog_search}}'s \code{datetime} column), or a
#'   \code{"start/end"} range string.
#' @param collection Sentinel Hub collection ID. Default \code{"landsat-ot-l2"}.
#' @param resolution Output pixel size in metres. Default 30 (native Landsat).
#' @param evalscript Evalscript string. Default \code{\link{shub_evalscript_lst}()}.
#' @param dest File path to write the output GeoTIFF to. Default a temp file.
#' @return A single-layer \code{SpatRaster} (degrees Celsius, with the
#'   default evalscript).
#' @export
shub_request_lst <- function(token, aoi, datetime, collection = "landsat-ot-l2",
                              resolution = 30, evalscript = shub_evalscript_lst(),
                              dest = tempfile(fileext = ".tif")) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("Package 'httr2' is required for Sentinel Hub functions.", call. = FALSE)
  }
  bbox <- as.numeric(sf::st_bbox(sf::st_transform(aoi, 4326)))
  width <- max(1, round((bbox[3] - bbox[1]) * 111320 / resolution))
  height <- max(1, round((bbox[4] - bbox[2]) * 111320 / resolution))

  body <- list(
    input = list(
      bounds = list(bbox = bbox, properties = list(crs = "http://www.opengis.net/def/crs/EPSG/0/4326")),
      data = list(list(type = collection, dataFilter = list(timeRange = list(
        from = if (grepl("/", datetime)) strsplit(datetime, "/")[[1]][1] else datetime,
        to   = if (grepl("/", datetime)) strsplit(datetime, "/")[[1]][2] else datetime
      ))))
    ),
    output = list(width = width, height = height,
                  responses = list(list(identifier = "default", format = list(type = "image/tiff")))),
    evalscript = evalscript
  )

  req <- httr2::request("https://services.sentinel-hub.com/api/v1/process") |>
    httr2::req_auth_bearer_token(token) |>
    httr2::req_body_json(body)

  httr2::req_perform(req, path = dest)
  terra::rast(dest)
}

#' End-to-end Landsat LST acquisition via Sentinel Hub
#'
#' Authenticates, searches the catalog for the least-cloudy scene in the
#' date range, and requests its ST_B10-derived LST raster, in one call.
#'
#' @inheritParams shub_catalog_search
#' @param client_id,client_secret Sentinel Hub OAuth credentials (see
#'   \code{\link{shub_authenticate}}).
#' @param resolution Output pixel size in metres. Default 30.
#' @return A single-layer \code{SpatRaster} in degrees Celsius.
#' @export
get_landsat_lst_sentinelhub <- function(aoi, start_date, end_date,
                                         client_id = NULL, client_secret = NULL,
                                         collection = "landsat-ot-l2",
                                         max_cloud_cover = 20, resolution = 30) {
  auth <- shub_authenticate(client_id, client_secret)
  scenes <- shub_catalog_search(auth$token, aoi, start_date, end_date,
                                 collection = collection, max_cloud_cover = max_cloud_cover)
  if (nrow(scenes) == 0) {
    stop("No Landsat scenes found for the given AOI/date range/cloud threshold.", call. = FALSE)
  }
  best <- scenes[1, ]
  message(sprintf("Using scene %s (%s, %.1f%% cloud cover)", best$id, best$datetime, best$cloud_cover))
  shub_request_lst(auth$token, aoi, best$datetime, collection = collection, resolution = resolution)
}
