#' NASA AppEEARS authentication
#'
#' Authenticate against the NASA AppEEARS API (Earthdata Login credentials
#' required, see \url{https://urs.earthdata.nasa.gov/}) and obtain a bearer
#' token used by all other \code{lst_appeears_*} functions.
#'
#' @param username Earthdata Login username. If \code{NULL}, read from the
#'   \code{EARTHDATA_USER} environment variable.
#' @param password Earthdata Login password. If \code{NULL}, read from the
#'   \code{EARTHDATA_PASS} environment variable.
#' @return A list with \code{token} (character) and \code{expiration} (POSIXct).
#' @export
lst_appeears_login <- function(username = NULL, password = NULL) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("Package 'httr2' is required for AppEEARS functions. Install it with ",
         "install.packages('httr2').", call. = FALSE)
  }
  username <- username %||% Sys.getenv("EARTHDATA_USER")
  password <- password %||% Sys.getenv("EARTHDATA_PASS")

  if (identical(username, "") || identical(password, "")) {
    stop(
      "Earthdata credentials not supplied. Pass username/password or set ",
      "EARTHDATA_USER / EARTHDATA_PASS environment variables.",
      call. = FALSE
    )
  }

  req <- httr2::request("https://appeears.earthdatacloud.nasa.gov/api/login") |>
    httr2::req_auth_basic(username, password) |>
    httr2::req_method("POST")

  resp <- httr2::req_perform(req)
  body <- httr2::resp_body_json(resp)

  list(
    token = body$token,
    expiration = as.POSIXct(body$expiration, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )
}

#' Submit an AppEEARS area request for MODIS LST
#'
#' Submits an "area" task to AppEEARS for a MODIS Land Surface Temperature
#' product (default: MOD11A2.061, 8-day 1km composite) clipped to a
#' user-supplied polygon, over a date range.
#'
#' @param token Bearer token from \code{\link{lst_appeears_login}}.
#' @param aoi An \code{sf} polygon (single feature or union) defining the area
#'   of interest, in EPSG:4326.
#' @param start_date,end_date Character dates \code{"MM-DD-YYYY"}.
#' @param product MODIS product + version string. Default \code{"MOD11A2.061"}
#'   (Terra, 8-day). Use \code{"MYD11A2.061"} for Aqua, or \code{"MOD11A1.061"}
#'   for daily.
#' @param layers Character vector of layers to request. Default requests day
#'   and night LST plus their QC layers.
#' @param task_name Character label for the task.
#' @return Character task ID.
#' @export
lst_appeears_submit_task <- function(token,
                                      aoi,
                                      start_date,
                                      end_date,
                                      product = "MOD11A2.061",
                                      layers = c("LST_Day_1km", "QC_Day",
                                                 "LST_Night_1km", "QC_Night"),
                                      task_name = paste0("LST_", format(Sys.time(), "%Y%m%d_%H%M%S"))) {

  aoi_geojson <- geojsonsf_from_sf(aoi)

  layer_list <- lapply(layers, function(l) list(layer = l, product = product))

  task_body <- list(
    task_type = "area",
    task_name = task_name,
    params = list(
      dates = list(list(startDate = start_date, endDate = end_date)),
      layers = layer_list,
      geo = aoi_geojson,
      output = list(format = list(type = "geotiff"), projection = "geographic")
    )
  )

  req <- httr2::request("https://appeears.earthdatacloud.nasa.gov/api/task") |>
    httr2::req_auth_bearer_token(token) |>
    httr2::req_body_json(task_body) |>
    httr2::req_method("POST")

  resp <- httr2::req_perform(req)
  body <- httr2::resp_body_json(resp)
  body$task_id
}

#' Check status of an AppEEARS task
#'
#' @param task_id Character task ID from \code{\link{lst_appeears_submit_task}}.
#' @param token Bearer token.
#' @return Character status: one of \code{"pending"}, \code{"processing"},
#'   \code{"done"}, \code{"error"}.
#' @export
lst_appeears_check_status <- function(task_id, token) {
  req <- httr2::request(paste0("https://appeears.earthdatacloud.nasa.gov/api/task/", task_id)) |>
    httr2::req_auth_bearer_token(token)
  resp <- httr2::req_perform(req)
  body <- httr2::resp_body_json(resp)
  body$status
}

#' Poll an AppEEARS task until completion
#'
#' @param task_id Character task ID.
#' @param token Bearer token.
#' @param poll_interval Seconds between status checks. Default 60.
#' @param timeout Maximum seconds to wait before giving up. Default 7200 (2h).
#' @param verbose Print status updates.
#' @return Invisibly, the final status string.
#' @export
lst_appeears_wait <- function(task_id, token, poll_interval = 60,
                               timeout = 7200, verbose = TRUE) {
  start <- Sys.time()
  repeat {
    status <- lst_appeears_check_status(task_id, token)
    if (verbose) message(sprintf("[%s] task %s: %s", format(Sys.time(), "%H:%M:%S"), task_id, status))
    if (status %in% c("done", "error")) return(invisible(status))
    if (as.numeric(difftime(Sys.time(), start, units = "secs")) > timeout) {
      stop("Timed out waiting for AppEEARS task ", task_id, call. = FALSE)
    }
    Sys.sleep(poll_interval)
  }
}

#' Download completed AppEEARS task output
#'
#' Downloads all GeoTIFF bundle files for a completed task to a local
#' directory.
#'
#' @param task_id Character task ID.
#' @param token Bearer token.
#' @param outdir Directory to write files to. Created if it doesn't exist.
#' @return Character vector of downloaded file paths, invisibly.
#' @export
lst_appeears_download <- function(task_id, token, outdir) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  bundle_req <- httr2::request(paste0("https://appeears.earthdatacloud.nasa.gov/api/bundle/", task_id)) |>
    httr2::req_auth_bearer_token(token)
  bundle <- httr2::resp_body_json(httr2::req_perform(bundle_req))

  files <- Filter(function(f) grepl("\\.tif$", f$file_name), bundle$files)
  if (length(files) == 0) stop("No GeoTIFF files found in bundle for task ", task_id, call. = FALSE)

  out_paths <- character(length(files))
  for (i in seq_along(files)) {
    f <- files[[i]]
    dest <- file.path(outdir, basename(f$file_name))
    dl_req <- httr2::request(
      paste0("https://appeears.earthdatacloud.nasa.gov/api/bundle/", task_id, "/", f$file_id)
    ) |> httr2::req_auth_bearer_token(token)
    httr2::req_perform(dl_req, path = dest)
    out_paths[i] <- dest
  }
  invisible(out_paths)
}

#' End-to-end MODIS LST acquisition
#'
#' Convenience wrapper: login, submit, wait, download, and return a
#' processed LST raster (see \code{\link{read_lst_local}}) in a single call.
#'
#' @inheritParams lst_appeears_submit_task
#' @param username,password Earthdata credentials (see \code{\link{lst_appeears_login}}).
#' @param outdir Local directory to store downloaded files.
#' @param qc_filter Whether to mask low-quality pixels using the QC layer.
#'   Default \code{TRUE}.
#' @return A \code{SpatRaster} (or stack if multiple dates), in degrees Celsius.
#' @export
get_modis_lst <- function(aoi, start_date, end_date,
                           username = NULL, password = NULL,
                           product = "MOD11A2.061",
                           outdir = tempfile("lst_"),
                           qc_filter = TRUE) {
  auth <- lst_appeears_login(username, password)
  task_id <- lst_appeears_submit_task(auth$token, aoi, start_date, end_date, product = product)
  lst_appeears_wait(task_id, auth$token)
  files <- lst_appeears_download(task_id, auth$token, outdir)

  day_files <- grep("LST_Day", files, value = TRUE)
  qc_files  <- grep("QC_Day", files, value = TRUE)

  if (length(day_files) == 0) stop("No LST_Day files found among downloads.", call. = FALSE)

  lst_stack <- terra::rast(day_files)
  lst_stack <- lst_stack * 0.02 - 273.15  # MODIS LST scale factor -> Kelvin -> Celsius

  if (qc_filter && length(qc_files) > 0) {
    qc_stack <- terra::rast(qc_files)
    good_mask <- decode_modis_qc(qc_stack)
    lst_stack <- terra::mask(lst_stack, good_mask, maskvalues = 0)
  }

  names(lst_stack) <- tools::file_path_sans_ext(basename(day_files))
  lst_stack
}

#' @keywords internal
`%||%` <- function(a, b) if (is.null(a) || identical(a, "")) b else a

#' @keywords internal
geojsonsf_from_sf <- function(aoi) {
  gj <- sf::st_as_sf(aoi) |> sf::st_geometry() |> sf::st_as_sfc() |> sf::st_as_text()
  # AppEEARS expects a GeoJSON FeatureCollection
  fc <- list(
    type = "FeatureCollection",
    features = list(list(
      type = "Feature",
      properties = list(),
      geometry = jsonlite::fromJSON(sf_geometry_to_geojson(sf::st_geometry(aoi)), simplifyVector = FALSE)
    ))
  )
  fc
}

#' @keywords internal
sf_geometry_to_geojson <- function(geom) {
  # thin wrapper kept separate so it can be swapped for geojsonsf::sfc_geojson()
  # if that package is available, without changing the public API.
  if (requireNamespace("geojsonsf", quietly = TRUE)) {
    geojsonsf::sfc_geojson(geom)
  } else {
    sf::st_write(sf::st_sf(geometry = geom), dsn = tempfile(fileext = ".geojson"), quiet = TRUE)
  }
}
