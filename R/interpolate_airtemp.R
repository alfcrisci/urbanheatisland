#' Interpolate station air temperature to a continuous raster surface
#'
#' Unified interface to several spatial interpolation methods for producing
#' a continuous air-temperature raster from point station observations.
#' Method names and defaults mirror those used in \code{R.IBE.airqino}, so
#' outputs are directly comparable across packages.
#'
#' @param stations An \code{sf} POINT object with station locations.
#' @param value_col Character name of the column in \code{stations} holding
#'   the temperature values to interpolate.
#' @param template A \code{SpatRaster} defining the output grid (extent,
#'   resolution, CRS). Values are ignored, only geometry is used.
#' @param method One of \code{"idw"}, \code{"kriging"}, \code{"tps"} (thin
#'   plate spline), \code{"tin"}, \code{"nn"} (nearest neighbour).
#' @param idp Inverse distance power for IDW. Default \code{2}.
#' @param nmax Maximum number of neighbours used for IDW/kriging local
#'   prediction. Default \code{Inf} (global).
#' @param variogram_model Variogram model for kriging, passed to
#'   \code{gstat::vgm}. Default \code{"Sph"} (spherical).
#' @param ... Additional arguments passed through to the underlying engine
#'   (\code{gstat::gstat}, \code{fields::Tps}, or \code{terra::interpNear}).
#' @return A single-layer \code{SpatRaster} with the interpolated surface.
#' @export
interpolate_airtemp <- function(stations, value_col, template,
                                 method = c("idw", "kriging", "tps", "tin", "nn"),
                                 idp = 2, nmax = Inf,
                                 variogram_model = "Sph", ...) {
  method <- match.arg(method)

  if (!inherits(stations, "sf")) stop("`stations` must be an sf object.", call. = FALSE)
  if (!value_col %in% names(stations)) stop("`value_col` not found in `stations`.", call. = FALSE)

  switch(method,
    idw     = .interp_idw(stations, value_col, template, idp, nmax, ...),
    kriging = .interp_kriging(stations, value_col, template, variogram_model, nmax, ...),
    tps     = .interp_tps(stations, value_col, template, ...),
    tin     = .interp_tin(stations, value_col, template, ...),
    nn      = .interp_nn(stations, value_col, template, ...)
  )
}

#' @keywords internal
.interp_idw <- function(stations, value_col, template, idp, nmax, ...) {
  sp_pts <- sf::st_coordinates(stations)
  df <- data.frame(x = sp_pts[, 1], y = sp_pts[, 2], z = stations[[value_col]])
  sf_pts <- sf::st_as_sf(df, coords = c("x", "y"), crs = sf::st_crs(stations))

  grid_sf <- .grid_cells_as_sf(template, sf::st_crs(stations))

  g <- gstat::gstat(formula = z ~ 1, locations = sf_pts, set = list(idp = idp), nmax = nmax)
  pred <- stats::predict(g, newdata = grid_sf, ...)

  out <- template[[1]]
  terra::values(out) <- pred$var1.pred
  names(out) <- paste0(value_col, "_idw")
  out
}

#' @keywords internal
.interp_kriging <- function(stations, value_col, template, variogram_model, nmax, ...) {
  sp_pts <- sf::st_coordinates(stations)
  df <- data.frame(x = sp_pts[, 1], y = sp_pts[, 2], z = stations[[value_col]])
  sf_pts <- sf::st_as_sf(df, coords = c("x", "y"), crs = sf::st_crs(stations))

  vgm_emp <- gstat::variogram(z ~ 1, sf_pts)
  vgm_fit <- tryCatch(
    gstat::fit.variogram(vgm_emp, gstat::vgm(variogram_model)),
    error = function(e) gstat::vgm(psill = stats::var(sf_pts$z), model = variogram_model,
                                    range = max(vgm_emp$dist) / 3, nugget = 0)
  )

  grid_sf <- .grid_cells_as_sf(template, sf::st_crs(stations))

  g <- gstat::gstat(formula = z ~ 1, locations = sf_pts, model = vgm_fit,
                     nmax = if (is.infinite(nmax)) Inf else nmax)
  pred <- stats::predict(g, newdata = grid_sf, ...)

  out <- template[[1]]
  terra::values(out) <- pred$var1.pred
  names(out) <- paste0(value_col, "_kriging")
  attr(out, "variogram") <- vgm_fit
  out
}

#' @keywords internal
.interp_tps <- function(stations, value_col, template, ...) {
  if (!requireNamespace("fields", quietly = TRUE)) {
    stop("Package 'fields' required for method = 'tps'.", call. = FALSE)
  }
  coords <- sf::st_coordinates(stations)
  z <- stations[[value_col]]
  fit <- fields::Tps(coords, z, ...)

  full_coords <- terra::xyFromCell(template, seq_len(terra::ncell(template)))
  pred_vals <- fields::predict.Tps(fit, full_coords)

  out <- template[[1]]
  terra::values(out) <- as.numeric(pred_vals)
  names(out) <- paste0(value_col, "_tps")
  out
}

#' Build an sf POINT layer of every raster cell centre, in a target CRS
#' @keywords internal
.grid_cells_as_sf <- function(template, target_crs) {
  xy <- terra::xyFromCell(template, seq_len(terra::ncell(template)))
  grid_sf <- sf::st_as_sf(as.data.frame(xy), coords = c("x", "y"),
                           crs = terra::crs(template))
  sf::st_transform(grid_sf, target_crs)
}

#' @keywords internal
.interp_tin <- function(stations, value_col, template, ...) {
  if (!requireNamespace("terra", quietly = TRUE)) stop("terra required.", call. = FALSE)
  coords <- sf::st_coordinates(stations)
  # terra::interpNear uses nearest-neighbour by default; for a true TIN
  # (Delaunay-based linear interpolation) use terra::rasterize with a
  # linear method or interp::interp() if available.
  if (requireNamespace("interp", quietly = TRUE)) {
    li <- interp::interp(x = coords[, 1], y = coords[, 2], z = stations[[value_col]],
                          xo = terra::xFromCol(template, 1:terra::ncol(template)),
                          yo = terra::yFromRow(template, 1:terra::nrow(template)),
                          linear = TRUE, extrap = FALSE)
    out <- terra::rast(li$z[, ncol(li$z):1] |> t(), extent = terra::ext(template),
                        crs = terra::crs(template))
    out <- terra::resample(out, template[[1]])
  } else {
    warning("Package 'interp' not available; falling back to nearest-neighbour for method = 'tin'.")
    out <- .interp_nn(stations, value_col, template, ...)
  }
  names(out) <- paste0(value_col, "_tin")
  out
}

#' @keywords internal
.interp_nn <- function(stations, value_col, template, radius = NULL, ...) {
  coords <- sf::st_coordinates(stations)
  df <- data.frame(x = coords[, 1], y = coords[, 2], z = stations[[value_col]])
  v <- terra::vect(df, geom = c("x", "y"), crs = terra::crs(template))

  if (is.null(radius)) {
    e <- terra::ext(template)
    radius <- max(e[2] - e[1], e[4] - e[3])  # covers the full extent
  }

  out <- terra::interpNear(template[[1]], v, field = "z", radius = radius, ...)
  names(out) <- paste0(value_col, "_nn")
  out
}

#' Cross-validate interpolation methods (leave-one-out)
#'
#' Runs leave-one-out cross-validation for one or more interpolation methods
#' and reports RMSE / MAE, to help choose the best method for a given station
#' network and date (mirrors the diagnostic approach used in
#' \code{R.IBE.airqino}).
#'
#' @inheritParams interpolate_airtemp
#' @param methods Character vector of methods to compare.
#' @return A data.frame with one row per method: \code{method}, \code{rmse}, \code{mae}.
#' @export
cv_interpolate_airtemp <- function(stations, value_col, template,
                                    methods = c("idw", "kriging", "nn")) {
  n <- nrow(stations)
  results <- lapply(methods, function(m) {
    preds <- numeric(n)
    for (i in seq_len(n)) {
      train <- stations[-i, ]
      test <- stations[i, ]
      r <- tryCatch(
        interpolate_airtemp(train, value_col, template, method = m),
        error = function(e) NULL
      )
      if (is.null(r)) { preds[i] <- NA; next }
      preds[i] <- terra::extract(r, sf::st_transform(test, terra::crs(template)))[, 2]
    }
    obs <- stations[[value_col]]
    err <- obs - preds
    data.frame(method = m,
               rmse = sqrt(mean(err^2, na.rm = TRUE)),
               mae = mean(abs(err), na.rm = TRUE))
  })
  do.call(rbind, results)
}
