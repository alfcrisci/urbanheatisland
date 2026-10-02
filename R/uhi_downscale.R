#' Compute urban heat island intensity from LST
#'
#' Computes cell-wise UHI intensity as the difference between an LST (or
#' interpolated air temperature) raster and a reference "rural" value,
#' either a fixed reference raster/mask, or the mean over a rural reference
#' zone.
#'
#' @param temp_rast A \code{SpatRaster} of surface or air temperature.
#' @param rural_ref Riferimento rurale: un poligono \code{sf}/\code{SpatVector}
#'   (usa la media al suo interno), un \code{SpatRaster} binario (1 = rurale,
#'   0/NA = escluso — es. da \code{\link{define_rural_reference}}), oppure un
#'   valore numerico fisso di temperatura di riferimento.
#' @return A \code{SpatRaster} of UHI intensity (same units as \code{temp_rast}).
#' @export
compute_uhi <- function(temp_rast, rural_ref) {
  if (is.numeric(rural_ref)) {
    baseline <- rural_ref
  } else if (inherits(rural_ref, "SpatRaster")) {
    mask_aligned <- if (!terra::compareGeom(temp_rast, rural_ref, stopOnError = FALSE)) {
      align_rasters(temp_rast, rural_ref)
    } else rural_ref
    masked <- terra::mask(temp_rast, mask_aligned, maskvalues = c(0, NA), inverse = FALSE)
    baseline <- terra::global(masked, "mean", na.rm = TRUE)[1, 1]
  } else {
    v <- if (inherits(rural_ref, "sf")) terra::vect(rural_ref) else rural_ref
    baseline <- terra::extract(temp_rast, v, fun = mean, na.rm = TRUE)[, 2]
    baseline <- mean(baseline, na.rm = TRUE)
  }
  uhi <- temp_rast - baseline
  names(uhi) <- paste0(names(temp_rast), "_UHI")
  uhi
}

#' Align two rasters to a common grid
#'
#' Resamples \code{r2} to match the extent, resolution, and CRS of \code{r1}.
#'
#' @param r1 Reference \code{SpatRaster}.
#' @param r2 \code{SpatRaster} to align to \code{r1}.
#' @param method Resampling method, passed to \code{terra::resample}.
#'   Default \code{"bilinear"}.
#' @return \code{r2} resampled onto \code{r1}'s grid.
#' @export
align_rasters <- function(r1, r2, method = "bilinear") {
  if (!terra::same.crs(r1, r2)) r2 <- terra::project(r2, terra::crs(r1))
  terra::resample(r2, r1, method = method)
}

#' Downscale interpolated air temperature using LST as a covariate
#'
#' Produces a fine-resolution air temperature surface by combining a coarser
#' interpolated air-temperature raster with a finer-resolution LST raster,
#' using either a simple linear regression, regression-kriging of the
#' residuals, or random forest.
#'
#' @param airtemp_rast Coarse interpolated air temperature \code{SpatRaster}
#'   (e.g. from \code{\link{interpolate_airtemp}}).
#' @param lst_rast Fine-resolution LST \code{SpatRaster}, used as the
#'   downscaling covariate.
#' @param stations Optional \code{sf} station points with observed air
#'   temperature, required for \code{method = "rk"} (regression-kriging of
#'   residuals) to get point-level residuals.
#' @param value_col Column name of observed temperature in \code{stations}
#'   (required for \code{method = "rk"}).
#' @param method One of \code{"regression"} (global linear model of
#'   air temp ~ LST, applied pixel-wise), \code{"rk"} (regression-kriging:
#'   linear trend + kriged residuals), or \code{"rf"} (random forest).
#' @return A \code{SpatRaster} at the resolution of \code{lst_rast}.
#' @export
downscale_airtemp <- function(airtemp_rast, lst_rast, stations = NULL,
                               value_col = NULL,
                               method = c("regression", "rk", "rf")) {
  method <- match.arg(method)
  lst_al <- if (!terra::compareGeom(airtemp_rast, lst_rast, stopOnError = FALSE)) {
    align_rasters(lst_rast, airtemp_rast)
  } else airtemp_rast

  airtemp_on_lst <- align_rasters(lst_rast, airtemp_rast)
  stack <- c(airtemp_on_lst, lst_rast)
  names(stack) <- c("airtemp", "lst")
  df <- terra::as.data.frame(stack, xy = TRUE, na.rm = TRUE)

  if (method == "regression") {
    fit <- stats::lm(airtemp ~ lst, data = df)
    lst_named <- lst_rast
    names(lst_named) <- "lst"
    pred <- terra::predict(lst_named, fit)
    names(pred) <- "airtemp_downscaled"
    return(pred)
  }

  if (method == "rk") {
    if (is.null(stations) || is.null(value_col)) {
      stop("`stations` and `value_col` are required for method = 'rk'.", call. = FALSE)
    }
    fit <- stats::lm(airtemp ~ lst, data = df)
    lst_named <- lst_rast
    names(lst_named) <- "lst"
    trend_pred <- terra::predict(lst_named, fit)

    st_lst <- terra::extract(lst_rast, sf::st_transform(stations, terra::crs(lst_rast)))[, 2]
    trend_at_station <- fit$coefficients[1] + fit$coefficients[2] * st_lst
    resid_val <- stations[[value_col]] - trend_at_station

    resid_sf <- sf::st_sf(resid = resid_val, geometry = sf::st_geometry(stations))
    resid_rast <- interpolate_airtemp(resid_sf, "resid", lst_rast, method = "kriging")

    out <- trend_pred + resid_rast
    names(out) <- "airtemp_downscaled_rk"
    return(out)
  }

  if (method == "rf") {
    if (!requireNamespace("randomForest", quietly = TRUE)) {
      stop("Package 'randomForest' required for method = 'rf'.", call. = FALSE)
    }
    fit <- randomForest::randomForest(airtemp ~ lst, data = df)
    lst_named <- lst_rast
    names(lst_named) <- "lst"
    pred <- terra::predict(lst_named, fit)
    names(pred) <- "airtemp_downscaled_rf"
    return(pred)
  }
}
