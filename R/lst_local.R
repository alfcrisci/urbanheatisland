#' Read a local LST raster file
#'
#' Reads a GeoTIFF (e.g. MODIS MOD11A1/A2, or any pre-downloaded LST product)
#' and optionally applies the MODIS scale factor and Kelvin-to-Celsius
#' conversion, and/or QC-based masking.
#'
#' @param path Path to LST GeoTIFF.
#' @param qc_path Optional path to matching QC GeoTIFF (same product family).
#' @param scale Multiplicative scale factor applied to raw digital numbers.
#'   Default \code{0.02} (standard MODIS LST scale factor). Set to \code{1}
#'   if the file already stores physical values.
#' @param kelvin_to_celsius Logical, subtract 273.15 after scaling. Default \code{TRUE}.
#' @param keep_qc_classes Character vector of QC classes to keep when
#'   \code{qc_path} is supplied. One or more of \code{"good"}, \code{"marginal"},
#'   \code{"cloud"}, \code{"other"}. Default \code{c("good", "marginal")}.
#' @return A \code{SpatRaster} in degrees Celsius (or raw units if
#'   \code{scale = 1, kelvin_to_celsius = FALSE}).
#' @export
read_lst_local <- function(path, qc_path = NULL, scale = 0.02,
                            kelvin_to_celsius = TRUE,
                            keep_qc_classes = c("good", "marginal")) {
  r <- terra::rast(path)
  r <- r * scale
  if (kelvin_to_celsius) r <- r - 273.15

  if (!is.null(qc_path)) {
    qc <- terra::rast(qc_path)
    good_mask <- decode_modis_qc(qc, keep_classes = keep_qc_classes)
    r <- terra::mask(r, good_mask, maskvalues = 0)
  }
  r
}

#' Decode MODIS LST QC bit flags
#'
#' Decodes the 8-bit QC layer accompanying MOD11A1/A2 (and MYD equivalents)
#' LST products, following the official MODIS LST QC bit table (bits 0-1:
#' mandatory QA flag).
#'
#' @param qc A \code{SpatRaster} of raw QC byte values.
#' @param keep_classes Character vector of classes to mark as \code{1}
#'   (keep); everything else becomes \code{0}. One or more of \code{"good"}
#'   (bits 0-1 = 00, "LST produced, good quality"), \code{"marginal"} (01,
#'   "LST produced, other quality"), \code{"cloud"} (10, "not produced due to
#'   cloud"), \code{"other"} (11, "not produced, other reason").
#' @return A binary \code{SpatRaster} mask (1 = keep, 0 = discard, NA preserved).
#' @export
decode_modis_qc <- function(qc, keep_classes = c("good", "marginal")) {
  class_codes <- c(good = 0, marginal = 1, cloud = 2, other = 3)
  keep_codes <- class_codes[keep_classes]
  if (anyNA(keep_codes)) stop("Unknown QC class in keep_classes.", call. = FALSE)

  # bits 0-1 = value %% 4
  mandatory_qa <- qc %% 4
  mask <- terra::app(mandatory_qa, function(x) as.numeric(x %in% keep_codes))
  mask
}

#' Mask an LST raster using a QC raster
#'
#' @param lst A LST \code{SpatRaster}.
#' @param qc A raw QC \code{SpatRaster}, same extent/resolution as \code{lst}.
#' @param keep_classes See \code{\link{decode_modis_qc}}.
#' @return Masked \code{SpatRaster}.
#' @export
mask_lst_by_qc <- function(lst, qc, keep_classes = c("good", "marginal")) {
  good_mask <- decode_modis_qc(qc, keep_classes = keep_classes)
  terra::mask(lst, good_mask, maskvalues = 0)
}

#' Mosaic multiple LST tiles/dates
#'
#' @param rast_list A list of \code{SpatRaster} objects, or a character
#'   vector of file paths, to mosaic together.
#' @param fun Function used to combine overlapping cells. Default \code{"mean"}.
#' @return A single mosaicked \code{SpatRaster}.
#' @export
mosaic_lst <- function(rast_list, fun = "mean") {
  if (is.character(rast_list)) rast_list <- lapply(rast_list, terra::rast)
  sprc <- terra::sprc(rast_list)
  terra::mosaic(sprc, fun = fun)
}

#' Composite a time series of LST rasters (e.g. cloud-gap filling)
#'
#' Reduces a multi-layer LST \code{SpatRaster} (one layer per date) to a
#' single composite using a per-cell summary function, useful for filling
#' cloud gaps across an 8-day or monthly window.
#'
#' @param lst_stack A multi-layer \code{SpatRaster}.
#' @param fun Summary function applied per cell across layers. Default
#'   \code{"mean"}. \code{"median"} is more robust to residual cloud
#'   contamination.
#' @param min_obs Minimum number of non-NA observations required to keep a
#'   cell; cells with fewer valid observations are set to \code{NA}.
#' @return A single-layer \code{SpatRaster} composite.
#' @export
composite_lst <- function(lst_stack, fun = "mean", min_obs = 1) {
  comp <- terra::app(lst_stack, fun, na.rm = TRUE)
  n_obs <- terra::app(lst_stack, function(x) sum(!is.na(x)))
  terra::mask(comp, n_obs, maskvalues = (0:(min_obs - 1)))
}
