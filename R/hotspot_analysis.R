#' Getis-Ord Gi* hotspot analysis on a raster
#'
#' Computes the local Getis-Ord Gi* statistic for every non-NA cell of a
#' raster, identifying statistically significant clusters of high values
#' (hotspots) and low values (coldspots).
#'
#' For rasters with fewer than \code{max_cells_direct} valid cells, an exact
#' \code{spdep}-based computation is used (converts to points, builds a
#' neighbour list, computes \code{spdep::localG}). For larger rasters, a
#' faster \code{terra::focal}-based approximation is used with a fixed-radius
#' window, which is accurate for regular grids but does not use exact
#' great-circle/projected distances at tile edges.
#'
#' @param r A single-layer \code{SpatRaster}.
#' @param d Neighbourhood distance (map units of \code{r}'s CRS). Default is
#'   \code{1.5 * } the raster's cell resolution (~1 ring of 8 neighbours).
#' @param style Weight style passed to \code{spdep::nb2listw} (exact mode
#'   only). Default \code{"B"} (binary), standard for Gi*.
#' @param max_cells_direct Maximum number of valid cells for which the exact
#'   \code{spdep} method is used. Above this, the faster focal-window
#'   approximation is used. Default \code{200000}.
#' @return A single-layer \code{SpatRaster} of Gi* z-scores, same grid as \code{r}.
#' @export
gi_star_raster <- function(r, d = NULL, style = "B", max_cells_direct = 200000) {
  stopifnot(inherits(r, "SpatRaster"), terra::nlyr(r) == 1)
  if (is.null(d)) d <- terra::res(r)[1] * 1.5

  n_valid <- terra::global(r, fun = "notNA")[1, 1]

  if (n_valid <= max_cells_direct) {
    .gi_star_exact(r, d, style)
  } else {
    .gi_star_focal(r, d)
  }
}

#' @keywords internal
.gi_star_exact <- function(r, d, style) {
  pts <- terra::as.points(r, na.rm = TRUE)
  df <- terra::as.data.frame(pts, geom = "XY")
  val_col <- names(r)[1]
  coords <- as.matrix(df[, c("x", "y")])
  values_vec <- df[[val_col]]

  nb <- spdep::dnearneigh(coords, d1 = 0, d2 = d)
  nb <- spdep::include.self(nb)
  lw <- spdep::nb2listw(nb, style = style, zero.policy = TRUE)

  gi <- spdep::localG(values_vec, lw, zero.policy = TRUE)

  out <- r
  terra::values(out) <- NA_real_
  cell_idx <- terra::cellFromXY(r, coords)
  vals <- rep(NA_real_, terra::ncell(out))
  vals[cell_idx] <- as.numeric(gi)
  terra::values(out) <- vals
  names(out) <- paste0(val_col, "_Gi")
  out
}

#' @keywords internal
.gi_star_focal <- function(r, d) {
  cell_res <- terra::res(r)[1]
  radius_cells <- max(1, round(d / cell_res))
  w_size <- 2 * radius_cells + 1

  w <- matrix(1, nrow = w_size, ncol = w_size)
  # circular window mask
  center <- radius_cells + 1
  for (i in seq_len(w_size)) for (j in seq_len(w_size)) {
    if (sqrt((i - center)^2 + (j - center)^2) > radius_cells + 0.5) w[i, j] <- 0
  }
  n_j <- sum(w)

  vals <- terra::values(r)[, 1]
  x_bar <- mean(vals, na.rm = TRUE)
  s <- stats::sd(vals, na.rm = TRUE)
  n <- sum(!is.na(vals))

  local_sum <- terra::focal(r, w = w, fun = "sum", na.rm = TRUE, fillvalue = NA)
  w_sum <- n_j
  w_sq_sum <- n_j

  numerator <- local_sum - x_bar * w_sum
  denom <- s * sqrt((n * w_sq_sum - w_sum^2) / (n - 1))
  gi <- numerator / denom

  names(gi) <- paste0(names(r)[1], "_Gi")
  terra::mask(gi, r)
}

#' Local Moran's I hotspot / cluster analysis on a raster
#'
#' Computes Anselin's Local Moran's I for every non-NA cell, classifying
#' each cell into High-High, Low-Low (spatial clusters), High-Low, Low-High
#' (spatial outliers), or Not significant.
#'
#' @param r A single-layer \code{SpatRaster}.
#' @param d Neighbourhood distance (map units). Default \code{1.5 *} cell
#'   resolution.
#' @param style Weight style. Default \code{"W"} (row-standardised), standard
#'   for Local Moran's I.
#' @param alpha Significance threshold for cluster classification. Default \code{0.05}.
#' @return A \code{SpatRaster} stack with layers \code{Ii} (Local Moran's I),
#'   \code{p_value}, and \code{cluster} (categorical: 1=HH, 2=LL, 3=HL, 4=LH,
#'   0=not significant).
#' @export
local_moran_raster <- function(r, d = NULL, style = "W", alpha = 0.05) {
  stopifnot(inherits(r, "SpatRaster"), terra::nlyr(r) == 1)
  if (is.null(d)) d <- terra::res(r)[1] * 1.5

  pts <- terra::as.points(r, na.rm = TRUE)
  df <- terra::as.data.frame(pts, geom = "XY")
  val_col <- names(r)[1]
  coords <- as.matrix(df[, c("x", "y")])
  values_vec <- df[[val_col]]

  nb <- spdep::dnearneigh(coords, d1 = 0, d2 = d)
  lw <- spdep::nb2listw(nb, style = style, zero.policy = TRUE)

  lm_res <- spdep::localmoran(values_vec, lw, zero.policy = TRUE)

  z <- scale(values_vec)[, 1]
  lag_z <- spdep::lag.listw(lw, z, zero.policy = TRUE)

  cluster <- integer(length(values_vec))
  sig <- lm_res[, "Pr(z != E(Ii))"] < alpha
  cluster[sig & z > 0 & lag_z > 0] <- 1L  # HH
  cluster[sig & z < 0 & lag_z < 0] <- 2L  # LL
  cluster[sig & z > 0 & lag_z < 0] <- 3L  # HL
  cluster[sig & z < 0 & lag_z > 0] <- 4L  # LH

  cell_idx <- terra::cellFromXY(r, coords)

  make_layer <- function(values) {
    out <- r
    vals <- rep(NA_real_, terra::ncell(out))
    vals[cell_idx] <- values
    terra::values(out) <- vals
    out
  }

  ii_rast <- make_layer(lm_res[, "Ii"])
  p_rast <- make_layer(lm_res[, "Pr(z != E(Ii))"])
  cl_rast <- make_layer(cluster)

  names(ii_rast) <- "Ii"
  names(p_rast) <- "p_value"
  names(cl_rast) <- "cluster"
  levels(cl_rast) <- data.frame(
    id = 0:4,
    cluster = c("Not significant", "High-High", "Low-Low", "High-Low", "Low-High")
  )

  c(ii_rast, p_rast, cl_rast)
}

#' Classify a Gi* z-score raster into hot/cold spot confidence bands
#'
#' @param gi_rast A Gi* z-score \code{SpatRaster}, as returned by
#'   \code{\link{gi_star_raster}}.
#' @param apply_fdr Logical, apply a false-discovery-rate correction to
#'   two-sided p-values before classification (recommended for large numbers
#'   of cells / multiple testing). Default \code{TRUE}.
#' @return A categorical \code{SpatRaster} with 7 classes: Cold 99\%/95\%/90\%,
#'   Not significant, Hot 90\%/95\%/99\%.
#' @export
classify_hotspots <- function(gi_rast, apply_fdr = TRUE) {
  vals <- terra::values(gi_rast)[, 1]
  p_vals <- 2 * (1 - stats::pnorm(abs(vals)))

  if (apply_fdr) {
    valid <- !is.na(p_vals)
    p_adj <- rep(NA_real_, length(p_vals))
    p_adj[valid] <- stats::p.adjust(p_vals[valid], method = "fdr")
    p_vals <- p_adj
  }

  cls <- ifelse(is.na(vals), NA_integer_,
    ifelse(vals > 0 & p_vals < 0.01, 6L,
    ifelse(vals > 0 & p_vals < 0.05, 5L,
    ifelse(vals > 0 & p_vals < 0.10, 4L,
    ifelse(vals < 0 & p_vals < 0.01, 0L,
    ifelse(vals < 0 & p_vals < 0.05, 1L,
    ifelse(vals < 0 & p_vals < 0.10, 2L, 3L)))))))

  out <- gi_rast
  terra::values(out) <- cls
  levels(out) <- data.frame(
    id = 0:6,
    class = c("Cold 99%", "Cold 95%", "Cold 90%", "Not significant",
              "Hot 90%", "Hot 95%", "Hot 99%")
  )
  names(out) <- "hotspot_class"
  out
}
