test_that("gi_star_raster detects a synthetic hotspot", {
  skip_if_not_installed("terra")
  skip_if_not_installed("spdep")

  set.seed(1)
  r <- terra::rast(nrows = 30, ncols = 30, xmin = 0, xmax = 30, ymin = 0, ymax = 30)
  vals <- matrix(stats::rnorm(900, mean = 20, sd = 1), nrow = 30)
  # inject a clear hotspot block
  vals[10:15, 10:15] <- vals[10:15, 10:15] + 15
  terra::values(r) <- as.vector(t(vals))
  names(r) <- "temp"

  gi <- gi_star_raster(r, d = 1.8, max_cells_direct = 1e6)

  expect_s4_class(gi, "SpatRaster")
  expect_equal(terra::ncell(gi), terra::ncell(r))

  # cell near hotspot centre should have a strongly positive Gi*
  center_xy <- terra::xyFromCell(r, terra::cellFromRowCol(r, 12, 12))
  gi_center <- terra::extract(gi, center_xy)[1, 1]
  expect_true(gi_center > 1.5)
})

test_that("classify_hotspots returns expected categorical levels", {
  skip_if_not_installed("terra")

  r <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 5)
  terra::values(r) <- c(-3, -2, -1.7, 0, 1.7, 2, 3, rep(0, 18))
  names(r) <- "z"

  cls <- classify_hotspots(r, apply_fdr = FALSE)
  lv <- terra::levels(cls)[[1]]

  expect_true(all(c("Cold 99%", "Hot 99%", "Not significant") %in% lv$class))
})

test_that("decode_modis_qc keeps only requested classes", {
  skip_if_not_installed("terra")

  qc <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2)
  terra::values(qc) <- c(0, 1, 2, 3)  # good, marginal, cloud, other

  mask_good_only <- decode_modis_qc(qc, keep_classes = "good")
  expect_equal(terra::values(mask_good_only)[, 1], c(1, 0, 0, 0))

  mask_good_marg <- decode_modis_qc(qc, keep_classes = c("good", "marginal"))
  expect_equal(terra::values(mask_good_marg)[, 1], c(1, 1, 0, 0))
})
