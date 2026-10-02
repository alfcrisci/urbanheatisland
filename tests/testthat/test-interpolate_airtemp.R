test_that("interpolate_airtemp IDW and NN produce full-grid rasters", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")
  skip_if_not_installed("gstat")

  set.seed(42)
  coords <- cbind(runif(15, 0, 10), runif(15, 0, 10))
  stations <- sf::st_as_sf(
    data.frame(coords, temp = rnorm(15, 22, 2)),
    coords = c("X1", "X2"), crs = 4326
  )

  template <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 10, ymin = 0, ymax = 10, crs = "EPSG:4326")

  idw_out <- interpolate_airtemp(stations, "temp", template, method = "idw")
  expect_s4_class(idw_out, "SpatRaster")
  expect_equal(terra::ncell(idw_out), terra::ncell(template))
  expect_true(sum(!is.na(terra::values(idw_out))) > 0)

  nn_out <- interpolate_airtemp(stations, "temp", template, method = "nn")
  expect_s4_class(nn_out, "SpatRaster")
})

test_that("interpolate_airtemp rejects invalid inputs", {
  skip_if_not_installed("sf")
  skip_if_not_installed("terra")

  df <- data.frame(x = 1:3, y = 1:3, temp = c(20, 21, 22))
  not_sf <- df
  template <- terra::rast(nrows = 5, ncols = 5)

  expect_error(interpolate_airtemp(not_sf, "temp", template, method = "idw"),
               "must be an sf object")

  stations <- sf::st_as_sf(df, coords = c("x", "y"), crs = 4326)
  expect_error(interpolate_airtemp(stations, "nonexistent_col", template, method = "idw"),
               "not found")
})
