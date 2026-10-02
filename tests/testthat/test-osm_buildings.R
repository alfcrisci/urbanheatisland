test_that("rasterize_buildings calcola coverage_fraction correttamente", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  template <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "EPSG:32632")

  # edificio interamente dentro la cella (0,0)-(10,10): area 36 su 100 -> frazione 0.36
  b1 <- sf::st_polygon(list(rbind(c(2, 2), c(8, 2), c(8, 8), c(2, 8), c(2, 2))))
  buildings <- sf::st_sf(id = 1, geometry = sf::st_sfc(b1, crs = 32632))

  cov <- rasterize_buildings(buildings, template, metric = "coverage_fraction")
  val <- as.numeric(terra::extract(cov, cbind(5, 5))[1, 1])
  expect_equal(val, 0.36, tolerance = 0.01)
  expect_true(all(as.numeric(terra::values(cov)) >= 0 & as.numeric(terra::values(cov)) <= 1))
})

test_that("rasterize_buildings supporta le metriche binary e count", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  template <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 40, ymin = 0, ymax = 40, crs = "EPSG:32632")
  b1 <- sf::st_polygon(list(rbind(c(2, 2), c(8, 2), c(8, 8), c(2, 8), c(2, 2))))
  b2 <- sf::st_polygon(list(rbind(c(15, 15), c(25, 15), c(25, 25), c(15, 25), c(15, 15))))
  buildings <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(b1, b2, crs = 32632))

  bin <- rasterize_buildings(buildings, template, metric = "binary")
  expect_true(all(as.numeric(terra::values(bin)) %in% c(0, 1)))
  expect_true(sum(as.numeric(terra::values(bin))) > 0)

  cnt <- rasterize_buildings(buildings, template, metric = "count")
  expect_equal(sum(as.numeric(terra::values(cnt)), na.rm = TRUE), 2)
})

test_that("rasterize_buildings gestisce correttamente l'assenza di edifici", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  template <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 30, ymin = 0, ymax = 30, crs = "EPSG:32632")
  b1 <- sf::st_polygon(list(rbind(c(2, 2), c(8, 2), c(8, 8), c(2, 8), c(2, 2))))
  buildings <- sf::st_sf(id = 1, geometry = sf::st_sfc(b1, crs = 32632))
  empty <- buildings[0, ]

  out <- rasterize_buildings(empty, template, metric = "coverage_fraction")
  expect_true(all(as.numeric(terra::values(out)) == 0))
  expect_equal(terra::ncell(out), terra::ncell(template))
})

test_that("get_osm_buildings segnala se osmdata non e' installato", {
  # simuliamo l'assenza del pacchetto sostituendo temporaneamente requireNamespace
  testthat::local_mocked_bindings(
    requireNamespace = function(pkg, ...) if (pkg == "osmdata") FALSE else TRUE,
    .package = "base"
  )
  aoi <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(11.25, 43.77)), crs = 4326)) |>
    sf::st_buffer(0.01)
  expect_error(get_osm_buildings(aoi), "osmdata")
})
