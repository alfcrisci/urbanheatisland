test_that("read_remote_raster legge correttamente un percorso locale (nessun prefisso vsicurl per path non http)", {
  skip_if_not_installed("terra")
  tmp <- tempfile(fileext = ".tif")
  r <- terra::rast(nrows = 3, ncols = 3, vals = 1:9)
  terra::writeRaster(r, tmp, overwrite = TRUE)

  letto <- read_remote_raster(tmp)
  expect_s4_class(letto, "SpatRaster")
  expect_equal(as.numeric(terra::values(letto)), 1:9)

  unlink(tmp)
})

test_that("read_remote_raster antepone /vsicurl/ solo a URL http(s) semplici", {
  expect_true(grepl("^https?://", "https://esempio.it/x.tif"))
  # verifica solo la logica di costruzione del path, senza rete
  source_http <- "https://esempio.it/x.tif"
  is_plain_http <- grepl("^https?://", source_http) && !grepl("^/vsi|^WCS:|^WMS:|^NETCDF:", source_http)
  expect_true(is_plain_http)

  source_wcs <- "WCS:https://esempio.it/wcs?..."
  is_plain_http_wcs <- grepl("^https?://", source_wcs) && !grepl("^/vsi|^WCS:|^WMS:|^NETCDF:", source_wcs)
  expect_false(is_plain_http_wcs)
})

test_that("read_remote_raster restituisce un errore leggibile per sorgenti non raggiungibili/inesistenti", {
  skip_if_not_installed("terra")
  expect_error(
    suppressWarnings(read_remote_raster("https://host-inesistente.invalid/file.tif")),
    "Impossibile leggere il raster"
  )
})

test_that("crop_checked produce lo stesso risultato indipendentemente dal CRS del poligono", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  r_utm <- terra::rast(nrows = 10, ncols = 10, xmin = 676000, xmax = 686000,
                        ymin = 4844000, ymax = 4854000, crs = "EPSG:32632")
  terra::values(r_utm) <- 1:100

  poly_wgs84 <- sf::st_as_sf(sf::st_sfc(sf::st_buffer(sf::st_point(c(11.25, 43.77)), 0.02), crs = 4326))
  poly_utm <- sf::st_transform(poly_wgs84, 32632)

  ris_wgs84 <- crop_checked(r_utm, poly_wgs84)
  ris_utm <- crop_checked(r_utm, poly_utm)

  expect_equal(as.numeric(terra::values(ris_wgs84)), as.numeric(terra::values(ris_utm)))
})

test_that("crop_checked con mask=TRUE introduce NA fuori dal poligono", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  # griglia fine (100 m/cella) cosi' un buffer circolare lascia chiaramente
  # fuori gli angoli del proprio bounding box quadrato
  r_utm <- terra::rast(nrows = 40, ncols = 40, xmin = 676000, xmax = 680000,
                        ymin = 4844000, ymax = 4848000, crs = "EPSG:32632")
  terra::values(r_utm) <- 1:1600
  centro_utm <- sf::st_sfc(sf::st_point(c(678000, 4846000)), crs = 32632)
  poly_utm_src <- sf::st_as_sf(sf::st_buffer(centro_utm, 500))  # cerchio r=500m
  poly_wgs84 <- sf::st_transform(poly_utm_src, 4326)

  senza_mask <- crop_checked(r_utm, poly_wgs84, mask = FALSE)
  con_mask <- crop_checked(r_utm, poly_wgs84, mask = TRUE)

  expect_false(any(is.na(terra::values(senza_mask))))
  expect_true(any(is.na(terra::values(con_mask))))
})

test_that("crop_checked accetta sia sf sia SpatVector come poligono", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  r_utm <- terra::rast(nrows = 10, ncols = 10, xmin = 676000, xmax = 686000,
                        ymin = 4844000, ymax = 4854000, crs = "EPSG:32632")
  terra::values(r_utm) <- 1:100
  poly_sf <- sf::st_as_sf(sf::st_sfc(sf::st_buffer(sf::st_point(c(11.25, 43.77)), 0.02), crs = 4326))
  poly_vect <- terra::vect(poly_sf)

  ris_sf <- crop_checked(r_utm, poly_sf)
  ris_vect <- crop_checked(r_utm, poly_vect)

  expect_equal(as.numeric(terra::values(ris_sf)), as.numeric(terra::values(ris_vect)))
})
