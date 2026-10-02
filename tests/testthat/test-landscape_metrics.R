test_that("classify_landcover assegna correttamente le classi", {
  skip_if_not_installed("terra")
  ndvi <- terra::rast(nrows = 2, ncols = 2, vals = c(0.6, 0.1, 0.5, 0.05))
  edificato <- terra::rast(nrows = 2, ncols = 2, vals = c(0.02, 0.02, 0.3, 0.3))

  lc <- classify_landcover(ndvi, ndvi_veg_min = 0.3, edificato = edificato, edificato_built_min = 0.1)
  # cella1: NDVI alto, edificato basso -> vegetazione (1)
  # cella2: NDVI basso, edificato basso -> altro (0)
  # cella3: NDVI alto MA edificato alto -> edificato prevale (2)
  # cella4: NDVI basso, edificato alto -> edificato (2)
  expect_equal(as.numeric(terra::values(lc)[, 1]), c(1, 0, 2, 2))
})

test_that("classify_landcover richiede almeno una fonte per l'edificato", {
  skip_if_not_installed("terra")
  ndvi <- terra::rast(nrows = 2, ncols = 2, vals = c(0.5, 0.2, 0.4, 0.1))
  expect_error(classify_landcover(ndvi), "edificato.*ndbi|ndbi.*edificato")
})

test_that("compute_landscape_metrics richiede window_size dispari", {
  skip_if_not_installed("terra")
  skip_if_not_installed("landscapemetrics")
  lc <- terra::rast(nrows = 10, ncols = 10, vals = sample(0:2, 100, replace = TRUE))
  expect_error(compute_landscape_metrics(lc, window_size = 4), "dispari")
})

test_that("compute_landscape_metrics multi-classe produce un raster valido", {
  skip_if_not_installed("terra")
  skip_if_not_installed("landscapemetrics")
  set.seed(1)
  lc <- terra::rast(nrows = 15, ncols = 15, xmin = 0, xmax = 150, ymin = 0, ymax = 150, crs = "EPSG:32632")
  terra::values(lc) <- sample(0:2, 225, replace = TRUE, prob = c(0.3, 0.4, 0.3))
  lc <- terra::as.int(lc)

  shdi <- compute_landscape_metrics(lc, window_size = 5, what = "lsm_l_shdi")
  expect_s4_class(shdi, "SpatRaster")
  expect_equal(terra::ncell(shdi), terra::ncell(lc))
  expect_equal(names(shdi), "lsm_l_shdi")
})

test_that("compute_landscape_metrics per classe produce layer nominati correttamente", {
  skip_if_not_installed("terra")
  skip_if_not_installed("landscapemetrics")
  set.seed(2)
  lc <- terra::rast(nrows = 15, ncols = 15, xmin = 0, xmax = 150, ymin = 0, ymax = 150, crs = "EPSG:32632")
  terra::values(lc) <- sample(0:2, 225, replace = TRUE, prob = c(0.3, 0.4, 0.3))
  lc <- terra::as.int(lc)

  out <- compute_landscape_metrics(lc, window_size = 5, what = c("lsm_l_ed", "lsm_l_pd"),
                                    classes_of_interest = c(1, 2),
                                    class_labels = c("vegetazione", "edificato"))
  expect_equal(names(out), c("vegetazione_ed", "vegetazione_pd", "edificato_ed", "edificato_pd"))
  expect_equal(terra::nlyr(out), 4)
  # densita' margini/patch non devono essere tutte zero (verifica della codifica
  # corretta dello sfondo come classe 0, non NA)
  expect_true(any(as.numeric(terra::values(out[["vegetazione_ed"]])) > 0, na.rm = TRUE))
})

test_that("compute_landscape_metrics segnala class_labels di lunghezza errata", {
  skip_if_not_installed("terra")
  skip_if_not_installed("landscapemetrics")
  lc <- terra::rast(nrows = 10, ncols = 10, vals = sample(0:2, 100, replace = TRUE))
  lc <- terra::as.int(lc)
  expect_error(
    compute_landscape_metrics(lc, classes_of_interest = c(1, 2), class_labels = "solo_uno"),
    "stessa lunghezza"
  )
})

test_that("get_landscape_covariates restituisce tutti i layer attesi", {
  skip_if_not_installed("terra")
  skip_if_not_installed("landscapemetrics")
  set.seed(4)
  n <- 15
  template <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n * 30, ymin = 0, ymax = n * 30, crs = "EPSG:32632")
  ndvi <- template
  terra::values(ndvi) <- runif(n * n, 0, 1)
  edificato <- template
  terra::values(edificato) <- runif(n * n, 0, 0.5)

  cov <- get_landscape_covariates(ndvi, edificato = edificato, window_size = 5)
  expect_equal(names(cov), c("vegetazione_ed", "vegetazione_pd", "edificato_ed", "edificato_pd", "shdi"))
})

test_that("le covariate di landscape ecology si integrano con model_lst_drivers", {
  skip_if_not_installed("terra")
  skip_if_not_installed("landscapemetrics")
  set.seed(6)
  n <- 15
  template <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n * 30, ymin = 0, ymax = n * 30, crs = "EPSG:32632")
  ndvi <- template; terra::values(ndvi) <- runif(n * n, 0, 1)
  edificato <- template; terra::values(edificato) <- runif(n * n, 0, 0.5)
  lst <- template; terra::values(lst) <- rnorm(n * n, 30, 2); names(lst) <- "lst"

  cov <- get_landscape_covariates(ndvi, edificato = edificato, window_size = 5)

  driver <- model_lst_drivers(lst, covariates = list(shdi = cov[["shdi"]], veg_ed = cov[["vegetazione_ed"]]),
                               method = "lm", predict_gaps = FALSE)
  expect_s3_class(driver$model, "lm")
  expect_true(is.numeric(driver$r_squared))
})
