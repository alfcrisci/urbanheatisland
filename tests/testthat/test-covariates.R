test_that("compute_ndvi e compute_ndbi calcolano correttamente", {
  skip_if_not_installed("terra")
  nir <- terra::rast(nrows = 2, ncols = 2, vals = c(0.5, 0.4, 0.3, 0.6))
  red <- terra::rast(nrows = 2, ncols = 2, vals = c(0.1, 0.2, 0.25, 0.1))
  swir <- terra::rast(nrows = 2, ncols = 2, vals = c(0.3, 0.5, 0.2, 0.15))

  ndvi <- compute_ndvi(nir, red)
  expected_ndvi <- (c(0.5, 0.4, 0.3, 0.6) - c(0.1, 0.2, 0.25, 0.1)) /
                    (c(0.5, 0.4, 0.3, 0.6) + c(0.1, 0.2, 0.25, 0.1))
  expect_equal(terra::values(ndvi)[, 1], expected_ndvi, tolerance = 1e-6)
  expect_equal(names(ndvi), "NDVI")

  ndbi <- compute_ndbi(swir, nir)
  expected_ndbi <- (c(0.3, 0.5, 0.2, 0.15) - c(0.5, 0.4, 0.3, 0.6)) /
                    (c(0.3, 0.5, 0.2, 0.15) + c(0.5, 0.4, 0.3, 0.6))
  expect_equal(terra::values(ndbi)[, 1], expected_ndbi, tolerance = 1e-6)
})

test_that("define_rural_reference combina correttamente le soglie NDVI e impervious", {
  skip_if_not_installed("terra")
  ndvi <- terra::rast(nrows = 2, ncols = 2, vals = c(0.6, 0.2, 0.5, 0.1))       # rurale: >=0.4
  imperv <- terra::rast(nrows = 2, ncols = 2, vals = c(0.05, 0.05, 0.5, 0.05))  # rurale: <=0.1

  mask_ndvi_only <- define_rural_reference(ndvi = ndvi, ndvi_min = 0.4)
  expect_equal(terra::values(mask_ndvi_only)[, 1], c(1, 0, 1, 0))

  mask_combined <- define_rural_reference(ndvi = ndvi, ndvi_min = 0.4,
                                           impervious = imperv, impervious_max = 0.1)
  # solo la cella 1 soddisfa entrambe le condizioni (cella 3 ha NDVI ok ma impervious alto)
  expect_equal(terra::values(mask_combined)[, 1], c(1, 0, 0, 0))
})

test_that("define_rural_reference richiede almeno una covariata", {
  expect_error(define_rural_reference(), "almeno una covariata")
})

test_that("compute_uhi accetta una maschera SpatRaster come rural_ref", {
  skip_if_not_installed("terra")
  temp <- terra::rast(nrows = 2, ncols = 2, vals = c(30, 32, 20, 22))
  rural_mask <- terra::rast(nrows = 2, ncols = 2, vals = c(0, 0, 1, 1))  # celle 3-4 rurali

  uhi <- compute_uhi(temp, rural_mask)
  # baseline = media(20,22) = 21
  expect_equal(terra::values(uhi)[, 1], c(30, 32, 20, 22) - 21, tolerance = 1e-6)
})

test_that("model_lst_drivers stima un modello lineare coerente e ricostruisce i gap", {
  skip_if_not_installed("terra")
  set.seed(5)
  # costruisce una LST fortemente correlata a NDVI (relazione nota) con un buco
  ndvi_vals <- runif(100, 0, 1)
  lst_vals <- 35 - 10 * ndvi_vals + rnorm(100, sd = 0.3)
  lst_vals[c(10, 50)] <- NA  # simula pixel mascherati da nuvole

  ndvi_r <- terra::rast(nrows = 10, ncols = 10, vals = ndvi_vals)
  lst_r  <- terra::rast(nrows = 10, ncols = 10, vals = lst_vals)
  names(lst_r) <- "lst"

  res <- model_lst_drivers(lst_r, covariates = list(ndvi = ndvi_r), method = "lm")

  expect_true(res$r_squared > 0.8)  # relazione forte per costruzione
  expect_s3_class(res$model, "lm")
  expect_s4_class(res$predicted, "SpatRaster")
  expect_s4_class(res$gap_filled, "SpatRaster")

  # i pixel precedentemente NA ora devono avere un valore stimato
  gap_vals <- terra::values(res$gap_filled)[, 1]
  expect_false(any(is.na(gap_vals[c(10, 50)])))
})

test_that("model_lst_drivers segnala covariate non nominate", {
  skip_if_not_installed("terra")
  lst_r <- terra::rast(nrows = 5, ncols = 5, vals = runif(25, 20, 30))
  ndvi_r <- terra::rast(nrows = 5, ncols = 5, vals = runif(25))
  expect_error(model_lst_drivers(lst_r, covariates = list(ndvi_r)), "nominata")
})
