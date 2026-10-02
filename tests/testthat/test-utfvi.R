test_that("compute_utfvi applica correttamente la formula in Kelvin", {
  skip_if_not_installed("terra")
  r <- terra::rast(nrows = 2, ncols = 2, vals = c(28, 30, 32, 34))
  names(r) <- "LST"

  utfvi <- compute_utfvi(r, input_unit = "celsius")
  mean_k <- mean(c(28, 30, 32, 34)) + 273.15
  expected <- ((c(28, 30, 32, 34) + 273.15) - mean_k) / mean_k

  expect_equal(as.numeric(terra::values(utfvi)[, 1]), expected, tolerance = 1e-8)
})

test_that("compute_utfvi in Celsius diretto (senza conversione) da' un risultato diverso", {
  skip_if_not_installed("terra")
  r <- terra::rast(nrows = 2, ncols = 2, vals = c(28, 30, 32, 34))
  utfvi_corretto <- compute_utfvi(r, input_unit = "celsius")
  utfvi_kelvin_falso <- compute_utfvi(r, input_unit = "kelvin")  # trattato erroneamente come se fosse gia' Kelvin
  expect_false(isTRUE(all.equal(
    as.numeric(terra::values(utfvi_corretto)[, 1]),
    as.numeric(terra::values(utfvi_kelvin_falso)[, 1])
  )))
})

test_that("compute_utfvi con reference='per_layer' azzera la media di ogni layer", {
  skip_if_not_installed("terra")
  r1 <- terra::rast(nrows = 3, ncols = 3, vals = c(28, 30, 32, 29, 31, 33, 27, 30, 32))
  r2 <- r1 + 5
  stack <- c(r1, r2)
  names(stack) <- c("a", "b")

  utfvi <- compute_utfvi(stack, reference = "per_layer")
  expect_equal(terra::nlyr(utfvi), 2)
  for (i in 1:2) {
    expect_equal(terra::global(utfvi[[i]], "mean", na.rm = TRUE)[1, 1], 0, tolerance = 1e-6)
  }
})

test_that("compute_utfvi con reference numerico fisso usa quella baseline", {
  skip_if_not_installed("terra")
  r <- terra::rast(nrows = 2, ncols = 2, vals = c(28, 30, 32, 34))
  utfvi <- compute_utfvi(r, reference = 30)
  mean_k <- 30 + 273.15
  expected <- ((c(28, 30, 32, 34) + 273.15) - mean_k) / mean_k
  expect_equal(as.numeric(terra::values(utfvi)[, 1]), expected, tolerance = 1e-8)
})

test_that("classify_utfvi assegna le classi standard correttamente", {
  skip_if_not_installed("terra")
  # valori scelti per cadere chiaramente in ciascuna delle 6 classi
  vals <- c(-0.01, 0.002, 0.007, 0.012, 0.017, 0.025)
  r <- terra::rast(nrows = 2, ncols = 3, vals = vals)
  names(r) <- "UTFVI"

  cls <- classify_utfvi(r)
  lv <- terra::levels(cls)[[1]]
  expect_equal(lv$classe, c("Eccellente", "Buona", "Normale", "Scarsa", "Peggiore", "Pessima"))

  freq_tab <- terra::freq(cls)
  expect_equal(sort(freq_tab$value),
               sort(c("Eccellente", "Buona", "Normale", "Scarsa", "Peggiore", "Pessima")))
})

test_that("classify_utfvi segnala breaks/labels incoerenti", {
  skip_if_not_installed("terra")
  r <- terra::rast(nrows = 2, ncols = 2, vals = c(-0.01, 0.01, 0.02, 0.03))
  expect_error(classify_utfvi(r, breaks = c(-Inf, 0, Inf), labels = c("a", "b", "c")),
               "elemento in piu")
})

test_that("summarize_utfvi produce percentuali che sommano a 100 (singolo layer)", {
  skip_if_not_installed("terra")
  vals <- c(-0.01, 0.002, 0.007, 0.012, 0.017, 0.025, -0.01, 0.002, 0.007)
  r <- terra::rast(nrows = 3, ncols = 3, vals = vals)
  cls <- classify_utfvi(r)

  summ <- summarize_utfvi(cls)
  expect_true(all(c("classe", "n_celle", "area_km2", "percentuale") %in% names(summ)))
  expect_equal(sum(summ$percentuale), 100, tolerance = 1e-6)
  expect_equal(sum(summ$n_celle), 9)
})

test_that("summarize_utfvi gestisce correttamente stack multi-layer", {
  skip_if_not_installed("terra")
  vals <- c(-0.01, 0.002, 0.007, 0.012, 0.017, 0.025, -0.01, 0.002, 0.007)
  r1 <- terra::rast(nrows = 3, ncols = 3, vals = vals)
  r2 <- terra::rast(nrows = 3, ncols = 3, vals = rev(vals))
  stack <- c(r1, r2)
  names(stack) <- c("data1", "data2")

  cls_stack <- classify_utfvi(stack)
  expect_equal(terra::nlyr(cls_stack), 2)

  summ <- summarize_utfvi(cls_stack)
  expect_true("layer" %in% names(summ))
  expect_equal(length(unique(summ$layer)), 2)
  # ogni layer deve sommare a 100%
  by_layer <- tapply(summ$percentuale, summ$layer, sum)
  expect_true(all(abs(by_layer - 100) < 1e-6))
})
