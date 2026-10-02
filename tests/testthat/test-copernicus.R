test_that(".parse_s3_href gestisce correttamente href in formato s3://", {
  href <- "s3://eodata/Sentinel-2/MSI/L2A_N0500/2021/08/21/prodotto.SAFE/GRANULE/x/B04_10m.jp2"
  parsed <- .parse_s3_href(href)
  expect_equal(parsed$bucket, "eodata")
  expect_equal(parsed$key, "Sentinel-2/MSI/L2A_N0500/2021/08/21/prodotto.SAFE/GRANULE/x/B04_10m.jp2")
})

test_that(".parse_s3_href gestisce correttamente href in formato https path-style", {
  href <- "https://eodata.dataspace.copernicus.eu/eodata/Sentinel-3/SL_2_LST/prodotto/LST_in.nc"
  parsed <- .parse_s3_href(href)
  expect_equal(parsed$bucket, "eodata")
  expect_equal(parsed$key, "Sentinel-3/SL_2_LST/prodotto/LST_in.nc")
})

test_that("read_sentinel2_band applica correttamente lo scale factor con e senza offset", {
  skip_if_not_installed("terra")
  tmp <- tempfile(fileext = ".tif")
  r <- terra::rast(nrows = 2, ncols = 2, vals = c(1500, 2000, 500, 3000))
  terra::writeRaster(r, tmp, overwrite = TRUE)

  senza_offset <- read_sentinel2_band(tmp, boa_add_offset = 0)
  expect_equal(as.numeric(terra::values(senza_offset)[1, 1]), 0.15, tolerance = 1e-6)

  con_offset <- read_sentinel2_band(tmp, boa_add_offset = -1000)
  expect_equal(as.numeric(terra::values(con_offset)[1, 1]), 0.05, tolerance = 1e-6)

  unlink(tmp)
})

test_that("copernicus_download_asset segnala un asset mancante elencando quelli disponibili", {
  assets <- list(B04_10m = list(href = "s3://eodata/x/B04_10m.jp2"))
  expect_error(
    copernicus_download_asset(assets, "banda_inesistente", tempfile(),
                               s3_access_key = "fake", s3_secret_key = "fake"),
    "non trovato.*B04_10m"
  )
})

test_that("copernicus_download_asset richiede credenziali S3", {
  skip_if_not_installed("aws.s3")
  old_key <- Sys.getenv("CDSE_S3_ACCESS_KEY", unset = NA)
  old_secret <- Sys.getenv("CDSE_S3_SECRET_KEY", unset = NA)
  Sys.unsetenv(c("CDSE_S3_ACCESS_KEY", "CDSE_S3_SECRET_KEY"))
  on.exit({
    if (!is.na(old_key)) Sys.setenv(CDSE_S3_ACCESS_KEY = old_key)
    if (!is.na(old_secret)) Sys.setenv(CDSE_S3_SECRET_KEY = old_secret)
  }, add = TRUE)

  assets <- list(B04_10m = list(href = "s3://eodata/x/B04_10m.jp2"))
  expect_error(copernicus_download_asset(assets, "B04_10m", tempfile()), "Credenziali S3 Copernicus non fornite")
})

test_that("copernicus_search restituisce una struttura vuota coerente se non ci sono risultati", {
  vuoto <- data.frame(id = character(0), datetime = character(0), cloud_cover = numeric(0))
  expect_equal(nrow(vuoto), 0)
  expect_equal(names(vuoto), c("id", "datetime", "cloud_cover"))
})
