test_that("decode_landsat_qa_pixel correctly identifies clear vs cloudy pixels", {
  skip_if_not_installed("terra")

  qa <- terra::rast(nrows = 2, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 2)

  # bit values: clear land (bit6=1) = 64
  #             cloud (bit3=1) = 8
  #             cloud shadow (bit4=1) = 16
  #             cirrus (bit2=1) = 4
  #             dilated cloud (bit1=1) = 2
  #             fill (bit0=1) = 1
  terra::values(qa) <- c(64, 8, 16, 4, 2, 1)

  mask <- decode_landsat_qa_pixel(qa)  # default excludes fill/dilated/cirrus/cloud/shadow
  vals <- terra::values(mask)[, 1]

  expect_equal(vals, c(1, 0, 0, 0, 0, 0))
})

test_that("decode_landsat_qa_pixel respects a custom exclude list", {
  skip_if_not_installed("terra")

  qa <- terra::rast(nrows = 1, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 1)
  terra::values(qa) <- c(64, 32)  # clear, snow

  # only exclude snow -> clear pixel kept, snow pixel discarded
  mask <- decode_landsat_qa_pixel(qa, exclude = "snow")
  expect_equal(terra::values(mask)[, 1], c(1, 0))
})

test_that("decode_landsat_qa_pixel rejects unknown flags", {
  skip_if_not_installed("terra")
  qa <- terra::rast(nrows = 1, ncols = 1)
  terra::values(qa) <- 64
  expect_error(decode_landsat_qa_pixel(qa, exclude = "not_a_flag"), "Unknown flag")
})

test_that("read_landsat_lst applies the correct Collection 2 scale factors", {
  skip_if_not_installed("terra")

  tmp <- tempfile(fileext = ".tif")
  r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2)
  # raw DN chosen so that DN * 0.00341802 + 149.0 = 300 K = 26.85 C
  dn <- (300 - 149.0) / 0.00341802
  terra::values(r) <- rep(dn, 4)
  terra::writeRaster(r, tmp, overwrite = TRUE)

  out <- read_landsat_lst(tmp)
  expect_equal(as.numeric(terra::values(out)[1, 1]), 300 - 273.15, tolerance = 1e-3)
  expect_equal(names(out), "LST_celsius")

  unlink(tmp)
})

test_that(".strip_scene_suffix removes only the trailing _SR/_ST collection suffix", {
  expect_equal(.strip_scene_suffix("LC08_L2SP_192030_20250715_02_T1_SR"), "LC08_L2SP_192030_20250715_02_T1")
  expect_equal(.strip_scene_suffix("LC08_L2SP_192030_20250715_02_T1_ST"), "LC08_L2SP_192030_20250715_02_T1")
  # non deve toccare id senza suffisso di collezione
  expect_equal(.strip_scene_suffix("LC08_L2SP_192030_20250715_02_T1"), "LC08_L2SP_192030_20250715_02_T1")
})

test_that(".merge_landsat_scenes combines SR and ST hrefs for the same scene", {
  sr_df <- data.frame(
    base_id = "LC08_L2SP_192030_20250715_02_T1",
    datetime = "2025-07-15T10:00:00Z", cloud_cover = 5.2,
    st_b10_href = NA_character_, qa_pixel_href = "https://sr/QA_PIXEL.tif",
    sr_b4_href = "https://sr/SR_B4.tif", sr_b5_href = "https://sr/SR_B5.tif",
    sr_b6_href = "https://sr/SR_B6.tif", stringsAsFactors = FALSE
  )
  st_df <- data.frame(
    base_id = "LC08_L2SP_192030_20250715_02_T1",
    datetime = "2025-07-15T10:00:00Z", cloud_cover = 5.2,
    st_b10_href = "https://st/ST_B10.tif", qa_pixel_href = "https://st/QA_PIXEL.tif",
    sr_b4_href = NA_character_, sr_b5_href = NA_character_, sr_b6_href = NA_character_,
    stringsAsFactors = FALSE
  )

  merged <- .merge_landsat_scenes(sr_df, st_df)

  expect_equal(nrow(merged), 1)
  expect_equal(merged$id, "LC08_L2SP_192030_20250715_02_T1")
  # il bug corretto: sia la banda termica (ST) sia le bande SR devono essere presenti
  expect_false(is.na(merged$st_b10_href))
  expect_false(is.na(merged$sr_b4_href))
  expect_false(is.na(merged$sr_b5_href))
  expect_false(is.na(merged$sr_b6_href))
  expect_equal(merged$st_b10_href, "https://st/ST_B10.tif")
  expect_equal(merged$sr_b4_href, "https://sr/SR_B4.tif")
})

test_that(".merge_landsat_scenes keeps scenes found in only one collection, with NA for the missing bands", {
  sr_df <- data.frame(
    base_id = "SCENE_A", datetime = "2025-07-01T10:00:00Z", cloud_cover = 3,
    st_b10_href = NA_character_, qa_pixel_href = "https://sr/QA_PIXEL.tif",
    sr_b4_href = "https://sr/SR_B4.tif", sr_b5_href = "https://sr/SR_B5.tif",
    sr_b6_href = "https://sr/SR_B6.tif", stringsAsFactors = FALSE
  )
  st_df <- data.frame(
    base_id = "SCENE_B", datetime = "2025-07-10T10:00:00Z", cloud_cover = 8,
    st_b10_href = "https://st/ST_B10.tif", qa_pixel_href = "https://st/QA_PIXEL.tif",
    sr_b4_href = NA_character_, sr_b5_href = NA_character_, sr_b6_href = NA_character_,
    stringsAsFactors = FALSE
  )

  merged <- .merge_landsat_scenes(sr_df, st_df)

  expect_equal(nrow(merged), 2)
  expect_equal(sort(merged$id), c("SCENE_A", "SCENE_B"))

  scene_a <- merged[merged$id == "SCENE_A", ]
  expect_true(is.na(scene_a$st_b10_href))
  expect_false(is.na(scene_a$sr_b4_href))

  scene_b <- merged[merged$id == "SCENE_B", ]
  expect_false(is.na(scene_b$st_b10_href))
  expect_true(is.na(scene_b$sr_b4_href))
})

test_that(".merge_landsat_scenes orders by cloud_cover ascending", {
  sr_df <- data.frame(
    base_id = c("A", "B"), datetime = c("2025-07-01", "2025-07-05"), cloud_cover = c(20, 5),
    st_b10_href = NA_character_, qa_pixel_href = NA_character_,
    sr_b4_href = "x", sr_b5_href = "x", sr_b6_href = "x", stringsAsFactors = FALSE
  )
  st_df <- sr_df[0, ]
  merged <- .merge_landsat_scenes(sr_df, st_df)
  expect_equal(merged$id, c("B", "A"))  # B ha meno nuvole (5 < 20)
})

test_that(".merge_landsat_scenes handles both inputs empty", {
  empty <- data.frame(base_id = character(0), datetime = character(0), cloud_cover = numeric(0),
                       st_b10_href = character(0), qa_pixel_href = character(0),
                       sr_b4_href = character(0), sr_b5_href = character(0), sr_b6_href = character(0),
                       stringsAsFactors = FALSE)
  merged <- .merge_landsat_scenes(empty, empty)
  expect_equal(nrow(merged), 0)
})


test_that(".find_asset_href resolves the real live server's common asset key names", {
  # struttura reale osservata dall'utente su landsatlook.usgs.gov/stac-server:
  # la chiave e' il nome comune ('lwir11'), non 'ST_B10' (che compare solo nel filename)
  assets <- list(
    lwir11 = list(href = "https://landsatlook.usgs.gov/data/.../LC09_L2SP_192029_20250625_20250627_02_T1_ST_B10.TIF"),
    red    = list(href = "https://landsatlook.usgs.gov/data/.../LC09_L2SP_192029_20250625_20250627_02_T1_SR_B4.TIF"),
    nir08  = list(href = "https://landsatlook.usgs.gov/data/.../LC09_L2SP_192029_20250625_20250627_02_T1_SR_B5.TIF"),
    swir16 = list(href = "https://landsatlook.usgs.gov/data/.../LC09_L2SP_192029_20250625_20250627_02_T1_SR_B6.TIF"),
    qa_pixel = list(href = "https://landsatlook.usgs.gov/data/.../LC09_L2SP_192029_20250625_20250627_02_T1_QA_PIXEL.TIF")
  )

  expect_true(grepl("ST_B10", .find_asset_href(assets, c("ST_B10", "st_b10", "lwir11"))))
  expect_true(grepl("SR_B4", .find_asset_href(assets, c("SR_B4", "sr_b4", "red"))))
  expect_true(grepl("SR_B5", .find_asset_href(assets, c("SR_B5", "sr_b5", "nir08"))))
  expect_true(grepl("SR_B6", .find_asset_href(assets, c("SR_B6", "sr_b6", "swir16"))))
  expect_true(grepl("QA_PIXEL", .find_asset_href(assets, c("QA_PIXEL", "qa_pixel"))))
})

test_that(".find_asset_href still resolves the native uppercase key names (older/alternate schema)", {
  assets <- list(
    ST_B10 = list(href = "https://example.com/ST_B10.tif"),
    SR_B4  = list(href = "https://example.com/SR_B4.tif")
  )
  expect_equal(.find_asset_href(assets, c("ST_B10", "st_b10", "lwir11")), "https://example.com/ST_B10.tif")
  expect_equal(.find_asset_href(assets, c("SR_B4", "sr_b4", "red")), "https://example.com/SR_B4.tif")
})

test_that(".find_asset_href returns NULL when no candidate key is present", {
  assets <- list(some_other_band = list(href = "https://example.com/x.tif"))
  expect_null(.find_asset_href(assets, c("ST_B10", "st_b10", "lwir11")))
})

test_that("shub_evalscript_lst returns a non-empty evalscript string", {
  es <- shub_evalscript_lst()
  expect_type(es, "character")
  expect_true(nchar(es) > 50)
  expect_true(grepl("ST_B10", es))
  expect_true(grepl("0.00341802", es))
})
