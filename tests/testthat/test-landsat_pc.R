.crea_item_pc_finto <- function(id, cloud_cover, datetime, assets = NULL) {
  if (is.null(assets)) {
    assets <- list(
      lwir11 = list(href = paste0("https://esempio.it/", id, "_ST_B10.TIF")),
      red    = list(href = paste0("https://esempio.it/", id, "_SR_B4.TIF")),
      nir08  = list(href = paste0("https://esempio.it/", id, "_SR_B5.TIF")),
      swir16 = list(href = paste0("https://esempio.it/", id, "_SR_B6.TIF")),
      qa_pixel = list(href = paste0("https://esempio.it/", id, "_QA_PIXEL.TIF"))
    )
  }
  list(id = id, properties = list(datetime = datetime, `eo:cloud_cover` = cloud_cover), assets = assets)
}

test_that("pc_scene_table produce una tabella ordinata per nuvolosita' crescente", {
  items <- list(features = list(
    .crea_item_pc_finto("A", 45.0, "2025-07-10T10:00:00Z"),
    .crea_item_pc_finto("B", 5.2, "2025-07-01T10:00:00Z"),
    .crea_item_pc_finto("C", 12.0, "2025-07-20T10:00:00Z")
  ))

  tabella <- pc_scene_table(items)
  expect_equal(nrow(tabella), 3)
  expect_equal(tabella$id, c("B", "C", "A"))  # ordine per cloud_cover crescente
  expect_equal(tabella$cloud_cover, c(5.2, 12.0, 45.0))
})

test_that("pc_scene_table gestisce una lista di scene vuota", {
  vuoto <- pc_scene_table(list(features = list()))
  expect_equal(nrow(vuoto), 0)
  expect_equal(names(vuoto), c("indice", "id", "datetime", "cloud_cover"))
})

test_that("pc_asset_href estrae correttamente l'href per nome asset", {
  items <- list(features = list(.crea_item_pc_finto("SCENA1", 5.2, "2025-07-01T10:00:00Z")))

  expect_equal(pc_asset_href(items, 1, "lwir11"), "https://esempio.it/SCENA1_ST_B10.TIF")
  expect_equal(pc_asset_href(items, 1, "red"), "https://esempio.it/SCENA1_SR_B4.TIF")
})

test_that("pc_asset_href segnala un asset non disponibile elencando quelli presenti", {
  items <- list(features = list(.crea_item_pc_finto("SCENA1", 5.2, "2025-07-01T10:00:00Z")))
  expect_error(pc_asset_href(items, 1, "banda_inesistente"), "non trovato.*lwir11")
})

test_that("pc_asset_href segnala un indice fuori range", {
  items <- list(features = list(.crea_item_pc_finto("SCENA1", 5.2, "2025-07-01T10:00:00Z")))
  expect_error(pc_asset_href(items, 99, "lwir11"), "fuori range")
  expect_error(pc_asset_href(items, 0, "lwir11"), "fuori range")
})

test_that("il filtro client-side per copertura nuvolosa funziona come atteso", {
  # replica la logica di filtro usata internamente da pc_search_landsat
  features <- list(
    .crea_item_pc_finto("BASSA", 5.2, "2025-07-01"),
    .crea_item_pc_finto("ALTA", 45.0, "2025-07-10"),
    .crea_item_pc_finto("MEDIA", 18.0, "2025-07-20")
  )
  max_cloud_cover <- 20
  keep <- vapply(features, function(f) {
    cc <- f$properties[["eo:cloud_cover"]]
    !is.null(cc) && cc <= max_cloud_cover
  }, logical(1))
  filtrati <- features[keep]

  expect_equal(length(filtrati), 2)
  expect_equal(sort(sapply(filtrati, function(f) f$id)), c("BASSA", "MEDIA"))
})
