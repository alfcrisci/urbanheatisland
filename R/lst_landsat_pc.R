#' Cerca scene Landsat su Microsoft Planetary Computer (via rstac)
#'
#' Alternativa a \code{\link{landsat_stac_search}} (STAC USGS, che richiede
#' M2M per il download) e all'accesso M2M USGS (che richiede
#' approvazione manuale dell'account): Planetary Computer
#' ospita gli stessi dati Landsat Collection 2 Livello 2, in un'unica
#' collezione che include sia le bande SR che la banda termica ST_B10
#' (a differenza dello STAC USGS, dove sono separate), e non richiede
#' registrazione per l'uso base — solo una "firma" temporanea degli URL
#' (SAS token), gestita automaticamente da questa funzione.
#'
#' @param aoi Poligono \code{sf} che definisce l'area di interesse.
#' @param start_date,end_date Caratteri \code{"YYYY-MM-DD"}.
#' @param max_cloud_cover Copertura nuvolosa massima (0-100). Default 20.
#' @param collection Collezione STAC su Planetary Computer. Default
#'   \code{"landsat-c2-l2"} (Landsat Collection 2 Livello 2, tutte le
#'   bande in un'unica collezione).
#' @param stac_url URL base dell'API STAC di Planetary Computer. Default
#'   \code{"https://planetarycomputer.microsoft.com/api/stac/v1"}.
#' @param limit Numero massimo di risultati richiesti al server prima del
#'   filtro per nuvolosita' (il filtro stesso e' applicato lato client).
#'   Default 50.
#' @return Un oggetto \code{doc_items} di rstac (classe \code{STACItemCollection}),
#'   con gli asset gia' firmati (URL temporanei pronti per il download
#'   diretto) e gia' filtrato per \code{max_cloud_cover}. Passalo a
#'   \code{\link{pc_scene_table}} per una tabella riepilogativa leggibile, o
#'   a \code{\link{pc_asset_href}} per estrarre l'URL di una banda specifica.
#' @export
pc_search_landsat <- function(aoi, start_date, end_date, max_cloud_cover = 20,
                               collection = "landsat-c2-l2",
                               stac_url = "https://planetarycomputer.microsoft.com/api/stac/v1",
                               limit = 50) {
  if (!requireNamespace("rstac", quietly = TRUE)) {
    stop("Il pacchetto 'rstac' e' necessario. Installalo con install.packages('rstac').", call. = FALSE)
  }
  bbox <- as.numeric(sf::st_bbox(sf::st_transform(aoi, 4326)))
  datetime_range <- paste0(start_date, "T00:00:00Z/", end_date, "T23:59:59Z")

  items <- rstac::stac(stac_url) |>
    rstac::stac_search(collections = collection, bbox = bbox, datetime = datetime_range, limit = limit) |>
    rstac::get_request()

  # filtro lato client per copertura nuvolosa (piu' robusto del supporto
  # server-side della Query Extension, che puo' variare)
  keep <- vapply(items$features, function(f) {
    cc <- f$properties[["eo:cloud_cover"]]
    !is.null(cc) && cc <= max_cloud_cover
  }, logical(1))
  items$features <- items$features[keep]

  if (length(items$features) == 0) {
    warning("Nessuna scena trovata con copertura nuvolosa <= ", max_cloud_cover, "%.", call. = FALSE)
    return(items)
  }

  # firma gli asset: senza questo passaggio gli href puntano a blob privati
  # e restituiscono 404 al download
  rstac::items_sign_planetary_computer(items)
}

#' Tabella riepilogativa delle scene trovate su Planetary Computer
#'
#' @param items Oggetto \code{doc_items}, da \code{\link{pc_search_landsat}}.
#' @return Un data.frame con \code{indice} (da usare in
#'   \code{\link{pc_asset_href}}), \code{id}, \code{datetime},
#'   \code{cloud_cover}, ordinato per nuvolosita' crescente.
#' @export
pc_scene_table <- function(items) {
  if (length(items$features) == 0) {
    return(data.frame(indice = integer(0), id = character(0),
                       datetime = character(0), cloud_cover = numeric(0)))
  }
  tabella <- do.call(rbind, lapply(seq_along(items$features), function(i) {
    f <- items$features[[i]]
    data.frame(indice = i, id = f$id, datetime = f$properties$datetime %||% NA_character_,
               cloud_cover = f$properties[["eo:cloud_cover"]] %||% NA_real_, stringsAsFactors = FALSE)
  }))
  tabella[order(tabella$cloud_cover), ]
}

#' Estrae l'URL (gia' firmato) di un asset da una scena Planetary Computer
#'
#' @param items Oggetto \code{doc_items}, da \code{\link{pc_search_landsat}}.
#' @param indice Indice della scena (colonna \code{indice} di
#'   \code{\link{pc_scene_table}}).
#' @param asset_name Nome dell'asset. Per \code{landsat-c2-l2}: \code{"lwir11"}
#'   (banda termica, equivalente a ST_B10), \code{"red"} (SR_B4),
#'   \code{"nir08"} (SR_B5), \code{"swir16"} (SR_B6), \code{"qa_pixel"}
#'   (maschera qualita').
#' @return URL firmato (stringa), pronto per il download diretto con
#'   \code{\link{landsat_download_asset}}.
#' @export
pc_asset_href <- function(items, indice, asset_name) {
  if (indice > length(items$features) || indice < 1) {
    stop("Indice ", indice, " fuori range (scene disponibili: ", length(items$features), ").", call. = FALSE)
  }
  href <- items$features[[indice]]$assets[[asset_name]]$href
  if (is.null(href)) {
    disponibili <- names(items$features[[indice]]$assets)
    stop("Asset '", asset_name, "' non trovato. Asset disponibili: ",
         paste(disponibili, collapse = ", "), call. = FALSE)
  }
  href
}

#' Acquisizione end-to-end di LST Landsat via Planetary Computer
#'
#' Alternativa a \code{\link{get_landsat_lst}} (STAC USGS, download rotto
#' dal 2025) e all'accesso M2M USGS (richiede approvazione manuale
#' dell'account): nessuna registrazione necessaria, solo la firma
#' automatica degli URL gestita da \code{\link{pc_search_landsat}}.
#'
#' @inheritParams pc_search_landsat
#' @param outdir Cartella locale per i file scaricati.
#' @param qc_filter Applica il mascheramento QA_PIXEL. Default \code{TRUE}.
#' @return Un \code{SpatRaster} LST in gradi Celsius.
#' @export
get_landsat_lst_pc <- function(aoi, start_date, end_date, max_cloud_cover = 20,
                                collection = "landsat-c2-l2",
                                outdir = tempfile("landsat_pc_"), qc_filter = TRUE) {
  items <- pc_search_landsat(aoi, start_date, end_date, max_cloud_cover = max_cloud_cover,
                              collection = collection)
  tabella <- pc_scene_table(items)
  if (nrow(tabella) == 0) {
    stop("Nessuna scena trovata su Planetary Computer per i criteri indicati.", call. = FALSE)
  }
  best <- tabella[1, ]
  message(sprintf("Uso scena %s (%.1f%% nuvole)", best$id, best$cloud_cover))

  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  st_href <- pc_asset_href(items, best$indice, "lwir11")
  st_dest <- file.path(outdir, paste0(best$id, "_ST_B10.tif"))
  landsat_download_asset(st_href, st_dest)

  qa_dest <- NULL
  if (qc_filter) {
    qa_href <- tryCatch(pc_asset_href(items, best$indice, "qa_pixel"), error = function(e) NULL)
    if (!is.null(qa_href)) {
      qa_dest <- file.path(outdir, paste0(best$id, "_QA_PIXEL.tif"))
      landsat_download_asset(qa_href, qa_dest)
    }
  }

  read_landsat_lst(st_dest, qa_path = qa_dest)
}

#' Acquisizione end-to-end di NDVI/NDBI via Planetary Computer
#'
#' Equivalente a \code{\link{get_landsat_indices}} (STAC USGS) ma su
#' Planetary Computer — nessuna registrazione necessaria.
#'
#' @inheritParams pc_search_landsat
#' @param outdir Cartella locale per i file scaricati.
#' @return Un \code{SpatRaster} a due livelli: \code{NDVI}, \code{NDBI}.
#' @export
get_planetary_computer_indices <- function(aoi, start_date, end_date, max_cloud_cover = 20,
                                            collection = "landsat-c2-l2",
                                            outdir = tempfile("landsat_pc_idx_")) {
  items <- pc_search_landsat(aoi, start_date, end_date, max_cloud_cover = max_cloud_cover,
                              collection = collection)
  tabella <- pc_scene_table(items)
  if (nrow(tabella) == 0) {
    stop("Nessuna scena trovata su Planetary Computer per i criteri indicati.", call. = FALSE)
  }
  best <- tabella[1, ]
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  red_href  <- pc_asset_href(items, best$indice, "red")
  nir_href  <- pc_asset_href(items, best$indice, "nir08")
  swir_href <- pc_asset_href(items, best$indice, "swir16")

  red_dest  <- file.path(outdir, paste0(best$id, "_SR_B4.tif"))
  nir_dest  <- file.path(outdir, paste0(best$id, "_SR_B5.tif"))
  swir_dest <- file.path(outdir, paste0(best$id, "_SR_B6.tif"))
  landsat_download_asset(red_href, red_dest)
  landsat_download_asset(nir_href, nir_dest)
  landsat_download_asset(swir_href, swir_dest)

  ndvi <- compute_ndvi(read_landsat_sr(nir_dest), read_landsat_sr(red_dest))
  ndbi <- compute_ndbi(read_landsat_sr(swir_dest), read_landsat_sr(nir_dest))
  c(ndvi, ndbi)
}
