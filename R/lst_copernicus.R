#' Cerca prodotti sul catalogo STAC di Copernicus Data Space Ecosystem
#'
#' La ricerca STAC su Copernicus Data Space e' pubblica (nessuna
#' autenticazione richiesta) — a differenza del download effettivo dei
#' dati, che passa da uno storage S3-compatibile e richiede credenziali
#' (vedi \code{\link{copernicus_download_asset}}).
#'
#' @param aoi Poligono \code{sf} che definisce l'area di interesse.
#' @param start_date,end_date Caratteri \code{"YYYY-MM-DD"}.
#' @param collection ID della collezione STAC. Esempi utili:
#'   \code{"sentinel-3-sl-2-lst-ntc"} (Sentinel-3 SLSTR Land Surface
#'   Temperature, risoluzione 1 km, prodotto consolidato non-time-critical),
#'   \code{"sentinel-3-sl-2-lst-nrt"} (stessa banda, versione near-real-time,
#'   piu' recente ma meno accurata), \code{"sentinel-2-l2a"} (ottico,
#'   10-20 m, per covariate NDVI/NDBI ad alta risoluzione).
#' @param max_cloud_cover Copertura nuvolosa massima (0-100). \code{NULL}
#'   (default) per non filtrare — utile per \code{sentinel-3-sl-2-lst-ntc},
#'   che non sempre espone \code{eo:cloud_cover}.
#' @param stac_url URL base dell'API STAC. Default
#'   \code{"https://stac.dataspace.copernicus.eu/v1"}.
#' @param limit Numero massimo di risultati. Default 50.
#' @return Un data.frame con \code{id}, \code{datetime}, \code{cloud_cover}
#'   (\code{NA} se non disponibile per la collezione), \code{assets} (lista
#'   colonna con gli asset grezzi di ogni item, da passare a
#'   \code{\link{copernicus_download_asset}}). Ordinato per copertura
#'   nuvolosa crescente se disponibile, altrimenti per data.
#' @export
copernicus_search <- function(aoi, start_date, end_date,
                               collection,
                               max_cloud_cover = NULL,
                               stac_url = "https://stac.dataspace.copernicus.eu/v1",
                               limit = 50) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("Il pacchetto 'httr2' e' necessario. Installalo con install.packages('httr2').", call. = FALSE)
  }
  bbox <- as.numeric(sf::st_bbox(sf::st_transform(aoi, 4326)))

  body <- list(
    collections = list(collection),
    bbox = bbox,
    datetime = paste0(start_date, "T00:00:00Z/", end_date, "T23:59:59Z"),
    limit = limit
  )
  if (!is.null(max_cloud_cover)) {
    body$query <- list(`eo:cloud_cover` = list(lt = max_cloud_cover))
  }

  req <- httr2::request(paste0(stac_url, "/search")) |> httr2::req_body_json(body)
  resp <- httr2::req_perform(req)
  result <- httr2::resp_body_json(resp)

  if (length(result$features) == 0) {
    return(data.frame(id = character(0), datetime = character(0), cloud_cover = numeric(0)))
  }

  df <- do.call(rbind, lapply(result$features, function(f) {
    data.frame(
      id = f$id,
      datetime = f$properties$datetime %||% NA_character_,
      cloud_cover = f$properties[["eo:cloud_cover"]] %||% NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  df$assets <- lapply(result$features, function(f) f$assets)

  if (all(is.na(df$cloud_cover))) {
    df[order(df$datetime, decreasing = TRUE), ]
  } else {
    df[order(df$cloud_cover), ]
  }
}

#' Scarica un asset da Copernicus Data Space (storage S3-compatibile)
#'
#' A differenza dello STAC USGS/landsatlook o di Planetary Computer, gli
#' asset di Copernicus Data Space non sono scaricabili con un semplice GET
#' HTTPS pubblico: risiedono su uno storage S3-compatibile
#' (\code{eodata.dataspace.copernicus.eu}) e richiedono credenziali S3 —
#' generabili gratuitamente e immediatamente dal tuo account Copernicus
#' (\url{https://documentation.dataspace.copernicus.eu/APIs/S3.html},
#' sezione "Manage S3 keys"), senza il processo di approvazione manuale
#' che invece serve per l'accesso M2M di USGS.
#'
#' @param assets Lista degli asset di un item (colonna \code{assets} del
#'   data.frame restituito da \code{\link{copernicus_search}}, un elemento
#'   per riga: \code{tabella$assets[[i]]}).
#' @param asset_name Nome dell'asset da scaricare (varia per collezione;
#'   ispeziona \code{names(assets)} per scoprire quelli disponibili).
#' @param dest File di destinazione.
#' @param s3_access_key,s3_secret_key Credenziali S3 Copernicus. Se
#'   \code{NULL}, lette dalle variabili d'ambiente
#'   \code{CDSE_S3_ACCESS_KEY} / \code{CDSE_S3_SECRET_KEY}.
#' @param s3_endpoint Endpoint S3-compatibile. Default
#'   \code{"eodata.dataspace.copernicus.eu"}.
#' @return \code{dest}, invisibile.
#' @export
copernicus_download_asset <- function(assets, asset_name, dest,
                                       s3_access_key = NULL, s3_secret_key = NULL,
                                       s3_endpoint = "eodata.dataspace.copernicus.eu") {
  href <- assets[[asset_name]]$href
  if (is.null(href)) {
    stop("Asset '", asset_name, "' non trovato. Asset disponibili: ",
         paste(names(assets), collapse = ", "), call. = FALSE)
  }

  if (!requireNamespace("aws.s3", quietly = TRUE)) {
    stop("Il pacchetto 'aws.s3' e' necessario per scaricare da Copernicus Data Space. ",
         "Installalo con install.packages('aws.s3').", call. = FALSE)
  }
  s3_access_key <- s3_access_key %||% Sys.getenv("CDSE_S3_ACCESS_KEY")
  s3_secret_key <- s3_secret_key %||% Sys.getenv("CDSE_S3_SECRET_KEY")
  if (identical(s3_access_key, "") || identical(s3_secret_key, "")) {
    stop(
      "Credenziali S3 Copernicus non fornite. Passa s3_access_key/s3_secret_key oppure ",
      "imposta le variabili d'ambiente CDSE_S3_ACCESS_KEY / CDSE_S3_SECRET_KEY. ",
      "Generale da https://eodata-s3keysmanager.dataspace.copernicus.eu/.",
      call. = FALSE
    )
  }

  parsed <- .parse_s3_href(href)

  Sys.setenv(AWS_S3_ENDPOINT = s3_endpoint)
  aws.s3::save_object(
    object = parsed$key, bucket = parsed$bucket, file = dest,
    region = "", key = s3_access_key, secret = s3_secret_key,
    base_url = s3_endpoint, use_https = TRUE
  )
  invisible(dest)
}

#' @keywords internal
.parse_s3_href <- function(href) {
  if (grepl("^s3://", href)) {
    # s3://bucket/key/con/path
    no_scheme <- sub("^s3://", "", href)
    bucket <- sub("/.*$", "", no_scheme)
    key <- sub("^[^/]+/", "", no_scheme)
  } else {
    # https://host/bucket/key/con/path (path-style) — il primo segmento del
    # path dopo l'host e' interpretato come nome del bucket
    no_scheme <- sub("^https?://[^/]+/", "", href)
    bucket <- sub("/.*$", "", no_scheme)
    key <- sub("^[^/]+/", "", no_scheme)
  }
  list(bucket = bucket, key = key)
}

#' Acquisizione end-to-end di LST da Sentinel-3 SLSTR (Copernicus)
#'
#' Risoluzione nativa 1 km — molto piu' grossolana di Landsat (30 m):
#' adatta ad analisi di pattern termico a scala di area urbana/regionale,
#' non a scala di isolato. Vantaggio: rivisitazione quasi giornaliera
#' (Sentinel-3A + 3B combinati), utile per compositi meno sensibili alle
#' nuvole rispetto a Landsat.
#'
#' @inheritParams copernicus_search
#' @param outdir Cartella locale per i file scaricati.
#' @param s3_access_key,s3_secret_key Vedi \code{\link{copernicus_download_asset}}.
#' @return Un \code{SpatRaster} LST in gradi Celsius.
#' @export
get_sentinel3_lst <- function(aoi, start_date, end_date,
                               collection = "sentinel-3-sl-2-lst-ntc",
                               outdir = tempfile("sentinel3_"),
                               s3_access_key = NULL, s3_secret_key = NULL) {
  tabella <- copernicus_search(aoi, start_date, end_date, collection = collection)
  if (nrow(tabella) == 0) {
    stop("Nessun prodotto Sentinel-3 trovato per i criteri indicati.", call. = FALSE)
  }
  best <- tabella[1, ]
  message("Uso prodotto ", best$id, " (", best$datetime, ")")

  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  lst_dest <- file.path(outdir, paste0(best$id, "_LST.tif"))

  # il nome dell'asset LST puo' variare; proviamo le varianti piu' comuni
  asset_candidates <- c("LST_in", "lst", "LST")
  asset_trovato <- asset_candidates[asset_candidates %in% names(tabella$assets[[1]])]
  if (length(asset_trovato) == 0) {
    stop("Nessun asset LST riconosciuto tra: ", paste(names(tabella$assets[[1]]), collapse = ", "),
         ". Verifica il nome esatto e passalo manualmente a copernicus_download_asset().", call. = FALSE)
  }

  copernicus_download_asset(tabella$assets[[1]], asset_trovato[1], lst_dest,
                             s3_access_key = s3_access_key, s3_secret_key = s3_secret_key)

  lst_k <- terra::rast(lst_dest)
  lst_c <- k_to_c(lst_k)
  names(lst_c) <- "LST_celsius"
  lst_c
}

#' Legge una banda Sentinel-2 L2A (JPEG2000) con lo scale factor corretto
#'
#' \strong{Attenzione a un dettaglio spesso trascurato}: dalla baseline di
#' processing ESA 04.00 (fine gennaio 2022), i prodotti Sentinel-2 L2A
#' includono un offset additivo (\code{BOA_ADD_OFFSET}, tipicamente
#' \code{-1000}) prima della divisione per 10000. Prodotti processati con
#' baseline precedenti non hanno questo offset. Questa funzione applica
#' \code{(DN + boa_add_offset) / 10000} — se non sei sicuro della baseline
#' del prodotto che stai usando, verifica nei metadati
#' (\code{MTD_MSIL2A.xml} nel pacchetto originale, o la proprieta' STAC
#' \code{s2:processing_baseline} se esposta dalla collezione) prima di
#' usare valori diversi da zero.
#'
#' @param path Percorso del file banda (\code{.jp2}), es. da
#'   \code{\link{copernicus_download_asset}}.
#' @param boa_add_offset Offset additivo prima della divisione per 10000.
#'   Default \code{0} (baseline precedente alla 04.00, o quando non e'
#'   nota la baseline effettiva). Usa \code{-1000} per prodotti con
#'   baseline >= 04.00 con offset gia' applicato nel DN.
#' @return Un \code{SpatRaster} di riflettanza fisica (0-1 circa).
#' @export
read_sentinel2_band <- function(path, boa_add_offset = 0) {
  r <- terra::rast(path)
  (r + boa_add_offset) / 10000
}

#' Acquisizione end-to-end di NDVI/NDBI da Sentinel-2 (Copernicus, 10-20 m)
#'
#' Risoluzione molto piu' fine di Landsat (10 m per NDVI, 20 m per le bande
#' coinvolte nell'NDBI) — utile come covariata di dettaglio anche quando la
#' LST viene da una fonte piu' grossolana (es. Sentinel-3, 1 km, vedi
#' \code{\link{get_sentinel3_lst}}), aggregata poi alla risoluzione della
#' LST con \code{\link{align_rasters}}.
#'
#' @inheritParams copernicus_search
#' @param outdir Cartella locale per i file scaricati.
#' @param boa_add_offset Vedi \code{\link{read_sentinel2_band}}. Default \code{0}.
#' @param s3_access_key,s3_secret_key Vedi \code{\link{copernicus_download_asset}}.
#' @return Un \code{SpatRaster} a due livelli: \code{NDVI} (10 m),
#'   \code{NDBI} (ricampionato a 10 m dalla risoluzione nativa di B11, 20 m).
#' @export
get_sentinel2_indices <- function(aoi, start_date, end_date,
                                   collection = "sentinel-2-l2a",
                                   max_cloud_cover = 20,
                                   outdir = tempfile("sentinel2_"),
                                   boa_add_offset = 0,
                                   s3_access_key = NULL, s3_secret_key = NULL) {
  tabella <- copernicus_search(aoi, start_date, end_date, collection = collection,
                                max_cloud_cover = max_cloud_cover)
  if (nrow(tabella) == 0) {
    stop("Nessuna scena Sentinel-2 trovata per i criteri indicati.", call. = FALSE)
  }
  best <- tabella[1, ]
  message("Uso scena ", best$id, " (", best$cloud_cover, "% nuvole)")

  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  assets <- tabella$assets[[1]]

  asset_candidates <- list(red = c("B04_10m", "B04"), nir = c("B08_10m", "B08"), swir = c("B11_20m", "B11"))
  trovare_asset <- function(candidati) {
    hit <- candidati[candidati %in% names(assets)]
    if (length(hit) == 0) {
      stop("Nessuno tra gli asset ", paste(candidati, collapse = "/"), " trovato. Asset disponibili: ",
           paste(names(assets), collapse = ", "), call. = FALSE)
    }
    hit[1]
  }

  red_dest  <- file.path(outdir, paste0(best$id, "_B04.jp2"))
  nir_dest  <- file.path(outdir, paste0(best$id, "_B08.jp2"))
  swir_dest <- file.path(outdir, paste0(best$id, "_B11.jp2"))

  copernicus_download_asset(assets, trovare_asset(asset_candidates$red), red_dest,
                             s3_access_key = s3_access_key, s3_secret_key = s3_secret_key)
  copernicus_download_asset(assets, trovare_asset(asset_candidates$nir), nir_dest,
                             s3_access_key = s3_access_key, s3_secret_key = s3_secret_key)
  copernicus_download_asset(assets, trovare_asset(asset_candidates$swir), swir_dest,
                             s3_access_key = s3_access_key, s3_secret_key = s3_secret_key)

  red  <- read_sentinel2_band(red_dest, boa_add_offset = boa_add_offset)
  nir  <- read_sentinel2_band(nir_dest, boa_add_offset = boa_add_offset)
  swir <- read_sentinel2_band(swir_dest, boa_add_offset = boa_add_offset)

  ndvi <- compute_ndvi(nir, red)
  # B11 e' a 20 m: ricampiona a 10 m (griglia di NDVI) prima dell'NDBI
  swir_10m <- if (!terra::compareGeom(nir, swir, stopOnError = FALSE)) align_rasters(nir, swir) else swir
  ndbi <- compute_ndbi(swir_10m, nir)

  c(ndvi, ndbi)
}
