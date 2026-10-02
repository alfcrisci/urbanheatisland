#' Autenticazione USGS M2M API (token applicativo)
#'
#' Dal 2025 il download diretto da \code{landsatlook.usgs.gov/data/} (usato
#' da \code{\link{landsat_download_asset}}) richiede una sessione
#' autenticata: il viewer LandsatLook pubblico e' stato dismesso e le
#' richieste non autenticate vengono reindirizzate al login ERS. Questa
#' funzione autentica una sessione M2M usando un \emph{application token}
#' (non la password ERS — vedi
#' \url{https://www.usgs.gov/media/files/m2m-application-token-documentation}
#' per generarne uno dal tuo profilo ERS). Richiede inoltre che il tuo
#' account ERS abbia accesso M2M approvato (richiesta manuale via
#' \url{https://ers.cr.usgs.gov/profile/access}, tipicamente approvata in
#' un paio di giorni lavorativi).
#'
#' \strong{Nota implementativa}: questa funzione e le altre \code{landsat_m2m_*}
#' sono state scritte seguendo la documentazione ufficiale M2M, ma non e'
#' stato possibile testarle end-to-end in questo ambiente di sviluppo
#' (l'host \code{m2m.cr.usgs.gov} non e' raggiungibile dalla sandbox usata
#' per costruire il pacchetto). Se qualcosa non torna, apri una segnalazione
#' con l'errore esatto restituito dall'API — utile per correggere eventuali
#' scostamenti tra la documentazione e il comportamento reale del servizio.
#'
#' @param username Username ERS/USGS. Se \code{NULL}, letto dalla variabile
#'   d'ambiente \code{USGS_USERNAME}.
#' @param token Application token M2M (non la password). Se \code{NULL},
#'   letto dalla variabile d'ambiente \code{USGS_M2M_TOKEN}.
#' @param base_url URL base dell'API M2M. Default
#'   \code{"https://m2m.cr.usgs.gov/api/api/json/stable"}.
#' @return Una lista con \code{auth_token} (da passare come header
#'   \code{X-Auth-Token} alle altre funzioni \code{landsat_m2m_*}) e
#'   \code{base_url}.
#' @export
landsat_m2m_login <- function(username = NULL, token = NULL,
                               base_url = "https://m2m.cr.usgs.gov/api/api/json/stable") {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("Il pacchetto 'httr2' e' necessario. Installalo con install.packages('httr2').", call. = FALSE)
  }
  username <- username %||% Sys.getenv("USGS_USERNAME")
  token <- token %||% Sys.getenv("USGS_M2M_TOKEN")

  if (identical(username, "") || identical(token, "")) {
    stop(
      "Credenziali USGS M2M non fornite. Passa username/token oppure imposta ",
      "le variabili d'ambiente USGS_USERNAME / USGS_M2M_TOKEN. Il token e' un ",
      "'application token' generato dal profilo ERS, non la password.",
      call. = FALSE
    )
  }

  resp <- .m2m_request(base_url, "login-token", list(username = username, token = token))
  list(auth_token = resp$data, base_url = base_url)
}

#' Chiude la sessione M2M
#'
#' Buona pratica: chiude esplicitamente la sessione al termine del lavoro,
#' invece di lasciarla scadere per timeout.
#'
#' @param session Lista restituita da \code{\link{landsat_m2m_login}}.
#' @return Invisibile, \code{TRUE} se la chiusura ha avuto successo.
#' @export
landsat_m2m_logout <- function(session) {
  .m2m_request(session$base_url, "logout", list(), auth_token = session$auth_token)
  invisible(TRUE)
}

#' @keywords internal
.m2m_request <- function(base_url, endpoint, body, auth_token = NULL) {
  req <- httr2::request(paste0(base_url, "/", endpoint)) |> httr2::req_body_json(body)
  if (!is.null(auth_token)) {
    req <- httr2::req_headers(req, `X-Auth-Token` = auth_token)
  }
  # non lasciare che httr2 generi un errore R generico sui codici 4xx/5xx:
  # l'API M2M restituisce spesso un corpo JSON con errorCode/errorMessage
  # anche in caso di errore (es. 403), utile per capire la causa reale
  # (token scaduto, permessi insufficienti, dataset non valido, ecc.)
  req <- httr2::req_error(req, is_error = function(resp) FALSE)
  resp <- httr2::req_perform(req)
  status <- httr2::resp_status(resp)

  result <- tryCatch(httr2::resp_body_json(resp), error = function(e) NULL)

  if (is.null(result)) {
    body_text <- tryCatch(httr2::resp_body_string(resp), error = function(e) "<corpo non leggibile>")
    stop(sprintf("Risposta M2M API non interpretabile come JSON (HTTP %d) per l'endpoint '%s'.\nCorpo risposta (primi 500 caratteri): %s",
                  status, endpoint, substr(body_text, 1, 500)), call. = FALSE)
  }

  if (!is.null(result$errorCode) && !identical(result$errorCode, "")) {
    stop(sprintf("Errore M2M API [%s] (HTTP %d) sull'endpoint '%s': %s",
                  result$errorCode, status, endpoint, result$errorMessage %||% "nessun dettaglio"),
         call. = FALSE)
  }

  if (status >= 400) {
    stop(sprintf(
      "Errore HTTP %d dall'API M2M sull'endpoint '%s', senza errorCode nella risposta JSON. ",
      status, endpoint
    ), "Cause comuni di HTTP 403 su M2M: token di sessione scaduto (validita' tipica ~2 ore, ",
    "rifai il login), permessi M2M non ancora propagati dopo l'approvazione, o troppi download ",
    "in sospeso non ancora recuperati sull'account.", call. = FALSE)
  }

  result
}

#' Cerca scene Landsat nell'inventario USGS M2M
#'
#' A differenza di \code{\link{landsat_stac_search}} (che interroga il
#' catalogo STAC pubblico, senza autenticazione, ma con SR/ST come
#' collezioni separate), questa funzione interroga direttamente
#' l'inventario M2M. Nota: le funzioni di download bundle M2M non sono piu'
#' incluse nel pacchetto — per il download usa
#' \code{\link{get_landsat_lst_pc}} (Planetary Computer) o
#' \code{\link{get_sentinel3_lst}} (Copernicus Data Space).
#'
#' @param session Lista restituita da \code{\link{landsat_m2m_login}}.
#' @param aoi Poligono \code{sf} (in qualsiasi CRS proiettato; convertito
#'   automaticamente in EPSG:4326).
#' @param start_date,end_date Caratteri \code{"YYYY-MM-DD"}.
#' @param dataset Nome dataset M2M. Default \code{"landsat_ot_c2_l2"}
#'   (Landsat 8/9 OLI-TIRS Collection 2 Livello 2 — include sia le bande SR
#'   che la banda termica ST_B10 nello stesso prodotto).
#' @param max_cloud_cover Copertura nuvolosa massima (0-100). Default 20.
#' @param max_results Numero massimo di risultati. Default 50.
#' @return Un data.frame con \code{entity_id}, \code{display_id},
#'   \code{cloud_cover}, \code{acquisition_date}, ordinato per copertura
#'   nuvolosa crescente.
#' @export
landsat_m2m_scene_search <- function(session, aoi, start_date, end_date,
                                      dataset = "landsat_ot_c2_l2",
                                      max_cloud_cover = 20, max_results = 50) {
  bbox <- sf::st_bbox(sf::st_transform(aoi, 4326))

  body <- list(
    datasetName = dataset,
    sceneFilter = list(
      spatialFilter = list(
        filterType = "mbr",
        lowerLeft = list(latitude = unname(bbox["ymin"]), longitude = unname(bbox["xmin"])),
        upperRight = list(latitude = unname(bbox["ymax"]), longitude = unname(bbox["xmax"]))
      ),
      acquisitionFilter = list(start = start_date, end = end_date),
      cloudCoverFilter = list(min = 0, max = max_cloud_cover, includeUnknown = FALSE)
    ),
    maxResults = max_results
  )

  result <- .m2m_request(session$base_url, "scene-search", body, auth_token = session$auth_token)
  results <- result$data$results

  if (length(results) == 0) {
    return(data.frame(entity_id = character(0), display_id = character(0),
                       cloud_cover = numeric(0), acquisition_date = character(0)))
  }

  df <- do.call(rbind, lapply(results, function(r) {
    data.frame(
      entity_id = r$entityId,
      display_id = r$displayId %||% NA_character_,
      cloud_cover = r$cloudCover %||% NA_real_,
      acquisition_date = r$temporalCoverage$startDate %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }))
  df[order(df$cloud_cover), ]
}
