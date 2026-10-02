#' Estrae dai dati AirQino l'osservazione piu' vicina a un istante target
#'
#' Adatta l'output di \code{R.IBE.airqino::airqino_get_range()} (serie
#' storiche multi-riga per stazione) al formato "istantanea" richiesto da
#' \code{R.IBE.airqino::airqino_raster()} (una riga per stazione) —
#' selezionando, per ciascuna stazione, l'osservazione temporalmente più
#' vicina a un istante target (tipicamente l'orario di passaggio del
#' satellite, per confrontare la LST con le condizioni dell'aria nello
#' stesso momento anziché con una media giornaliera che le appiattirebbe).
#'
#' @param range_data Lista nominata \code{SMART{id} -> data.frame}, come
#'   restituita da \code{R.IBE.airqino::airqino_get_range()}. Ogni
#'   data.frame deve avere una colonna \code{time} (\code{POSIXct}).
#' @param target_datetime Istante target (\code{POSIXct}, o carattere in
#'   un formato interpretabile da \code{as.POSIXct()}). Se privo di fuso
#'   orario esplicito, viene assunto nel fuso indicato da \code{tz}.
#' @param tolerance_hours Scarto massimo (in ore) tra \code{target_datetime}
#'   e l'osservazione più vicina perché la stazione venga inclusa. Default
#'   \code{1.5} — oltre questa soglia la stazione viene scartata con un
#'   avviso, piuttosto che restituire un valore non rappresentativo.
#' @param tz Fuso orario da assumere se \code{target_datetime} non lo
#'   specifica esplicitamente. Default \code{"Europe/Rome"} (coerente con
#'   la conversione CET/CEST menzionata nella documentazione di
#'   \code{R.IBE.airqino}).
#' @return Lista nominata nello stesso formato atteso da
#'   \code{R.IBE.airqino::airqino_raster()}: \code{SMART{id} -> data.frame}
#'   a una riga. Le stazioni senza osservazioni entro \code{tolerance_hours}
#'   vengono omesse.
#' @export
airqino_snapshot_at <- function(range_data, target_datetime,
                                 tolerance_hours = 1.5, tz = "Europe/Rome") {
  if (is.character(target_datetime)) {
    target_datetime <- as.POSIXct(target_datetime, tz = tz)
  }
  tolerance_secs <- tolerance_hours * 3600

  risultato <- lapply(range_data, function(df) {
    if (is.null(df) || nrow(df) == 0 || !"time" %in% names(df)) return(NULL)
    scarti <- abs(as.numeric(difftime(df$time, target_datetime, units = "secs")))
    idx_migliore <- which.min(scarti)
    if (length(idx_migliore) == 0 || scarti[idx_migliore] > tolerance_secs) return(NULL)
    df[idx_migliore, , drop = FALSE]
  })

  risultato <- Filter(Negate(is.null), risultato)

  n_scartate <- length(range_data) - length(risultato)
  if (n_scartate > 0) {
    message(sprintf(
      "%d stazioni scartate: nessuna osservazione entro %.1f ore da %s.",
      n_scartate, tolerance_hours, format(target_datetime)
    ))
  }
  risultato
}

#' Interpola dati AirQino su una superficie continua (wrapper)
#'
#' Sottile wrapper attorno a \code{R.IBE.airqino::airqino_raster()} — vedi
#' quella funzione per la documentazione completa dei metodi disponibili
#' (\code{"idw"}, \code{"kriging"}, \code{"tin"}, \code{"spline"},
#' \code{"nearest"}). Presente qui principalmente per completezza
#' d'interfaccia all'interno di \code{urbanheatisland}; se preferisci puoi
#' chiamare direttamente \code{R.IBE.airqino::airqino_raster()}.
#'
#' @param all_raw Lista nominata \code{SMART{id} -> data.frame} (istantanea
#'   singola per stazione) — da \code{\link{airqino_snapshot_at}} per un
#'   confronto puntuale nel tempo con la LST, o da
#'   \code{R.IBE.airqino::airqino_get_last_value()} per l'ultimo valore
#'   disponibile.
#' @param ... Argomenti passati a \code{R.IBE.airqino::airqino_raster()}
#'   (\code{params}, \code{method}, \code{resolution}, \code{bbox}, ecc.).
#' @return Lista nominata \code{param -> SpatRaster}, in EPSG:4326.
#' @export
interpolate_airqino <- function(all_raw, ...) {
  if (!requireNamespace("R.IBE.airqino", quietly = TRUE)) {
    stop("Il pacchetto 'R.IBE.airqino' e' necessario. Installalo dal sorgente locale ",
         "(non e' su CRAN).", call. = FALSE)
  }
  R.IBE.airqino::airqino_raster(all_raw, ...)
}

#' Confronta la LST con una superficie AirQino interpolata (es. temperatura aria)
#'
#' Allinea una superficie interpolata AirQino (tipicamente prodotta da
#' \code{\link{interpolate_airqino}} o direttamente da
#' \code{R.IBE.airqino::airqino_raster()}) alla griglia della LST, poi
#' calcola correlazione, regressione lineare e — quando il parametro è una
#' temperatura — la differenza pixel per pixel.
#'
#' @param lst \code{SpatRaster} LST (gradi Celsius).
#' @param airqino_rast \code{SpatRaster} di un parametro AirQino interpolato
#'   (es. l'elemento \code{"temp"} della lista restituita da
#'   \code{\link{interpolate_airqino}}), in un CRS qualsiasi (tipicamente
#'   EPSG:4326).
#' @param param_name Nome del parametro (usato solo per etichettare gli
#'   output e decidere se calcolare la differenza — vedi
#'   \code{compute_difference}). Se \code{NULL}, usa \code{names(airqino_rast)}.
#' @param compute_difference Se \code{TRUE} (default), calcola anche
#'   \code{lst - airqino_aligned} pixel per pixel. Ha senso fisico solo se
#'   \code{airqino_rast} è anch'esso in gradi Celsius (es. temperatura
#'   dell'aria) — per altri parametri (NO2, PM10, ...) imposta
#'   \code{FALSE}, dato che LST e concentrazioni non sono nella stessa
#'   unità di misura e la "differenza" non sarebbe interpretabile.
#' @return Una lista con: \code{aligned} (\code{airqino_rast} riportato
#'   sulla griglia della LST), \code{difference} (\code{SpatRaster}, solo
#'   se \code{compute_difference = TRUE}), \code{correlation} (coefficiente
#'   di Pearson tra le celle non-NA di entrambi i raster), \code{model}
#'   (oggetto \code{lm}, \code{lst ~ airqino}), \code{n_cells} (numero di
#'   celle usate nel confronto).
#' @export
compare_lst_airqino <- function(lst, airqino_rast, param_name = NULL,
                                 compute_difference = TRUE) {
  if (is.null(param_name)) param_name <- names(airqino_rast)[1] %||% "airqino"

  aligned <- if (!terra::compareGeom(lst, airqino_rast, stopOnError = FALSE)) {
    align_rasters(lst, airqino_rast)
  } else {
    airqino_rast
  }
  names(aligned) <- param_name

  stack <- c(lst, aligned)
  names(stack) <- c("lst", param_name)
  df <- terra::as.data.frame(stack, na.rm = TRUE)

  n_cells <- nrow(df)
  if (n_cells < 3) {
    warning("Meno di 3 celle valide in comune tra LST e '", param_name,
             "' dopo l'allineamento — correlazione/regressione non affidabili.", call. = FALSE)
  }

  correlation <- if (n_cells >= 3) stats::cor(df$lst, df[[param_name]], use = "complete.obs") else NA_real_
  frm <- stats::as.formula(paste("lst ~", param_name))
  modello <- if (n_cells >= 3) stats::lm(frm, data = df) else NULL

  risultato <- list(
    aligned = aligned,
    correlation = correlation,
    model = modello,
    n_cells = n_cells
  )

  if (compute_difference) {
    differenza <- lst - aligned
    names(differenza) <- paste0("LST_meno_", param_name)
    risultato$difference <- differenza
  }

  risultato
}

#' Confronto end-to-end tra LST e parametri AirQino (temperatura, qualità dell'aria)
#'
#' Concatena l'intera catena: scarico dati storici AirQino nell'intervallo
#' indicato, estrazione dell'istantanea più vicina all'orario di passaggio
#' del satellite, interpolazione spaziale, e confronto con la LST — per uno
#' o più parametri in un'unica chiamata.
#'
#' @param lst \code{SpatRaster} LST (gradi Celsius) da confrontare.
#' @param lst_datetime Istante di acquisizione della LST (\code{POSIXct} o
#'   carattere), usato per selezionare le osservazioni AirQino
#'   temporalmente più vicine — vedi \code{\link{airqino_snapshot_at}}.
#' @param session Oggetto \code{airqino_session}, da
#'   \code{R.IBE.airqino::airqino_connect()}.
#' @param stations Vettore intero di numeri stazione.
#' @param start,end Intervallo di date per il download storico (vedi
#'   \code{R.IBE.airqino::airqino_get_range()}). Un intervallo di 1-2
#'   giorni centrato su \code{lst_datetime} è tipicamente sufficiente.
#' @param meta_df data.frame metadati stazioni, da
#'   \code{R.IBE.airqino::load_station_metadata()}.
#' @param params Vettore dei parametri da confrontare. Default
#'   \code{c("temp", "no2", "pm10", "pm25", "o3")} — i nomi standard di
#'   colonna usati da \code{R.IBE.airqino} (temperatura dell'aria e
#'   qualità dell'aria).
#' @param method Metodo di interpolazione spaziale (vedi
#'   \code{R.IBE.airqino::airqino_raster()}). Default \code{"idw"}.
#' @param tolerance_hours Vedi \code{\link{airqino_snapshot_at}}. Default \code{1.5}.
#' @param resolution Risoluzione della griglia di interpolazione, in gradi
#'   decimali (vedi \code{R.IBE.airqino::airqino_raster()}). Default \code{0.005}.
#' @param temperature_params Nomi tra \code{params} da trattare come
#'   temperature (per cui calcolare anche la differenza con la LST, non
#'   solo correlazione/regressione). Default \code{c("temp", "temp_int")}.
#' @param ... Altri argomenti passati a
#'   \code{R.IBE.airqino::airqino_raster()} (es. \code{idw_power},
#'   \code{bbox}).
#' @return Una lista nominata per parametro, ciascun elemento nel formato
#'   restituito da \code{\link{compare_lst_airqino}}. I parametri senza
#'   abbastanza stazioni valide vengono omessi con un avviso.
#' @export
compare_lst_with_airqino_range <- function(lst, lst_datetime, session, stations,
                                            start, end, meta_df,
                                            params = c("temp", "no2", "pm10", "pm25", "o3"),
                                            method = "idw", tolerance_hours = 1.5,
                                            resolution = 0.005,
                                            temperature_params = c("temp", "temp_int"),
                                            ...) {
  if (!requireNamespace("R.IBE.airqino", quietly = TRUE)) {
    stop("Il pacchetto 'R.IBE.airqino' e' necessario. Installalo dal sorgente locale ",
         "(non e' su CRAN).", call. = FALSE)
  }
  if (is.character(lst_datetime)) lst_datetime <- as.POSIXct(lst_datetime, tz = "Europe/Rome")

  range_data <- R.IBE.airqino::airqino_get_range(session, stations, start, end, meta_df)
  istantanea <- airqino_snapshot_at(range_data, lst_datetime, tolerance_hours = tolerance_hours)

  if (length(istantanea) < 3) {
    stop("Meno di 3 stazioni con un'osservazione entro ", tolerance_hours,
         " ore da ", format(lst_datetime), " — impossibile interpolare.", call. = FALSE)
  }

  superfici <- R.IBE.airqino::airqino_raster(
    istantanea, params = params, method = method, resolution = resolution, ...
  )

  risultati <- list()
  for (param in names(superfici)) {
    compute_diff <- param %in% temperature_params
    risultati[[param]] <- tryCatch(
      compare_lst_airqino(lst, superfici[[param]], param_name = param,
                           compute_difference = compute_diff),
      error = function(e) {
        message("Confronto LST-", param, " fallito: ", conditionMessage(e))
        NULL
      }
    )
  }

  Filter(Negate(is.null), risultati)
}
