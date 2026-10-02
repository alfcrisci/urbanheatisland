#' Calcola l'Urban Thermal Field Variance Index (UTFVI)
#'
#' L'UTFVI (Liu & Zhang, 2011) e' un indice ampiamente usato per valutare
#' l'effetto isola di calore urbana e la qualita' dell'ambiente termico:
#'
#' \deqn{UTFVI = \frac{LST - \overline{LST}}{\overline{LST}}}
#'
#' Valori positivi indicano pixel piu' caldi della media (peggiore qualita'
#' termica), valori negativi pixel piu' freschi (migliore qualita').
#'
#' \strong{Nota sulle unita' di misura}: la formula richiede LST in Kelvin —
#' non e' invariante per unita' di misura (a differenza, ad esempio, di un
#' UHI calcolato come semplice differenza). Poiche' in questo pacchetto la
#' LST e' convenzionalmente in gradi Celsius (vedi \code{\link{read_landsat_lst}},
#' \code{\link{get_modis_lst}}), questa funzione converte automaticamente in
#' Kelvin prima del calcolo — a meno che \code{input_unit = "kelvin"}.
#'
#' @param lst Un \code{SpatRaster} LST, a uno o piu' livelli (es. uno stack
#'   con un livello per data, o un singolo composito estivo).
#' @param input_unit Unita' di \code{lst}: \code{"celsius"} (default,
#'   convenzione del pacchetto) o \code{"kelvin"}.
#' @param reference Come calcolare la LST media di riferimento
#'   (\eqn{\overline{LST}}): \code{"per_layer"} (default — ogni livello
#'   usa la propria media spaziale, il calcolo standard in letteratura,
#'   adatto quando \code{lst} e' uno stack temporale e si vuole l'anomalia
#'   di ciascuna data rispetto a se stessa), \code{"global"} (una singola
#'   media calcolata su tutti i livelli insieme, utile per confrontare piu'
#'   date rispetto a una baseline comune), oppure un valore numerico fisso
#'   fornito direttamente (nell'unita' indicata da \code{input_unit}).
#' @return Un \code{SpatRaster} UTFVI (adimensionale), stessa struttura di
#'   \code{lst} (stesso numero di livelli).
#' @export
compute_utfvi <- function(lst, input_unit = c("celsius", "kelvin"),
                           reference = "per_layer") {
  input_unit <- match.arg(input_unit)

  lst_k <- if (input_unit == "celsius") c_to_k(lst) else lst

  if (is.numeric(reference)) {
    ref_k <- if (input_unit == "celsius") reference + 273.15 else reference
    utfvi <- (lst_k - ref_k) / ref_k
  } else {
    reference <- match.arg(reference, c("per_layer", "global"))
    if (reference == "per_layer") {
      means <- terra::global(lst_k, "mean", na.rm = TRUE)[, 1]
      utfvi <- lst_k
      for (i in seq_len(terra::nlyr(lst_k))) {
        utfvi[[i]] <- (lst_k[[i]] - means[i]) / means[i]
      }
    } else {
      ref_k <- terra::global(lst_k, "mean", na.rm = TRUE)
      ref_k <- mean(ref_k[, 1], na.rm = TRUE)
      utfvi <- (lst_k - ref_k) / ref_k
    }
  }

  names(utfvi) <- paste0(names(lst), "_UTFVI")
  utfvi
}

#' Classifica l'UTFVI in bande standard di qualita' termica
#'
#' Applica le soglie standard di letteratura (Liu & Zhang, 2011) per
#' classificare l'UTFVI in 6 classi di qualita' dell'ambiente termico /
#' intensita' dell'effetto isola di calore.
#'
#' @param utfvi Un \code{SpatRaster} UTFVI, come da \code{\link{compute_utfvi}}.
#' @param breaks Soglie di classificazione. Default gli standard di
#'   letteratura: \code{c(-Inf, 0, 0.005, 0.010, 0.015, 0.020, Inf)}.
#' @param labels Etichette delle 6 classi, dalla migliore alla peggiore
#'   qualita' termica. Default \code{c("Eccellente", "Buona", "Normale",
#'   "Scarsa", "Peggiore", "Pessima")}.
#' @return Un \code{SpatRaster} categorico con le classi di qualita' termica.
#' @export
classify_utfvi <- function(utfvi,
                            breaks = c(-Inf, 0, 0.005, 0.010, 0.015, 0.020, Inf),
                            labels = c("Eccellente", "Buona", "Normale",
                                       "Scarsa", "Peggiore", "Pessima")) {
  if (length(breaks) != length(labels) + 1) {
    stop("`breaks` deve avere un elemento in piu' di `labels`.", call. = FALSE)
  }

  rcl <- cbind(breaks[-length(breaks)], breaks[-1], seq_along(labels) - 1)
  classified <- terra::classify(utfvi, rcl, include.lowest = TRUE, right = TRUE)

  lv_df <- data.frame(id = seq_along(labels) - 1, classe = labels)
  levels(classified) <- rep(list(lv_df), terra::nlyr(classified))
  names(classified) <- paste0(names(utfvi), "_classe")
  classified
}

#' Riepilogo della distribuzione delle classi UTFVI
#'
#' Calcola l'area (in km\eqn{^2}, se il CRS e' proiettato in metri) e la
#' percentuale coperta da ciascuna classe di qualita' termica.
#'
#' @param utfvi_class Raster categorico, come da \code{\link{classify_utfvi}}.
#'   Se multi-livello, il riepilogo e' calcolato per ciascun livello
#'   separatamente e i risultati sono combinati con una colonna \code{layer}.
#' @return Un \code{data.frame} con colonne \code{layer} (se multi-livello),
#'   \code{classe}, \code{n_celle}, \code{area_km2}, \code{percentuale}.
#' @export
summarize_utfvi <- function(utfvi_class) {
  cell_area_km2 <- prod(terra::res(utfvi_class)) / 1e6

  summarize_layer <- function(r, layer_name) {
    freq_tab <- terra::freq(r)
    total_cells <- sum(freq_tab$count)
    out <- data.frame(
      classe = freq_tab$value,   # terra::freq() su raster categorico restituisce gia' l'etichetta
      n_celle = freq_tab$count,
      area_km2 = freq_tab$count * cell_area_km2,
      percentuale = 100 * freq_tab$count / total_cells
    )
    if (!is.null(layer_name)) out <- cbind(layer = layer_name, out)
    out
  }

  if (terra::nlyr(utfvi_class) == 1) {
    summarize_layer(utfvi_class, NULL)
  } else {
    do.call(rbind, lapply(seq_len(terra::nlyr(utfvi_class)), function(i) {
      summarize_layer(utfvi_class[[i]], names(utfvi_class)[i])
    }))
  }
}
