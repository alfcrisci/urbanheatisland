#' Classifica un uso del suolo semplificato da covariate spettrali
#'
#' Costruisce un raster categorico a 3 classi (0 = altro/suolo nudo/acqua,
#' 1 = vegetazione, 2 = edificato) da NDVI ed eventualmente edificato/NDBI,
#' pronto per essere usato con \code{\link{compute_landscape_metrics}}.
#' Classificazione volutamente semplice, a soglie — se disponibile, un
#' prodotto ufficiale di uso del suolo (Corine Land Cover, ISPRA Consumo di
#' Suolo) da' risultati piu' accurati: in quel caso salta questa funzione e
#' passa direttamente il tuo raster classificato a
#' \code{\link{compute_landscape_metrics}}.
#'
#' @param ndvi Raster NDVI (es. da \code{\link{compute_ndvi}}).
#' @param ndvi_veg_min Soglia NDVI minima per la classe "vegetazione".
#'   Default \code{0.3}.
#' @param edificato Raster di frazione edificato/impermeabilizzazione
#'   (es. da \code{\link{get_builtup_covariate}} o un layer ISPRA), 0-1.
#'   Opzionale — se assente si usa \code{ndbi} al suo posto.
#' @param edificato_built_min Soglia minima di frazione edificato per la
#'   classe "edificato". Default \code{0.1} (10\%).
#' @param ndbi Raster NDBI (es. da \code{\link{compute_ndbi}}), alternativa
#'   a \code{edificato} se non disponibile un layer di consumo di suolo.
#' @param ndbi_built_min Soglia NDBI minima per la classe "edificato" (solo
#'   se si usa \code{ndbi} invece di \code{edificato}). Default \code{0}.
#' @return Un \code{SpatRaster} categorico a valori interi (0/1/2), sulla
#'   griglia di \code{ndvi}. In caso di sovrapposizione tra vegetazione ed
#'   edificato (es. alberature su superfici parzialmente impermeabili),
#'   prevale la classe "edificato".
#' @export
classify_landcover <- function(ndvi, ndvi_veg_min = 0.3,
                                edificato = NULL, edificato_built_min = 0.1,
                                ndbi = NULL, ndbi_built_min = 0) {
  veg_mask <- ndvi >= ndvi_veg_min

  if (!is.null(edificato)) {
    built_mask <- edificato >= edificato_built_min
  } else if (!is.null(ndbi)) {
    ndbi_aligned <- if (!terra::compareGeom(ndvi, ndbi, stopOnError = FALSE)) {
      align_rasters(ndvi, ndbi)
    } else ndbi
    built_mask <- ndbi_aligned >= ndbi_built_min
  } else {
    stop("Fornire almeno uno tra `edificato` e `ndbi` per identificare la classe edificato.", call. = FALSE)
  }

  built_aligned <- if (!terra::compareGeom(ndvi, built_mask, stopOnError = FALSE)) {
    align_rasters(ndvi, built_mask, method = "near")
  } else built_mask

  landcover <- terra::ifel(built_aligned, 2, terra::ifel(veg_mask, 1, 0))
  landcover <- terra::as.int(landcover)
  names(landcover) <- "landcover"
  landcover
}

#' Metriche di landscape ecology in finestra mobile (covariate raster)
#'
#' Calcola una o piu' metriche di ecologia del paesaggio (pacchetto
#' \code{landscapemetrics}) in una finestra mobile centrata su ogni cella,
#' producendo raster di covariate compatibili con
#' \code{\link{model_lst_drivers}} e \code{\link{define_rural_reference}}.
#' Cattura la *configurazione spaziale* del paesaggio (frammentazione,
#' compattezza, densita' dei margini) — un'informazione complementare a
#' NDVI/edificato, che descrivono solo la *composizione* (quanto c'e' di
#' ogni classe) e non come e' disposta.
#'
#' Se \code{classes_of_interest} e' \code{NULL}, calcola le metriche
#' \code{what} direttamente sul raster multi-classe (adatto a metriche di
#' eterogeneita' complessiva come Shannon diversity o contagion). Se
#' \code{classes_of_interest} e' specificato, per ciascuna classe
#' binarizza il raster (classe di interesse = 1, tutto il resto = 0) e
#' calcola le stesse metriche — utile per metriche specifiche di classe
#' come "densita' dei margini della vegetazione" o "densita' delle patch
#' edificate", che il pacchetto \code{landscapemetrics} non supporta
#' direttamente in finestra mobile a livello di classe.
#'
#' @param landcover Raster categorico a valori interi (es. da
#'   \code{\link{classify_landcover}} o un tuo layer di uso del suolo).
#' @param window_size Lato della finestra mobile quadrata, in celle
#'   (deve essere dispari). Default \code{5} — per un raster Landsat a
#'   30 m corrisponde a un vicinato di 150x150 m.
#' @param what Vettore di nomi di metriche a livello landscape del
#'   pacchetto \code{landscapemetrics} (vedi
#'   \code{landscapemetrics::list_lsm(level = "landscape")}). Default
#'   \code{c("lsm_l_ed", "lsm_l_pd", "lsm_l_shdi")} (densita' dei margini,
#'   densita' delle patch, diversita' di Shannon).
#' @param classes_of_interest Vettore opzionale di valori di classe (es.
#'   \code{c(1, 2)}) per cui calcolare le metriche specifiche per classe
#'   (vedi sopra). Se \code{NULL} (default), le metriche sono calcolate sul
#'   raster multi-classe cosi' com'e'.
#' @param class_labels Etichette opzionali per \code{classes_of_interest}
#'   (es. \code{c("vegetazione", "edificato")}), usate nei nomi dei layer
#'   di output. Se \code{NULL}, si usano \code{"classe1"}, \code{"classe2"}, ecc.
#' @return Un \code{SpatRaster} multi-livello, un layer per ogni
#'   combinazione metrica/classe, sulla griglia di \code{landcover}.
#' @export
compute_landscape_metrics <- function(landcover, window_size = 5,
                                       what = c("lsm_l_ed", "lsm_l_pd", "lsm_l_shdi"),
                                       classes_of_interest = NULL,
                                       class_labels = NULL) {
  if (!requireNamespace("landscapemetrics", quietly = TRUE)) {
    stop("Il pacchetto 'landscapemetrics' e' necessario. Installalo con ",
         "install.packages('landscapemetrics').", call. = FALSE)
  }
  if (window_size %% 2 == 0) stop("`window_size` deve essere un numero dispari di celle.", call. = FALSE)

  landcover <- terra::as.int(landcover)
  w <- matrix(1, nrow = window_size, ncol = window_size)

  if (is.null(classes_of_interest)) {
    res <- landscapemetrics::window_lsm(landcover, window = w, level = "landscape", what = what)
    layer_res <- res[[1]]
    stack_out <- terra::rast(unname(layer_res))
    names(stack_out) <- names(layer_res)
    return(stack_out)
  }

  if (is.null(class_labels)) class_labels <- paste0("classe", classes_of_interest)
  if (length(class_labels) != length(classes_of_interest)) {
    stop("`class_labels` deve avere la stessa lunghezza di `classes_of_interest`.", call. = FALSE)
  }

  all_vals <- sort(unique(stats::na.omit(as.numeric(terra::values(landcover)))))

  out_list <- list()
  out_names <- character(0)
  for (i in seq_along(classes_of_interest)) {
    cls <- classes_of_interest[i]
    lbl <- class_labels[i]
    # sfondo codificato come classe reale 0 (non NA): landscapemetrics
    # calcola densita' dei margini/patch correttamente solo con classi
    # esplicite, non con celle NA
    from_to <- cbind(from = all_vals, to = ifelse(all_vals == cls, 1, 0))
    bin <- terra::classify(landcover, from_to)

    res <- landscapemetrics::window_lsm(bin, window = w, level = "landscape", what = what)
    layer_res <- res[[1]]
    for (m in names(layer_res)) {
      key <- paste0(lbl, "_", sub("^lsm_l_", "", m))
      out_list[[key]] <- layer_res[[m]]
      out_names <- c(out_names, key)
    }
  }

  stack_out <- terra::rast(unname(out_list))
  names(stack_out) <- out_names
  stack_out
}

#' Covariate di landscape ecology, end-to-end, per vegetazione ed edificato
#'
#' Wrapper di convenienza: classifica l'uso del suolo da NDVI/edificato e
#' calcola in un'unica chiamata le metriche di configurazione spaziale piu'
#' usate in letteratura UHI — densita' dei margini (ED) e densita' delle
#' patch (PD) per la classe vegetazione ed edificato separatamente, piu' la
#' diversita' di Shannon (SHDI) sul mosaico completo.
#'
#' @inheritParams classify_landcover
#' @param window_size Vedi \code{\link{compute_landscape_metrics}}. Default \code{5}.
#' @return Un \code{SpatRaster} multi-livello: \code{vegetazione_ed},
#'   \code{vegetazione_pd}, \code{edificato_ed}, \code{edificato_pd},
#'   \code{shdi} — pronto per \code{\link{model_lst_drivers}}.
#' @export
get_landscape_covariates <- function(ndvi, ndvi_veg_min = 0.3,
                                      edificato = NULL, edificato_built_min = 0.1,
                                      ndbi = NULL, ndbi_built_min = 0,
                                      window_size = 5) {
  landcover <- classify_landcover(ndvi, ndvi_veg_min = ndvi_veg_min,
                                   edificato = edificato, edificato_built_min = edificato_built_min,
                                   ndbi = ndbi, ndbi_built_min = ndbi_built_min)

  per_classe <- compute_landscape_metrics(
    landcover, window_size = window_size,
    what = c("lsm_l_ed", "lsm_l_pd"),
    classes_of_interest = c(1, 2),
    class_labels = c("vegetazione", "edificato")
  )

  shdi <- compute_landscape_metrics(landcover, window_size = window_size, what = "lsm_l_shdi")

  out <- c(per_classe, shdi)
  names(out)[terra::nlyr(out)] <- "shdi"
  out
}
