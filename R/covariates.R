#' Leggi una banda di riflettanza superficiale Landsat Collection 2 Livello 2
#'
#' Applica i fattori di scala ufficiali per le bande SR (Surface Reflectance)
#' di Landsat Collection 2 Livello 2 (diversi da quelli della banda termica
#' ST_B10).
#'
#' @param path Percorso del file GeoTIFF della banda SR (es. SR_B4, SR_B5, SR_B6).
#' @return Un \code{SpatRaster} con valori di riflettanza fisica (0-1 circa).
#' @export
read_landsat_sr <- function(path) {
  r <- terra::rast(path)
  r <- r * 0.0000275 - 0.2   # fattori di scala ufficiali Collection 2 SR
  r
}

#' Calcola l'NDVI (Normalized Difference Vegetation Index)
#'
#' @param nir Raster della banda infrarosso vicino (NIR), riflettanza fisica
#'   (es. da \code{\link{read_landsat_sr}} sulla banda SR_B5 per Landsat 8/9,
#'   o dalla banda B8 per Sentinel-2).
#' @param red Raster della banda rosso, riflettanza fisica (SR_B4 per
#'   Landsat 8/9, B4 per Sentinel-2).
#' @return Un \code{SpatRaster} NDVI, valori tra -1 e 1.
#' @export
compute_ndvi <- function(nir, red) {
  ndvi <- (nir - red) / (nir + red)
  names(ndvi) <- "NDVI"
  ndvi
}

#' Calcola l'NDBI (Normalized Difference Built-up Index)
#'
#' Indice proxy per il consumo di suolo/impermeabilizzazione, utile quando
#' non si dispone di un layer esterno di uso del suolo (es. ISPRA): valori
#' alti indicano superfici costruite/impermeabili, valori bassi vegetazione
#' o acqua.
#'
#' @param swir Raster della banda SWIR-1, riflettanza fisica (SR_B6 per
#'   Landsat 8/9, B11 per Sentinel-2).
#' @param nir Raster della banda NIR, riflettanza fisica (SR_B5 per
#'   Landsat 8/9, B8 per Sentinel-2).
#' @return Un \code{SpatRaster} NDBI, valori tra -1 e 1.
#' @export
compute_ndbi <- function(swir, nir) {
  ndbi <- (swir - nir) / (swir + nir)
  names(ndbi) <- "NDBI"
  ndbi
}

#' Acquisizione end-to-end di NDVI e NDBI da una scena Landsat
#'
#' Cerca la scena meno nuvolosa nel periodo indicato (via
#' \code{\link{landsat_stac_search}}), scarica le bande SR necessarie e
#' calcola NDVI e NDBI, allineati alla stessa griglia.
#'
#' @inheritParams landsat_stac_search
#' @param outdir Cartella locale per i file scaricati.
#' @return Un \code{SpatRaster} a due livelli: \code{NDVI}, \code{NDBI}.
#' @export
get_landsat_indices <- function(aoi, start_date, end_date, max_cloud_cover = 20,
                                 outdir = tempfile("landsat_idx_")) {
  scenes <- landsat_stac_search(aoi, start_date, end_date, max_cloud_cover = max_cloud_cover)
  if (nrow(scenes) == 0) {
    stop("Nessuna scena Landsat trovata per i criteri indicati.", call. = FALSE)
  }
  best <- scenes[1, ]
  if (any(is.na(c(best$sr_b4_href, best$sr_b5_href, best$sr_b6_href)))) {
    stop("La scena selezionata non ha tutte le bande SR necessarie (B4/B5/B6).", call. = FALSE)
  }
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  b4_dest <- file.path(outdir, paste0(best$id, "_SR_B4.tif"))
  b5_dest <- file.path(outdir, paste0(best$id, "_SR_B5.tif"))
  b6_dest <- file.path(outdir, paste0(best$id, "_SR_B6.tif"))
  landsat_download_asset(best$sr_b4_href, b4_dest)
  landsat_download_asset(best$sr_b5_href, b5_dest)
  landsat_download_asset(best$sr_b6_href, b6_dest)

  red  <- read_landsat_sr(b4_dest)
  nir  <- read_landsat_sr(b5_dest)
  swir <- read_landsat_sr(b6_dest)

  ndvi <- compute_ndvi(nir, red)
  ndbi <- compute_ndbi(swir, nir)
  c(ndvi, ndbi)
}

#' Definisci un riferimento rurale basato su covariate raster
#'
#' Alternativa a un semplice buffer di distanza: costruisce una maschera
#' "rurale" combinando soglie su una o piu' covariate raster (es. NDVI alto
#' E consumo di suolo/NDBI basso), indipendentemente dalla distanza dal
#' confine urbano. Il risultato puo' essere passato direttamente a
#' \code{\link{compute_uhi}} come \code{rural_ref}.
#'
#' @param ndvi Raster NDVI (opzionale se si usa solo \code{impervious}).
#' @param ndvi_min Soglia minima NDVI per considerare un pixel "rurale"
#'   (vegetazione sufficiente). Default \code{0.4}.
#' @param impervious Raster di frazione di suolo impermeabilizzato/consumato
#'   (es. layer ISPRA "Consumo di suolo", valori 0-1), oppure un raster NDBI
#'   se non si dispone di un layer esterno. Opzionale.
#' @param impervious_max Soglia massima di impermeabilizzazione/NDBI per
#'   considerare un pixel "rurale". Default \code{0.1}.
#' @param exclude_urban_mask Raster/maschera binaria opzionale (es. confine
#'   comunale bufferizzato, o un layer di area urbanizzata) da escludere a
#'   priori, per evitare che pixel con vegetazione dentro la citta' (parchi
#'   urbani) vengano conteggiati come "rurali".
#' @param min_area_cells Numero minimo di celle contigue "rurali" richieste
#'   per essere mantenute (rimuove pixel isolati/rumore). Default 1 (nessun
#'   filtro). Richiede il pacchetto \code{terra} con supporto per
#'   \code{patches()}.
#' @return Un \code{SpatRaster} binario (1 = rurale, 0 = non rurale, NA
#'   dove mancano i dati) sulla griglia della prima covariata fornita.
#' @export
define_rural_reference <- function(ndvi = NULL, ndvi_min = 0.4,
                                    impervious = NULL, impervious_max = 0.1,
                                    exclude_urban_mask = NULL,
                                    min_area_cells = 1) {
  if (is.null(ndvi) && is.null(impervious)) {
    stop("Fornire almeno una covariata tra `ndvi` e `impervious`.", call. = FALSE)
  }

  conditions <- list()
  if (!is.null(ndvi)) conditions[["ndvi"]] <- ndvi >= ndvi_min
  if (!is.null(impervious)) conditions[["impervious"]] <- impervious <= impervious_max

  ref <- conditions[[1]]
  if (length(conditions) > 1) {
    for (i in 2:length(conditions)) {
      cond_aligned <- if (!terra::compareGeom(ref, conditions[[i]], stopOnError = FALSE)) {
        align_rasters(ref, conditions[[i]], method = "near")
      } else conditions[[i]]
      ref <- ref & cond_aligned
    }
  }

  mask <- terra::ifel(ref, 1, 0)

  if (!is.null(exclude_urban_mask)) {
    urban_aligned <- if (!terra::compareGeom(mask, exclude_urban_mask, stopOnError = FALSE)) {
      align_rasters(mask, exclude_urban_mask)
    } else exclude_urban_mask
    mask <- terra::mask(mask, urban_aligned, maskvalues = 1, updatevalue = 0)
  }

  if (min_area_cells > 1) {
    patches <- terra::patches(mask, directions = 8, zeroAsNA = TRUE)
    sizes <- terra::freq(patches)
    small_patches <- sizes$value[sizes$count < min_area_cells]
    if (length(small_patches) > 0) {
      patches <- terra::subst(patches, from = small_patches, to = NA)
    }
    mask <- terra::ifel(is.na(patches), 0, 1)
  }

  names(mask) <- "rural_mask"
  mask
}

#' Modella i driver spaziali della LST usando covariate raster
#'
#' Analizza in che misura una o piu' covariate (NDVI, consumo di suolo,
#' altre variabili raster) spiegano il pattern spaziale della LST, pixel per
#' pixel. Utile sia a scopo esplicativo (quali fattori guidano l'isola di
#' calore) sia predittivo (ricostruzione/gap-filling di LST nei pixel
#' mascherati da nuvole).
#'
#' @param lst Raster LST (variabile risposta).
#' @param covariates Lista nominata di \code{SpatRaster}, uno per covariata
#'   (es. \code{list(ndvi = ndvi_rast, impervious = imperv_rast)}). Vengono
#'   allineati automaticamente alla griglia di \code{lst}.
#' @param method \code{"lm"} (regressione lineare, default, interpretabile)
#'   o \code{"rf"} (random forest, cattura relazioni non lineari, richiede
#'   il pacchetto \code{randomForest}).
#' @param predict_gaps Se \code{TRUE} (default), usa il modello per
#'   ricostruire i valori nei pixel di \code{lst} mancanti (es. mascherati
#'   da nuvole) dove le covariate sono disponibili.
#' @return Una lista con: \code{model} (oggetto lm o randomForest),
#'   \code{r_squared} (per \code{method = "lm"}), \code{importance} (per
#'   \code{method = "rf"}, importanza delle variabili),
#'   \code{predicted} (\code{SpatRaster} dei valori stimati dal modello su
#'   tutta l'area), \code{residuals} (\code{SpatRaster} osservato - stimato,
#'   solo dove LST era disponibile), \code{gap_filled} (\code{SpatRaster}
#'   LST originale con i buchi colmati dalla stima, se
#'   \code{predict_gaps = TRUE}).
#' @export
model_lst_drivers <- function(lst, covariates, method = c("lm", "rf"), predict_gaps = TRUE) {
  method <- match.arg(method)
  if (length(covariates) == 0 || is.null(names(covariates)) || any(names(covariates) == "")) {
    stop("`covariates` deve essere una lista nominata non vuota di SpatRaster.", call. = FALSE)
  }

  cov_aligned <- lapply(covariates, function(r) {
    if (!terra::compareGeom(lst, r, stopOnError = FALSE)) align_rasters(lst, r) else r
  })

  stack <- c(lst, terra::rast(cov_aligned))
  names(stack) <- c("lst", names(covariates))
  df <- terra::as.data.frame(stack, xy = TRUE, na.rm = FALSE)

  df_train <- df[stats::complete.cases(df[, c("lst", names(covariates))]), ]
  if (nrow(df_train) < length(covariates) + 2) {
    stop("Troppo pochi pixel con dati completi (LST + covariate) per stimare il modello.", call. = FALSE)
  }

  formula_str <- paste("lst ~", paste(names(covariates), collapse = " + "))
  frm <- stats::as.formula(formula_str)

  if (method == "lm") {
    fit <- stats::lm(frm, data = df_train)
    r_squared <- summary(fit)$r.squared
    importance <- NULL
  } else {
    if (!requireNamespace("randomForest", quietly = TRUE)) {
      stop("Il pacchetto 'randomForest' e' necessario per method = 'rf'. ",
           "Installalo con install.packages('randomForest').", call. = FALSE)
    }
    fit <- randomForest::randomForest(frm, data = df_train, importance = TRUE)
    r_squared <- NULL
    importance <- randomForest::importance(fit)
  }

  # predizione su tutti i pixel dove le covariate sono disponibili
  df_pred <- df[stats::complete.cases(df[, names(covariates), drop = FALSE]), ]
  pred_vals <- stats::predict(fit, newdata = df_pred)

  predicted <- lst
  terra::values(predicted) <- NA_real_
  cell_idx <- terra::cellFromXY(lst, as.matrix(df_pred[, c("x", "y")]))
  vals <- rep(NA_real_, terra::ncell(predicted))
  vals[cell_idx] <- as.numeric(pred_vals)
  terra::values(predicted) <- vals
  names(predicted) <- "LST_stimata"

  residuals_rast <- lst - predicted
  names(residuals_rast) <- "residui"

  result <- list(model = fit, r_squared = r_squared, importance = importance,
                  predicted = predicted, residuals = residuals_rast)

  if (predict_gaps) {
    gap_filled <- terra::cover(lst, predicted)
    names(gap_filled) <- "LST_gap_filled"
    result$gap_filled <- gap_filled
  }

  result
}
