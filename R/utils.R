#' Convert Kelvin to Celsius
#' @param k Numeric or SpatRaster in Kelvin.
#' @return Same type, in Celsius.
#' @export
k_to_c <- function(k) k - 273.15

#' Convert Celsius to Kelvin
#' @param c Numeric or SpatRaster in Celsius.
#' @return Same type, in Kelvin.
#' @export
c_to_k <- function(c) c + 273.15

#' Leggi un raster da una fonte di rete leggibile da GDAL
#'
#' Wrapper attorno a \code{terra::rast()} per sorgenti raster remote:
#' URL diretti a Cloud-Optimized GeoTIFF (COG), endpoint WMS/WCS, o
#' qualsiasi percorso che il file system virtuale di GDAL sappia aprire.
#' Utile per collegare il pacchetto a fonti dati non ancora coperte da un
#' modulo dedicato (es. un WCS regionale, un layer ISPRA via WMS, un COG
#' pubblico su un bucket) senza dover scrivere un downloader specifico.
#'
#' @param source URL o percorso della risorsa. Se non inizia gia' con uno
#'   dei prefissi del file system virtuale di GDAL (\code{/vsicurl/},
#'   \code{/vsis3/}, \code{WMS:}, ecc.) e sembra un URL HTTP(S) semplice,
#'   viene anteposto automaticamente \code{/vsicurl/} per permettere a GDAL
#'   di leggerlo in streaming senza scaricarlo prima su disco.
#' @param force_vsicurl Se \code{TRUE} (default), forza l'aggiunta del
#'   prefisso \code{/vsicurl/} per URL HTTP(S) semplici. Imposta a
#'   \code{FALSE} se \code{source} e' gia' un percorso GDAL completo (es.
#'   \code{"WCS:..."}, \code{"WMS:..."}) che non deve essere alterato.
#' @return Un \code{SpatRaster}.
#' @export
read_remote_raster <- function(source, force_vsicurl = TRUE) {
  is_plain_http <- grepl("^https?://", source) && !grepl("^/vsi|^WCS:|^WMS:|^NETCDF:", source)
  path <- if (force_vsicurl && is_plain_http) paste0("/vsicurl/", source) else source

  r <- tryCatch(
    terra::rast(path),
    error = function(e) {
      stop("Impossibile leggere il raster da '", source, "' (percorso GDAL tentato: '", path, "'). ",
           "Errore originale: ", conditionMessage(e), call. = FALSE)
    }
  )
  r
}

#' Ritaglia un raster su un poligono, controllando e allineando la proiezione
#'
#' Wrapper attorno a \code{terra::crop()} che verifica esplicitamente se il
#' poligono di ritaglio ha lo stesso CRS del raster, e lo riproietta prima
#' di procedere se necessario — invece di affidarsi alla riproiezione
#' implicita (e non sempre garantita a seconda della versione di terra) che
#' altrimenti avviene silenziosamente dentro \code{terra::crop()}.
#'
#' @param r Un \code{SpatRaster}.
#' @param polygon Un poligono \code{sf} o \code{SpatVector}, in qualsiasi CRS.
#' @param mask Se \code{TRUE} (default \code{FALSE}), applica anche
#'   \code{terra::mask()} dopo il crop (limita ai pixel dentro il poligono,
#'   non solo al suo bounding box).
#' @return Il \code{SpatRaster} ritagliato (ed eventualmente mascherato).
#' @export
crop_checked <- function(r, polygon, mask = FALSE) {
  v <- if (inherits(polygon, "sf")) terra::vect(polygon) else polygon

  if (!terra::same.crs(r, v)) {
    v <- terra::project(v, terra::crs(r))
  }

  out <- terra::crop(r, v)
  if (mask) out <- terra::mask(out, v)
  out
}

#' Quick diagnostic plot of a hotspot classification raster
#'
#' @param class_rast Output of \code{\link{classify_hotspots}}.
#' @param title Plot title.
#' @return Invisibly, the plotted \code{SpatRaster}.
#' @export
plot_hotspots <- function(class_rast, title = "Hotspot analysis (Gi*)") {
  cold_hot_palette <- c(
    "#2166AC", "#4393C3", "#92C5DE", "#F7F7F7",
    "#FDDBC7", "#D6604D", "#B2182B"
  )
  terra::plot(class_rast, col = cold_hot_palette, main = title, axes = FALSE)
  invisible(class_rast)
}
