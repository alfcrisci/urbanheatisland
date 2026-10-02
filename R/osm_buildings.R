#' Scarica i poligoni dell'edificato da OpenStreetMap
#'
#' Interroga l'API Overpass (via il pacchetto \code{osmdata}) per tutti gli
#' elementi con tag \code{building=*} all'interno del bounding box dell'area
#' di interesse, e ritaglia il risultato sul confine esatto dell'AOI.
#' Nessuna registrazione o token richiesto — Overpass e' un servizio
#' pubblico gratuito, come Nominatim e lo STAC USGS gia' usati altrove nel
#' pacchetto.
#'
#' @param aoi Poligono \code{sf} che definisce l'area di interesse (in
#'   qualsiasi CRS proiettato; viene automaticamente convertito in EPSG:4326
#'   per l'interrogazione, come richiesto da Overpass).
#' @param timeout Timeout in secondi per la query Overpass. Default 90 —
#'   aumentalo per aree molto estese o dense (es. centri urbani grandi).
#' @param overpass_url URL del server Overpass. Default
#'   \code{"https://overpass-api.de/api/interpreter"}. Se il server
#'   principale e' sovraccarico, un mirror alternativo comune e'
#'   \code{"https://overpass.kumi.systems/api/interpreter"}.
#' @return Un oggetto \code{sf} di poligoni (edifici), ritagliato sull'AOI,
#'   nello stesso CRS di \code{aoi}. Puo' essere vuoto se non ci sono
#'   edifici mappati nell'area (raro nei centri urbani, piu' comune in aree
#'   rurali con copertura OSM scarsa).
#' @export
get_osm_buildings <- function(aoi, timeout = 90,
                               overpass_url = "https://overpass-api.de/api/interpreter") {
  if (!requireNamespace("osmdata", quietly = TRUE)) {
    stop("Il pacchetto 'osmdata' e' necessario per scaricare l'edificato da OSM. ",
         "Installalo con install.packages('osmdata').", call. = FALSE)
  }
  if (!requireNamespace("sf", quietly = TRUE)) stop("Il pacchetto 'sf' e' necessario.", call. = FALSE)

  aoi_4326 <- sf::st_transform(aoi, 4326)
  bbox <- sf::st_bbox(aoi_4326)

  query <- osmdata::opq(bbox = c(bbox["xmin"], bbox["ymin"], bbox["xmax"], bbox["ymax"]),
                         timeout = timeout) |>
    osmdata::add_osm_feature(key = "building")

  # osmdata gestisce internamente l'URL del server Overpass tramite le sue
  # opzioni interne; qui esponiamo overpass_url per completezza futura, ma
  # l'implementazione corrente usa il default di osmdata a meno di override
  # esplicito via osmdata::set_overpass_url().
  if (!identical(overpass_url, "https://overpass-api.de/api/interpreter")) {
    osmdata::set_overpass_url(overpass_url)
  }

  result <- osmdata::osmdata_sf(query)

  poly_list <- list()
  if (!is.null(result$osm_polygons) && nrow(result$osm_polygons) > 0) {
    poly_list[["poly"]] <- sf::st_geometry(result$osm_polygons)
  }
  if (!is.null(result$osm_multipolygons) && nrow(result$osm_multipolygons) > 0) {
    poly_list[["multipoly"]] <- sf::st_geometry(result$osm_multipolygons)
  }

  if (length(poly_list) == 0) {
    message("Nessun edificio trovato nell'area richiesta (bounding box).")
    return(sf::st_sf(geometry = sf::st_sfc(crs = 4326))[0, ])
  }

  buildings <- sf::st_sf(geometry = do.call(c, poly_list))
  buildings <- sf::st_make_valid(buildings)

  # ritaglio esatto sull'AOI (la query Overpass restituisce tutto il bbox,
  # che e' piu' ampio del confine reale se l'AOI non e' rettangolare)
  buildings <- suppressWarnings(sf::st_intersection(buildings, sf::st_geometry(aoi_4326)))
  buildings <- sf::st_transform(buildings, sf::st_crs(aoi))
  buildings
}

#' Rasterizza i poligoni dell'edificato in una covariata raster
#'
#' Converte i poligoni degli edifici (es. da \code{\link{get_osm_buildings}})
#' in un raster allineato a una griglia di riferimento, con diverse metriche
#' possibili — la piu' utile come covariata di "consumo di suolo" e'
#' \code{"coverage_fraction"} (frazione di superficie edificata per cella,
#' 0-1), direttamente comparabile a un layer NDBI o a un layer esterno tipo
#' ISPRA.
#'
#' @param buildings Poligoni \code{sf} degli edifici (es. da
#'   \code{\link{get_osm_buildings}}).
#' @param template \code{SpatRaster} che definisce la griglia di output
#'   (estensione, risoluzione, CRS) — tipicamente lo stesso raster LST o
#'   NDVI usato nel resto dell'analisi.
#' @param metric Una tra: \code{"coverage_fraction"} (default, frazione 0-1
#'   di area coperta da edifici per cella — analoga a un indice di
#'   impermeabilizzazione), \code{"count"} (numero di edifici il cui
#'   baricentro cade nella cella), \code{"binary"} (1 se almeno un edificio
#'   interseca la cella, 0 altrimenti).
#' @return Un \code{SpatRaster} sulla griglia di \code{template}.
#' @export
rasterize_buildings <- function(buildings, template,
                                 metric = c("coverage_fraction", "count", "binary")) {
  metric <- match.arg(metric)
  if (nrow(buildings) == 0) {
    out <- terra::rast(template[[1]])
    terra::values(out) <- 0
    names(out) <- paste0("edificato_", metric)
    return(out)
  }

  v <- terra::vect(sf::st_transform(buildings, terra::crs(template)))

  out <- switch(metric,
    coverage_fraction = terra::rasterize(v, template[[1]], cover = TRUE, background = 0),
    binary            = terra::rasterize(v, template[[1]], field = 1, background = 0),
    count             = {
      centroids <- terra::centroids(v)
      terra::rasterize(centroids, template[[1]], fun = "length", background = 0)
    }
  )
  names(out) <- paste0("edificato_", metric)
  out
}

#' Acquisizione end-to-end della covariata "edificato" da OSM
#'
#' Scarica i poligoni degli edifici OSM per l'AOI e li rasterizza
#' direttamente sulla griglia desiderata, in un'unica chiamata.
#'
#' @inheritParams get_osm_buildings
#' @param template \code{SpatRaster} che definisce la griglia di output
#'   (vedi \code{\link{rasterize_buildings}}).
#' @param metric Vedi \code{\link{rasterize_buildings}}. Default
#'   \code{"coverage_fraction"}.
#' @return Un \code{SpatRaster} sulla griglia di \code{template}, pronto per
#'   essere usato come covariata in \code{\link{define_rural_reference}}
#'   (come \code{impervious}) o in \code{\link{model_lst_drivers}}.
#' @export
get_builtup_covariate <- function(aoi, template, timeout = 90,
                                   metric = c("coverage_fraction", "count", "binary"),
                                   overpass_url = "https://overpass-api.de/api/interpreter") {
  metric <- match.arg(metric)
  buildings <- get_osm_buildings(aoi, timeout = timeout, overpass_url = overpass_url)
  rasterize_buildings(buildings, template, metric = metric)
}
