# urbanheatisland

Land Surface Temperature (LST) acquisition, interpolated air temperature
surfaces, urban heat island (UHI) metrics, and local spatial hotspot
analysis on raster grids — built for CNR-IBE biometeorology and urban heat
island workflows (e.g. MIRifICUS).

## What it does

| Module | File | Key functions |
|---|---|---|
| **LST acquisition (remote)** | `R/lst_appeears.R` | `get_modis_lst()`, `lst_appeears_login()`, `lst_appeears_submit_task()`, `lst_appeears_wait()`, `lst_appeears_download()` |
| **LST processing (local)** | `R/lst_local.R` | `read_lst_local()`, `decode_modis_qc()`, `mask_lst_by_qc()`, `mosaic_lst()`, `composite_lst()` |
| **Air temperature interpolation** | `R/interpolate_airtemp.R` | `interpolate_airtemp()` (idw / kriging / tps / tin / nn), `cv_interpolate_airtemp()` |
| **UHI & downscaling** | `R/uhi_downscale.R` | `compute_uhi()`, `align_rasters()`, `downscale_airtemp()` (regression / regression-kriging / random forest) |
| **Hotspot analysis** | `R/hotspot_analysis.R` | `gi_star_raster()` (Getis-Ord Gi\*), `local_moran_raster()` (Anselin Local Moran's I), `classify_hotspots()` |
| **Utilities** | `R/utils.R` | `k_to_c()`, `c_to_k()`, `plot_hotspots()` |

## Installation

```r
# Core dependencies (available on CRAN)
install.packages(c("terra", "sf", "gstat", "spdep", "jsonlite"))

# Optional, enable extra features:
install.packages(c("httr2", "fields", "randomForest", "interp", "geojsonsf"))

# Then install this package from source
install.packages("urbanheatisland_0.6.1.tar.gz", repos = NULL, type = "source")
```

`httr2` is required only for the AppEEARS download functions (`get_modis_lst()`
and friends) — if you only work with local LST files, you don't need it.
`fields` is required only for `method = "tps"` interpolation. `randomForest`
is required only for `downscale_airtemp(method = "rf")`.

## Typical workflow

```r
library(urbanheatisland)

# 1. Get LST — either download...
lst <- get_modis_lst(
  aoi = my_aoi_polygon,           # sf polygon, EPSG:4326
  start_date = "06-01-2026",
  end_date   = "08-31-2026",
  product    = "MOD11A2.061"
)

# ...or read a file you already have
lst <- read_lst_local("MOD11A2.A2026153.LST_Day_1km.tif",
                       qc_path = "MOD11A2.A2026153.QC_Day.tif")

# 2. Interpolate station air temperature to a raster
air_temp <- interpolate_airtemp(stations, "temp_c", template = lst, method = "kriging")

# Compare interpolation methods first, if unsure:
cv_interpolate_airtemp(stations, "temp_c", template = lst,
                        methods = c("idw", "kriging", "nn"))

# 3. Urban heat island intensity relative to a rural reference zone
uhi <- compute_uhi(lst, rural_ref = rural_polygon)

# 4. Downscale coarse air temperature using LST as covariate
air_temp_fine <- downscale_airtemp(air_temp, lst, method = "regression")

# 5. Hotspot analysis
gi <- gi_star_raster(lst, d = 1500)          # e.g. 1.5 km neighbourhood
hotspots <- classify_hotspots(gi, apply_fdr = TRUE)
plot_hotspots(hotspots)

# Local Moran's I for cluster vs. outlier typing (HH/LL/HL/LH)
lisa <- local_moran_raster(lst, d = 1500)
```

## Authentication for MODIS download

`get_modis_lst()` uses the NASA AppEEARS REST API and requires a free
[Earthdata Login](https://urs.earthdata.nasa.gov/) account. Set credentials
via environment variables to avoid hardcoding them:

```r
Sys.setenv(EARTHDATA_USER = "your_username")
Sys.setenv(EARTHDATA_PASS = "your_password")
```

## Notes on hotspot analysis performance

`gi_star_raster()` automatically switches between an exact `spdep`-based
computation (accurate, used for ≤ 200,000 valid cells by default) and a
faster `terra::focal`-based circular-window approximation for larger
rasters. Tune the threshold with `max_cells_direct`.

## Testing

```r
# From the package source directory
testthat::test_dir("tests/testthat", package = "urbanheatisland")
```

All 12 unit tests pass against `terra`, `sf`, `spdep`, and `gstat` (tested
against terra 1.7.65, sf 1.0-15, spdep 1.3-1, gstat 2.1-1 on R 4.3.3).

## Status

v0.1.0 — functional core (LST acquisition/processing, interpolation, UHI,
hotspot analysis) with a passing unit test suite. `man/*.Rd` pages are not
yet generated — run `roxygen2::roxygenise()` after editing the roxygen
comments in `R/*.R` to build them, or use `devtools::document()`.
