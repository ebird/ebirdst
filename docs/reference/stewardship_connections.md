# Shared stewardship connections for migratory birds

Generate a shared stewardship connections raster mapping where migratory
species that occupy an origin region in one stationary season
(breeding/non-breeding) concentrate in the other stationary season
(non-breeding/breeding). For each species, the destination-season
proportion-of-population raster is multiplied by a species weight (the
fraction of the global population in the origin region), and the
weighted rasters are summed across all retained species.

## Usage

``` r
stewardship_connections(
  species = NULL,
  origin_region = NULL,
  region_code = NULL,
  species_weights = NULL,
  origin_season = c("breeding", "nonbreeding"),
  include_migration_seasons = FALSE,
  min_origin_threshold = 1e-04,
  min_migration_threshold = if (include_migration_seasons) 0.5 else 0.85,
  destination_mask = NULL,
  resolution = c("27km", "9km", "3km"),
  filename = NULL,
  path = ebirdst_data_dir(),
  force = FALSE,
  show_progress = interactive(),
  ...
)
```

## Arguments

- species:

  character or `NULL`; one or more species given as scientific names,
  common names, or six-letter species codes. Mixes of name types are
  accepted. The full list of valid species is in
  [ebirdst_runs](https://ebird.github.io/ebirdst/reference/ebirdst_runs.md).
  In custom polygon mode, `species` is required. In region-code mode,
  `species` optionally subsets the candidate list derived from the
  regional stats table. In weights-table mode, `species` optionally
  subsets the rows of `species_weights`.

- origin_region:

  `NULL`, [sf](https://r-spatial.github.io/sf/reference/sf.html), or
  [SpatVector](https://rspatial.github.io/terra/reference/SpatVector-class.html);
  polygon defining the focal origin region. Required for polygon mode;
  must be `NULL` when `region_code` or `species_weights` is supplied. If
  multiple features are present, they are unioned to a single polygon.

- region_code:

  character or `NULL`; a region code matching one present in
  [`ebirdst_regional_stats()`](https://ebird.github.io/ebirdst/reference/ebirdst_regional_stats.md)
  (e.g. `"USA"`, `"USA-NY"`, `"TWN"`). Required for region-code mode;
  must be `NULL` otherwise. Use
  [`ebirdst_regional_stats()`](https://ebird.github.io/ebirdst/reference/ebirdst_regional_stats.md)
  to browse valid codes via the `region_code` column.

- species_weights:

  data frame or `NULL`; a data frame with at least two columns:
  `species_code` (character, resolved via
  [`get_species()`](https://ebird.github.io/ebirdst/reference/get_species.md))
  and `weight` (numeric, finite, non-negative). Required for
  weights-table mode; must be `NULL` otherwise. Rows with duplicate
  resolved `species_code` values are not permitted. Threshold arguments
  (`min_origin_threshold`, `min_migration_threshold`) and
  `include_migration_seasons` are ignored in this mode.

- origin_season:

  character; the stationary season in which the focal region is
  occupied: `"breeding"` (default) or `"nonbreeding"`. The destination
  season is the opposite stationary season.

- include_migration_seasons:

  logical; if `TRUE`, the weight for each species also considers the
  fraction of the population within `origin_region` during the pre- and
  post-breeding migration seasons (taking the maximum across all
  origin-season layers present). The default `min_migration_threshold`
  is relaxed to `0.5` in this case. Applies only to polygon mode and
  region-code mode; ignored in weights-table mode. Defaults to `FALSE`.

- min_origin_threshold:

  numeric `[0, 1)`; minimum weight (fraction of global population in the
  origin region during origin season(s)) required to retain a species.
  Species with `weight < min_origin_threshold` are excluded before the
  destination-side computation. Applies to polygon mode and region-code
  mode; ignored in weights-table mode. Defaults to `0.0001` (1%). Set to
  `0` to disable.

- min_migration_threshold:

  numeric `[0, 1)`; minimum fraction by which a species' population must
  leave the origin region between the origin and destination stationary
  seasons. Applies to polygon mode and region-code mode; ignored in
  weights-table mode. Defaults to `0.85` when
  `include_migration_seasons = FALSE` and `0.5` when `TRUE`. Set to `0`
  to disable.

- destination_mask:

  `NULL`, [sf](https://r-spatial.github.io/sf/reference/sf.html),
  [sf](https://r-spatial.github.io/sf/reference/sfc.html),
  [SpatVector](https://rspatial.github.io/terra/reference/SpatVector-class.html),
  [SpatExtent](https://rspatial.github.io/terra/reference/SpatExtent-class.html),
  or a bounding box from
  [`sf::st_bbox()`](https://r-spatial.github.io/sf/reference/st_bbox.html);
  spatial object used to crop and mask the destination-season raster for
  each species before accumulation. If `NULL`, the full global raster is
  used — a `destination_mask` is strongly recommended for broadly
  distributed species to reduce memory and compute.

- resolution:

  character; resolution of the raster data: `"27km"` (default), `"9km"`,
  or `"3km"`. Note that the example dataset only provides 27 km data. If
  more than 25 species are requested at finer-than-27km resolution, or
  more than 100 at 27km, a message is issued about expected download
  size and compute time.

- filename:

  character or `NULL`; file path to write the connection raster to disk
  via
  [`terra::writeRaster()`](https://rspatial.github.io/terra/reference/writeRaster.html).
  If `NULL` (default), the raster is returned in memory only. When
  `filename` is supplied the return value is invisible and additional
  arguments in `...` are passed to
  [`terra::writeRaster()`](https://rspatial.github.io/terra/reference/writeRaster.html).

- path:

  character; directory to download the data to. All downloaded files
  will be placed in a sub-directory of this directory named for the data
  version year, e.g. "2020" for the 2020 Status Data Products. Each
  species' data package will then appear in a directory named with the
  eBird species code. Defaults to a persistent data directory, which can
  be found by calling
  [`ebirdst_data_dir()`](https://ebird.github.io/ebirdst/reference/ebirdst_data_dir.md).

- force:

  logical; if the data have already been downloaded, should a fresh copy
  be downloaded anyway.

- show_progress:

  logical; whether to print download progress information. Defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html), so
  downloads are silent in non-interactive sessions (e.g. scripts and R
  Markdown).

- ...:

  additional arguments passed to
  [`terra::writeRaster()`](https://rspatial.github.io/terra/reference/writeRaster.html)
  when `filename` is supplied.

## Value

A named list with two elements:

- `connection`: a single-layer
  [SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  giving the stewardship connection surface — the weighted sum of
  destination-season proportion-of-population rasters, with each species
  weighted by the fraction of its global population in the origin
  region.

- `weights`: a data frame with one row per retained species and columns:

  - `species_code`: the alphanumeric eBird species code.

  - `common_name`: the English common name.

  - `origin_season`: which origin-season produced the proportion of
    population `NA` in weights-table mode.

  - `weight`: fraction of global population in the origin region during
    the origin season(s). In region-code mode this is
    `total_pop_percent` from
    [`ebirdst_regional_stats()`](https://ebird.github.io/ebirdst/reference/ebirdst_regional_stats.md);
    in weights-table mode it is the user-supplied value.

  - `rel_pct_decrease`: fraction of the origin-season proportion of
    population that left the region by the destination season. `NA` in
    weights-table mode.

## Details

**Input modes:** The species and weights can be defined in three quite
distinct ways depending on whether `region_code`, `species_weights`, or
`origin_region` is supplied:

- **Region-code mode** (`region_code`): provide a region code (e.g.
  `"USA"`, `"USA-NY"`) matching one present in
  [`ebirdst_regional_stats()`](https://ebird.github.io/ebirdst/reference/ebirdst_regional_stats.md).
  The candidate species list is derived automatically from that region;
  `species` optionally subsets it. Weights and thresholds are derived
  from `total_pop_percent` in the regional stats table — no
  proportion-of-population rasters are downloaded for filtered-out
  species.

- **Custom polygon mode** (`origin_region`): provide an `sf` or
  `SpatVector` polygon defining the focal origin region. In this mode, a
  set of species must be supplied via the `species` argument. For each
  species the weight and inclusion thresholds are computed by
  calculating the seasonal proportion-of-population within the supplied
  polygon.

- **Weights-table mode** (`species_weights`): provide a data frame with
  columns `species_code` and `weight` (numeric, non-negative). Weights
  are used as supplied; threshold arguments and
  `include_migration_seasons` are ignored. `species` optionally subsets
  the rows.

**Species selection note:** There is significant download and compute
cost for processing large numbers of species. For the custom polygon
mode, the typical workflow is to use
[`ebirdst_regional_stats()`](https://ebird.github.io/ebirdst/reference/ebirdst_regional_stats.md),
or other datasets, to coarsely subset the species list before calling
this function. Region-code mode does this subsetting automatically.

This method was developed by the Cornell Lab of Ornithology and Partners
in Flight (PIF). The approach is described in:

- Partners in Flight (2023). Making Connections for Bird Conservation,
  2nd Edition. PIF Technical Series No. 10.
  <https://partnersinflight.org/resources/pif-technical-series-10-making-connections-appendices-and-maps/>

- Lello-Smith et al. (2026). Biological Conservation 313: 111551.
  <https://doi.org/10.1016/j.biocon.2025.111551>

**Important:** The resulting raster shows multispecies abundance
concentrations weighted by origin-region population representation — it
is **not** a migratory connectivity metric. The individuals present in
the destination season are not necessarily the same individuals that
occupy the focal region in the origin season.

## Examples

``` r
if (FALSE) { # \dontrun{
# Appalachian warblers with breeding populations in West Virginia
warblers <- c("cerwar", "woewar1", "louwat", "buwwar", "hoowar")

# West Virginia: rough bounding box in WGS84
origin_region <- sf::st_sf(geometry = sf::st_as_sfc(
  sf::st_bbox(
    c(xmin = -82.6, xmax = -77.7, ymin = 37.2, ymax = 40.6),
    crs = sf::st_crs(4326)
  )
))

# Central America destination mask (Central Mexico through Panama)
destination_mask <- sf::st_bbox(
    c(xmin = -99, xmax = -80, ymin = 10, ymax = 22),
    crs = sf::st_crs(4326)
)

# region-code mode: weights from ebirdst_regional_stats()
result_rc <- stewardship_connections(
  warblers,
  region_code = "USA-WV",
  origin_season = "breeding",
  destination_mask = destination_mask,
  resolution = "27km"
)
result_rc$weights

# polygon mode: weights derived from proportion-population rasters
result_poly <- stewardship_connections(
  warblers,
  origin_region = origin_region,
  origin_season = "breeding",
  destination_mask = destination_mask,
  resolution = "27km"
)
result_poly$weights

# weights-table mode: supply weights directly, skips threshold filtering
wt <- data.frame(
  species_code = warblers,
  weight = c(0.12, 0.08, 0.15, 0.20, 0.10)
)
result_wt <- stewardship_connections(
  species_weights = wt,
  origin_season = "breeding",
  destination_mask = destination_mask,
  resolution = "27km"
)

# display the connection surface
brks <- importance_breaks(result_poly$connection)
pal <- c("#e6e6e6", ebirdst_palettes(length(brks) - 2))
terra::plot(
  result_poly$connection,
  breaks = brks,
  col = pal,
  axes = FALSE,
  box = TRUE,
  legend = FALSE,
  mar = c(0.25, 0.25, 0.25, 0.25)
)
} # }
```
