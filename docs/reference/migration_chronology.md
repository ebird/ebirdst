# Compute weekly proportion of population within regions

Given one or more species and a region (defined as an
[sf](https://r-spatial.github.io/sf/reference/sf.html),
[sfc](https://r-spatial.github.io/sf/reference/sfc.html), or
[SpatVector](https://rspatial.github.io/terra/reference/SpatVector-class.html)
polygon object), compute the proportion of each species' global modeled
population (or the population within `reference_region`, if supplied)
that falls within each region polygon for each of the 52 weeks of the
year. The proportion-of-population raster data are downloaded
automatically on first use.

## Usage

``` r
migration_chronology(
  species,
  region,
  reference_region = NULL,
  resolution = c("27km", "9km", "3km"),
  bind_regions = FALSE,
  path = ebirdst_data_dir(),
  force = FALSE,
  show_progress = interactive()
)
```

## Arguments

- species:

  character; one or more species given as scientific names, common
  names, or six-letter species codes. Mixes of name types are accepted.
  The full list of valid species is in
  [ebirdst_runs](https://ebird.github.io/ebirdst/reference/ebirdst_runs.md).

- region:

  [sf](https://r-spatial.github.io/sf/reference/sf.html),
  [sfc](https://r-spatial.github.io/sf/reference/sfc.html), or
  [SpatVector](https://rspatial.github.io/terra/reference/SpatVector-class.html);
  polygon features defining the region(s) of interest. Each feature is
  treated as a separate region.

- reference_region:

  `NULL`, [sf](https://r-spatial.github.io/sf/reference/sf.html),
  [sfc](https://r-spatial.github.io/sf/reference/sfc.html), or
  [SpatVector](https://rspatial.github.io/terra/reference/SpatVector-class.html);
  optional polygon defining the reference population relative to which
  proportions are calculated. If multiple features are present, they are
  unioned into a single polygon. When supplied, proportions are
  calculated for both `region` and `reference_region`, and
  `proportion_population` is computed as the proportion in `region`
  divided by the proportion in `reference_region` (returning `0` when
  both are `0`, and `NA` when only the reference proportion is `0`).
  This can be used, for example, to calculate continental proportions,
  such as the proportion of the North American population within the
  given `region`. Defaults to `NULL` (proportions are relative to the
  global modeled population).

- resolution:

  character; resolution of the raster data: `"27km"` (default), `"9km"`,
  or `"3km"`. Note that the example dataset only provides 27 km data. If
  9 km or 3 km data are requested for more than 25 species that require
  downloading, a warning is issued; in an interactive session
  confirmation is required before the downloads begin.

- bind_regions:

  logical; if `TRUE`, all non-geometry attribute columns from `region`
  are appended to the output (ignored if `region` is an
  [sfc](https://r-spatial.github.io/sf/reference/sfc.html) object, which
  has no attributes). Defaults to `FALSE`.

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

## Value

A data frame with one row per combination of species, week, and region
polygon. Columns:

- `species_code`: the alphanumeric eBird species code.

- `season`: season the week falls within (`"breeding"`, `"nonbreeding"`,
  `"prebreeding_migration"`, `"postbreeding_migration"`, or
  `"resident"`), or `NA` if the week does not fall within a named season
  that passed expert review.

- `week`: midpoint date (`YYYY-MM-DD`) of the week.

- `proportion_population`: proportion of the species' global modeled
  population (or the population within `reference_region`, if supplied)
  within the region for that week. Computed as a weighted sum of
  proportion-of-population raster cells, where boundary cells are
  weighted by the fraction of the cell area covered by the polygon.

- `proportion_population_reference`: proportion of the species' global
  modeled population within `reference_region` for that week. Only
  present when `reference_region` is supplied.

- `.region_number`: integer (1-indexed) indicating which polygon in
  `region` the row corresponds to. Only present when `region` has more
  than one feature.

- Additional columns from `region` if `bind_regions = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
# download example data
ebirdst_download_status("yebsap-example")

# build a region polygon from the raster extent
rast <- load_raster("yebsap-example", "abundance", resolution = "27km")
region <- sf::st_as_sf(
  terra::as.polygons(terra::ext(rast), crs = terra::crs(rast))
)

chron <- migration_chronology("yebsap-example", region, resolution = "27km")
head(chron)
} # }
```
