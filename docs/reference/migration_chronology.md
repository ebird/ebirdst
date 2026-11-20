# Compute weekly proportion of population within regions

Given one or more species and a region (defined as an
[sf](https://r-spatial.github.io/sf/reference/sf.html) or
[SpatVector](https://rspatial.github.io/terra/reference/SpatVector-class.html)
polygon object), compute the proportion of each species' global modeled
population that falls within each region polygon for each of the 52
weeks of the year. The proportion-of-population raster data are
downloaded automatically on first use.

## Usage

``` r
migration_chronology(
  species,
  region,
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

  [sf](https://r-spatial.github.io/sf/reference/sf.html) or
  [SpatVector](https://rspatial.github.io/terra/reference/SpatVector-class.html);
  polygon features defining the region(s) of interest. Each row
  (feature) is treated as a separate region.

- resolution:

  character; resolution of the raster data: `"27km"` (default), `"9km"`,
  or `"3km"`. Note that the example dataset only provides 27 km data. If
  9 km or 3 km data are requested for more than 25 species that require
  downloading, a warning is issued; in an interactive session
  confirmation is required before the downloads begin.

- bind_regions:

  logical; if `TRUE`, all non-geometry attribute columns from `region`
  are appended to the output. Defaults to `FALSE`.

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
  population within the region for that week. Computed as a weighted sum
  of proportion-of-population raster cells, where boundary cells are
  weighted by the fraction of the cell area covered by the polygon.

- `.region_number`: integer (1-indexed) indicating which polygon in
  `region` the row corresponds to. Only present when `nrow(region) > 1`.

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
