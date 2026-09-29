# Assign the weeks of the year to seasons

The eBird Status Data Products provide estimates for each of the 52
weeks of the year. For migratory species, the full annual cycle is
divided into four seasons: breeding, non-breeding, pre-breeding
migration, and post-breeding migration; non-migratory species have a
single resident season. The start and end dates of these seasons are
species specific and, in addition, each season is assigned a quality
score from 0 (failed) to 3 (high quality) reflecting how much
extrapolation or omission occurs in that season's estimates. This
function identifies which season each week of the year falls within,
considering only those seasons meeting a minimum quality score. It's
intended to be used to identify the subset of weeks with sufficiently
reliable estimates for a given species, for example prior to summarizing
the weekly data products across the full annual cycle.

## Usage

``` r
assign_weeks_to_seasons(
  species,
  min_quality = 1,
  return_df = FALSE,
  path = ebirdst_data_dir(),
  force = FALSE,
  show_progress = interactive()
)
```

## Arguments

- species:

  character; the species to load data for, given as a scientific name,
  common name or six-letter species code (e.g. "woothr"). The full list
  of valid species is in the
  [ebirdst_runs](https://ebird.github.io/ebirdst/reference/ebirdst_runs.md)
  data frame included in this package. To download the example dataset,
  use `"yebsap-example"`.

- min_quality:

  integer; the minimum quality score (from 1 to 3) that a season must
  have for its weeks to be assigned to it. Weeks falling within a season
  with a lower quality score, or falling outside any season, are
  assigned `NA`.

- return_df:

  logical; if `TRUE`, return a data frame with one row per week and
  columns `week` (date), `season` (character), `quality` (integer, `0`
  for weeks falling outside any season), and `include` (logical, `TRUE`
  if the week's season quality is at least `min_quality`), rather than
  the default character vector.

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

By default, a character vector with 52 elements giving the season that
each week of the year falls within. The elements are in the same order
as the weekly layers of the data products, so this vector can be used
directly to subset the layers of a weekly raster cube. Weeks that don't
fall within a season meeting the minimum quality score are assigned
`NA`. If `return_df = TRUE`, a data frame with one row per week and
columns `week`, `season`, `quality`, and `include` is returned instead.

## Examples

``` r
if (FALSE) { # \dontrun{
# download example data if hasn't already been downloaded
ebirdst_download_status("yebsap-example")

# only weeks in seasons with the highest quality score
seasons <- assign_weeks_to_seasons("yebsap-example", min_quality = 3)

# use these weeks to subset a weekly raster cube
abd <- load_raster("yebsap-example", "abundance", resolution = "27km")
abd_high_quality <- abd[[!is.na(seasons)]]

# return a data frame instead
seasons_df <- assign_weeks_to_seasons(
  "yebsap-example",
  min_quality = 3,
  return_df = TRUE
)
} # }
```
