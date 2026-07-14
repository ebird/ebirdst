# Load regional summary statistics

Load seasonal summary statistics for regions consisting of countries and
states/provinces.

## Usage

``` r
load_regional_stats(
  species,
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

A data frame containing regional summary statistics with columns:

- `species_code`: alphanumeric eBird species code.

- `region_type`: `country` for countries or `state` for states,
  provinces, or other sub-national regions.

- `region_code`: alphanumeric code for the region.

- `region_name`: English name of the region.

- `continent_code`: alphanumeric code for continent that this region
  belongs to.

- `continent_name`: name of the continent that this region belongs to.

- `season`: name of the season that the summary statistics were
  calculated for.

- `abundance_mean`: mean relative abundance in the region.

- `total_pop_percent`: proportion of the seasonal modeled population
  falling within the region.

- `continent_pop_percent`: proportion of the seasonal modeled population
  for the continent (identified by `continent_name`) falling within the
  region.

- `max_week`: the week of the year with the highest proportion of the
  modeled population falling within the region.

- `max_week_percent_pop`: the proportion of the modeled population
  falling within the region in `max_week`, i.e. the maximum weekly
  value.

- `range_occupied_percent`: the proportion of the region occupied by the
  species during the given season.

- `range_total_percent`: the proportion of the species seasonal range
  falling within the region.

- `range_days_occupation`: number of days of the season that the region
  was occupied by this species.

## Examples

``` r
if (FALSE) { # \dontrun{
# download example data if hasn't already been downloaded
ebirdst_download_status("yebsap-example")

# load configuration parameters
regional <- load_regional_stats("yebsap-example")
} # }
```
