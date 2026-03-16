# Load regional summary statistics

Load seasonal summary statistics for regions consisting of countries and
states/provinces.

## Usage

``` r
load_regional_stats(species, path = ebirdst_data_dir())
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

- `range_percent_occupied`: the proportion of the region occupied by the
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
