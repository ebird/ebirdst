# Load eBird Trends estimates for a set of species

Load the relative abundance trend estimates for a single species or a
set of species. Trends are estimated on a 27 km by 27 km grid for a
single season per species (breeding, non-breeding, or resident). Note
that data must be downloaded using
[`ebirdst_download_trends()`](https://ebird.github.io/ebirdst/reference/ebirdst_download_trends.md)
prior to loading it using this function.

## Usage

``` r
load_trends(species, fold_estimates = FALSE, path = ebirdst_data_dir())
```

## Arguments

- species:

  character; one or more species given as scientific names, common names
  or six-letter species codes (e.g. "woothr"). The full list of valid
  species can be viewed in the
  [ebirdst_runs](https://ebird.github.io/ebirdst/reference/ebirdst_runs.md)
  data frame included in this package; species with trends estimates are
  indicated by the `has_trends` column. To access the example dataset,
  use "yebsap-example".

- fold_estimates:

  logical; by default, the trends summarized across the 100-fold
  ensemble are returned; however, by setting `fold_estimates = TRUE` the
  individual fold-level estimates are returned.

- path:

  character; directory to download the data to. All downloaded files
  will be placed in a sub-directory of this directory named for the data
  version year, e.g. "2020" for the 2020 Status Data Products. Each
  species' data package will then appear in a directory named with the
  eBird species code. Defaults to a persistent data directory, which can
  be found by calling
  [`ebirdst_data_dir()`](https://ebird.github.io/ebirdst/reference/ebirdst_data_dir.md).

## Value

A data frame containing the trends estimates for a set of species. The
following columns are included:

- `species_code`: the alphanumeric eBird species code uniquely
  identifying the species.

- `season`: season that the trend was estimated for: breeding,
  non-breeding, or resident.

- `start_year/end_year`: the start and end years of the trend time
  period.

- `start_date/end_date`: the start and end dates (`MM-DD` format) of the
  season for which the trend was estimated.

- `srd_id`: unique integer identifier for the grid cell.

- `longitude/latitude`: longitude and latitude of the grid cell center.

- `abd`: relative abundance estimate for the middle of the trend time
  period (e.g. 2014 for a 2007-2021 trend).

- `abd_ppy`: the median estimated percent per year change in relative
  abundance.

- `abd_ppy_lower/abd_ppy_upper`: the 80% confidence interval for the
  estimated percent per year change in relative abundance.

- `abd_ppy_nonzero`: a logical (TRUE/FALSE) value indicating if the 80%
  confidence limits overlap zero (FALSE) or don't overlap zero (TRUE)

- `abd_trend`: the median estimated cumulative change in relative
  abundance over the trend time period.

- `abd_trend_lower/abd_trend_upper`: the 80% confidence interval for the
  estimated cumulative change in relative abundance over the trend time
  period.

If `fold_estimates = TRUE`, a data frame of fold-level trend estimates
is returned with the following columns:

- `species_code`: the alphanumeric eBird species code uniquely
  identifying the species.

- `season`: season that the trend was estimated for: breeding,
  non-breeding, or resident.

- `srd_id`: unique integer identifier for the grid cell.

- `abd`: relative abundance estimate for the middle of the trend time
  period (e.g. 2014 for a 2007-2021 trend).

- `abd_ppy`: the estimated percent per year change in relative
  abundance.

## Details

The trends in relative abundance are estimated using a double machine
learning model. To quantify uncertainty, an ensemble of 100 estimates is
made at each location, each based on a random subsample of eBird data.
The estimated trend is the median across the ensemble, and the 80%
confidence intervals are the lower 10th and upper 90th percentiles
across the ensemble. To access estimates from the individual folds
making up the ensemble use `fold_estimates = TRUE`. These fold-level
estimates can be used to quantify uncertainty, for example, when
calculating the trend for a given region. For further details on the
methodology used to estimate trends consult Fink et al. 2023.

## References

Fink, D., Johnston, A., Strimas-Mackey, M., Auer, T., Hochachka, W. M.,
Ligocki, S., Oldham Jaromczyk, L., Robinson, O., Wood, C., Kelling, S.,
& Rodewald, A. D. (2023). A Double machine learning trend model for
citizen science data. Methods in Ecology and Evolution, 00, 1–14.
https://doi.org/10.1111/2041-210X.14186

## Examples

``` r
if (FALSE) { # \dontrun{
# download example trends data if it hasn't already been downloaded
ebirdst_download_trends("yebsap-example")

# load trends
trends <- load_trends("yebsap-example")

# load fold-level estimates
trends_folds <- load_trends("yebsap-example", fold_estimates = TRUE)
} # }
```
