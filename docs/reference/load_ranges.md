# Load seasonal eBird Status and Trends range polygons

Range polygons are defined as the boundaries of non-zero seasonal
relative abundance estimates, which are then (optionally) smoothed to
produce more aesthetically pleasing polygons using the `smoothr`
package.

## Usage

``` r
load_ranges(
  species,
  resolution = c("9km", "27km"),
  smoothed = TRUE,
  path = ebirdst_data_dir()
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

- resolution:

  character; the raster resolution from which the range polygons were
  derived.

- smoothed:

  logical; whether smoothed or unsmoothed ranges should be loaded.

- path:

  character; directory to download the data to. All downloaded files
  will be placed in a sub-directory of this directory named for the data
  version year, e.g. "2020" for the 2020 Status Data Products. Each
  species' data package will then appear in a directory named with the
  eBird species code. Defaults to a persistent data directory, which can
  be found by calling
  [`ebirdst_data_dir()`](https://ebird.github.io/ebirdst/reference/ebirdst_data_dir.md).

## Value

An `sf` update containing the seasonal range boundaries, with each
season provided as a different feature.

## Examples

``` r
if (FALSE) { # \dontrun{
# download example data if hasn't already been downloaded
ebirdst_download_status("yebsap-example")

# load smoothed ranges
# note that only 27 km data are provided for the example data
ranges <- load_ranges("yebsap-example", resolution = "27km")
} # }
```
