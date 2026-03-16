# Load eBird Status Data Products configuration file

Load the configuration file for an eBird Status run. This configuration
file is mostly for internal use and contains a variety of parameters
used in the modeling process.

## Usage

``` r
load_config(species, path = ebirdst_data_dir())
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

A list with the run configuration parameters.

## Examples

``` r
if (FALSE) { # \dontrun{
# download example data if hasn't already been downloaded
ebirdst_download_status("yebsap-example")

# load configuration parameters
p <- load_config("yebsap-example")
} # }
```
