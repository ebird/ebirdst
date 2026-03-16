# Get the path to the data package for a given species

This helper function can be used to get the path to a data package for a
given species.

## Usage

``` r
get_species_path(
  species,
  path = ebirdst_data_dir(),
  dataset = c("status", "trends"),
  check_downloaded = TRUE
)
```

## Arguments

- species:

  character; a single species given as a scientific name, common name or
  six-letter species code (e.g. "woothr"). The full list of valid
  species is in the
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

- dataset:

  character; whether the path to the Status or Trends data products
  should be returned.

- check_downloaded:

  logical; raise an error if no data have been downloaded for this
  species.

## Value

The path to the data package directory.

## Examples

``` r
if (FALSE) { # \dontrun{
# get the path
path <- get_species_path("yebsap-example")

# get the path to the full data package for yellow-bellied sapsucker
# common name, scientific name, or species code can be used
path <- get_species_path("Yellow-bellied Sapsucker")
path <- get_species_path("Sphyrapicus varius")
path <- get_species_path("yebsap")
} # }
```
