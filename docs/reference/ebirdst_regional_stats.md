# Regional summary statistics for all species

Load a single file of regional summary statistics covering all species
with eBird Status Data Products. Unlike most data products, which are
downloaded with a dedicated `ebirdst_download_*()` function and then
loaded with a `load_*()` function, this function downloads the file on
first use and loads it in a single step. If the file is not already
present in the data directory it will be downloaded, prompting for
confirmation first. This differs from
[`load_regional_stats()`](https://ebird.github.io/ebirdst/reference/load_regional_stats.md),
which loads the regional statistics for a single species from that
species' downloaded data package.

## Usage

``` r
ebirdst_regional_stats(
  path = ebirdst_data_dir(),
  force = FALSE,
  show_progress = TRUE
)
```

## Arguments

- path:

  character; directory that the data are stored in. Defaults to the
  persistent data directory returned by
  [`ebirdst_data_dir()`](https://ebird.github.io/ebirdst/reference/ebirdst_data_dir.md).

- force:

  logical; if the file has already been downloaded, should a fresh copy
  be downloaded anyway. Setting `force = TRUE` also skips the
  confirmation prompt, which is required to download in a
  non-interactive session.

- show_progress:

  logical; whether to print download progress information.

## Value

A data frame of regional summary statistics for all species. The columns
match those returned by
[`load_regional_stats()`](https://ebird.github.io/ebirdst/reference/load_regional_stats.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# download (if necessary) and load regional stats for all species
regional <- ebirdst_regional_stats()
} # }
```
