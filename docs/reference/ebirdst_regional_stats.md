# Regional summary statistics for all species

Load a single file of regional summary statistics covering all species
with eBird Status Data Products. This file is downloaded automatically
on first use and loaded in a single step; subsequent calls load the
already downloaded file directly. This differs from
[`load_regional_stats()`](https://ebird.github.io/ebirdst/reference/load_regional_stats.md),
which loads the regional statistics for a single species from that
species' downloaded data package.

## Usage

``` r
ebirdst_regional_stats(
  path = ebirdst_data_dir(),
  force = FALSE,
  show_progress = interactive()
)
```

## Arguments

- path:

  character; directory that the data are stored in. Defaults to the
  persistent data directory returned by
  [`ebirdst_data_dir()`](https://ebird.github.io/ebirdst/reference/ebirdst_data_dir.md).

- force:

  logical; if the file has already been downloaded, should a fresh copy
  be downloaded anyway.

- show_progress:

  logical; whether to print download progress information. Defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html), so
  downloads are silent in non-interactive sessions (e.g. scripts and R
  Markdown).

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
