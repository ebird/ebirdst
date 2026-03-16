# Download eBird Trends Data Products

Download eBird Trends Data Products for set of species, or for an
example species. Downloading Status and Trends data requires an access
key, consult
[`set_ebirdst_access_key()`](https://ebird.github.io/ebirdst/reference/set_ebirdst_access_key.md)
for instructions on how to obtain and store this key. The example data
consist of the results for Yellow-bellied Sapsucker subset to Michigan
and are much smaller than the full dataset, making these data quicker to
download and process. The example data are accessible without an access
key.

## Usage

``` r
ebirdst_download_trends(
  species,
  path = ebirdst_data_dir(),
  force = FALSE,
  show_progress = TRUE
)
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

  logical; whether to print download progress information.

## Value

Character vector of paths to the folders containing the downloaded data
packages for the given species. The trends data will be in the `trends/`
subdirectory.

## Examples

``` r
if (FALSE) { # \dontrun{
# download the example data
ebirdst_download_trends("yebsap-example")

# download the data package for wood thrush
ebirdst_download_trends("woothr")

# multiple species can be downloaded at once
ebirdst_download_trends(c("Sage Thrasher", "Abert's Towhee"))
} # }
```
