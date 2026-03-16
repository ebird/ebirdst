# Download eBird Status and Trends Data Coverage Products

In addition to the species-specific data products, the eBird Status data
products include two products providing estimates of weekly data
coverage at 3 km spatial resolution: site selection probability and
spatial coverage. This function downloads these data products in raster
GeoTIFF format.

## Usage

``` r
ebirdst_download_data_coverage(
  path = ebirdst_data_dir(),
  pattern = NULL,
  dry_run = FALSE,
  force = FALSE,
  show_progress = TRUE
)
```

## Arguments

- path:

  character; directory to download the data to. All downloaded files
  will be placed in a sub-directory of this directory named for the data
  version year, e.g. "2020" for the 2020 Status Data Products. Each
  species' data package will then appear in a directory named with the
  eBird species code. Defaults to a persistent data directory, which can
  be found by calling
  [`ebirdst_data_dir()`](https://ebird.github.io/ebirdst/reference/ebirdst_data_dir.md).

- pattern:

  character; regular expression pattern to supply to
  [str_detect()](https://stringr.tidyverse.org/reference/str_detect.html)
  to filter files to download. This filter will be applied in addition
  to any of the `download_` arguments. Note that some files are
  mandatory and will always be downloaded.

- dry_run:

  logical; whether to do a dry run, just listing files that will be
  downloaded. This can be useful when testing the use of `pattern` to
  filter the files to download.

- force:

  logical; if the data have already been downloaded, should a fresh copy
  be downloaded anyway.

- show_progress:

  logical; whether to print download progress information.

## Value

Path to the folder containing the downloaded data coverage products.

## Examples

``` r
if (FALSE) { # \dontrun{
# download all data coverage products
ebirdst_download_data_coverage()

# download just the spatial coverage products
ebirdst_download_data_coverage(pattern = "spatial-coverage")

# download a single week of data coverage products
ebirdst_download_data_coverage(pattern = "01-04")

# download all weeks in april
ebirdst_download_data_coverage(pattern = "04-")
} # }
```
