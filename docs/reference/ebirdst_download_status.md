# Download eBird Status Data Products

Download eBird Status Data Products for a single species, or for an
example species. Downloading Status and Trends data requires an access
key, consult
[`set_ebirdst_access_key()`](https://ebird.github.io/ebirdst/reference/set_ebirdst_access_key.md)
for instructions on how to obtain and store this key. The example data
consist of the results for Yellow-bellied Sapsucker subset to Michigan
and are much smaller than the full dataset, making these data quicker to
download and process. Only the low resolution (27 km) data are available
for the example data. In addition, the example data are accessible
without an access key.

## Usage

``` r
ebirdst_download_status(
  species,
  path = ebirdst_data_dir(),
  download_abundance = TRUE,
  download_occurrence = FALSE,
  download_count = FALSE,
  download_ranges = FALSE,
  download_regional = FALSE,
  download_pis = FALSE,
  download_ppms = FALSE,
  download_all = FALSE,
  pattern = NULL,
  dry_run = FALSE,
  force = FALSE,
  show_progress = TRUE
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

- download_abundance:

  whether to download estimates of abundance and proportion of
  population.

- download_occurrence:

  logical; whether to download estimates of occurrence.

- download_count:

  logical; whether to download estimates of count.

- download_ranges:

  logical; whether to download the range polygons.

- download_regional:

  logical; whether to download the regional summary stats, e.g. percent
  of population in regions.

- download_pis:

  logical; whether to download spatial estimates of predictor
  importance.

- download_ppms:

  logical; whether to download spatial predictive performance metrics.

- download_all:

  logical; download all files in the data package. Equivalent to setting
  all the `download_` arguments to `TRUE`.

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

Path to the folder containing the downloaded data package for the given
species. If `dry_run = TRUE` a list of files to download will be
returned.

## Details

The complete data package for each species contains a large number of
files, all of which are cataloged in the vignettes. Most users will only
require a small subset of these files, so by default this function only
downloads the most commonly used files: GeoTIFFs providing estimate of
relative abundance and proportion of population. For those interested in
additional data products, the arguments starting with `download_`
control the download of these other products. The `pattern` argument
provides even finer grained control over what gets downloaded.

## Examples

``` r
if (FALSE) { # \dontrun{
# download the example data
ebirdst_download_status("yebsap-example")

# download the data package for wood thrush
ebirdst_download_status("woothr")

# use pattern to only download low resolution (27 km) geotiff data
# dry_run can be used to see what files will be downloaded
ebirdst_download_status("lobcur", pattern = "_27km_", dry_run = TRUE)
# use pattern to only download high resolution (3 km) weekly abundance data
ebirdst_download_status("lobcur", pattern = "abundance_median_3km",
                        dry_run = TRUE)
} # }
```
