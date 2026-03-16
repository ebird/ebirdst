# Load eBird Status and Trends Data Coverage Products

The data coverage products are packaged as individual GeoTIFF files for
each product for each week of the year. This function loads one of the
available data products for one or more weeks into R as a
[SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
object. Note that data must be downloaded using
[`ebirdst_download_data_coverage()`](https://ebird.github.io/ebirdst/reference/ebirdst_download_data_coverage.md)
prior to loading it using this function.

## Usage

``` r
load_data_coverage(
  product = c("spatial-coverage", "selection-probability"),
  weeks,
  path = ebirdst_data_dir()
)
```

## Arguments

- product:

  character; data coverage raster product to load: spatial coverage or
  site selection probability.

- weeks:

  character; one or more weeks (expressed in `"MM-DD"` format) to load
  the raster layers for. If this argument is not specified, all
  downloaded weeks will be loaded. **Note that these rasters are quite
  large so it's recommended to only load a small number of weeks of data
  at the same time.**

- path:

  character; directory to download the data to. All downloaded files
  will be placed in a sub-directory of this directory named for the data
  version year, e.g. "2020" for the 2020 Status Data Products. Each
  species' data package will then appear in a directory named with the
  eBird species code. Defaults to a persistent data directory, which can
  be found by calling
  [`ebirdst_data_dir()`](https://ebird.github.io/ebirdst/reference/ebirdst_data_dir.md).

## Value

A
[SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
with between 1 and 52 layers for the given product for the given weeks,
where the layer names are the dates (`YYYY-MM-DD` format) of the
midpoint of each week.

## Details

In addition to the species-specific data products, the eBird Status data
products include two products providing estimates of weekly data
coverage at 3 km spatial resolution:

- `spatial-coverage`: a spatially smoothed estimate of the proportion of
  the area that was covered by eBird checklists for the given week.

- `selection-probability`: a modeled estimate of the probability that
  the given location and habitat was sampled by eBird data in the given
  week.

## Examples

``` r
if (FALSE) { # \dontrun{
# download example data if hasn't already been downloaded
ebirdst_download_data_coverage()

# load a single week of site selection probability data
load_data_coverage("selection-probability", weeks = "01-04")

# load all weeks of spatial coverage data
load_data_coverage("spatial-coverage", weeks = c("01-04", "01-11"))
} # }
```
