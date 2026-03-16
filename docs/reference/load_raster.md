# Load eBird Status Data Products raster data

Each of the eBird Status raster products is packaged as a GeoTIFF file
representing predictions on a regular grid. The core products are
occurrence, count, relative abundance, and proportion of population.
This function loads one of the available data products into R as a
[SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
object. Note that data must be downloaded using
[`ebirdst_download_status()`](https://ebird.github.io/ebirdst/reference/ebirdst_download_status.md)
prior to loading it using this function.

## Usage

``` r
load_raster(
  species,
  product = c("abundance", "count", "occurrence", "proportion-population"),
  period = c("weekly", "seasonal", "full-year"),
  metric = NULL,
  resolution = c("3km", "9km", "27km"),
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

- product:

  character; eBird Status raster product to load: occurrence, count,
  relative abundance, or proportion of population. See Details for a
  detailed explanation of each of these products.

- period:

  character; temporal period of the estimation. The eBird Status models
  make predictions for each week of the year; however, as a convenience,
  data are also provided summarized at the seasonal or annual
  ("full-year") level.

- metric:

  character; by default, the weekly products provide estimates of the
  median value (`metric = "median"`) and the summarized products are the
  cell-wise mean across the weeks within the season (`metric = "mean"`).
  However, additional variants exist for some of the products. For the
  weekly relative abundance, confidence intervals are provided: specify
  `metric = "lower"` to get the 10th quantile or `metric = "upper"` to
  get the 90th quantile. For the seasonal and annual products, the
  cell-wise maximum values across weeks can be obtained with
  `metric = "max"`.

- resolution:

  character; the resolution of the raster data to load. The default is
  to load the native 3 km resolution data; however, for some
  applications 9 km or 27 km data may be suitable.

- path:

  character; directory to download the data to. All downloaded files
  will be placed in a sub-directory of this directory named for the data
  version year, e.g. "2020" for the 2020 Status Data Products. Each
  species' data package will then appear in a directory named with the
  eBird species code. Defaults to a persistent data directory, which can
  be found by calling
  [`ebirdst_data_dir()`](https://ebird.github.io/ebirdst/reference/ebirdst_data_dir.md).

## Value

For the weekly cubes, a
[SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
with 52 layers for the given product, where the layer names are the
dates (`YYYY-MM-DD` format) of the midpoint of each week. Seasonal cubes
will have up to four layers named with the corresponding season. The
full-year products will have a single layer.

## Details

The core eBird Status data products provide weekly estimates across a
regular spatial grid. They are packaged as rasters with 52 layers, each
corresponding to estimates for a week of the year, and we refer to them
as "cubes" (e.g. the "relative abundance cube"). All estimates are the
median expected value for a standard 2 km, 1 hour eBird Traveling Count
by an expert eBird observer at the optimal time of day and for optimal
weather conditions to observe the given species. These products are:

- `occurrence`: the expected probability (0-1) of occurrence of a
  species.

- `count`: the expected count of a species, conditional on its
  occurrence at the given location.

- `abundance`: the expected relative abundance of a species, computed as
  the product of the probability of occurrence and the count conditional
  on occurrence.

- `proportion-population`: the proportion of the total relative
  abundance within each cell. This is a derived product calculated by
  dividing each cell value in the relative abundance raster by the total
  abundance summed across all cells.

In addition to these weekly data cubes, this function provides access to
data summarized over different periods. Seasonal cubes are produced by
taking the cell-wise mean or max across the weeks within each season.
The boundary dates for each season are species specific and are
available in `ebirdst_runs`, and if a season failed review no associated
layer will be included in the cube. In addition, full-year summaries
provide the mean or max across all weeks of the year that fall within a
season that passed review. Note that this is not necessarily all 52
weeks of the year. For example, if the estimates for the non-breeding
season failed expert review for a given species, the full-year summary
for that species will not include the weeks that would fall within the
non-breeding season.

## Examples

``` r
if (FALSE) { # \dontrun{
# download example data if hasn't already been downloaded
ebirdst_download_status("yebsap-example")

# weekly relative abundance
# note that only 27 km data are available for the example data
abd_weekly <- load_raster("yebsap-example", "abundance", resolution = "27km")

# the weeks for each layer are stored in the layer names
names(abd_weekly)
# they can be converted to date objects with as.Date
as.Date(names(abd_weekly))

# max seasonal abundance
abd_seasonal <- load_raster("yebsap-example", "abundance",
                            period = "seasonal", metric = "max",
                            resolution = "27km")
# available seasons in stack
names(abd_seasonal)
# subset to just breeding season abundance
abd_seasonal[["breeding"]]
} # }
```
