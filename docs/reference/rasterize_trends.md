# Convert eBird Trends Data Products to raster format

The eBird trends data are stored in a tabular format, where each row
gives the trend estimate for a single cell in a 27 km x 27 km equal area
grid. For many applications, an explicitly spatial format is more
useful. This function uses the cell center coordinates to convert the
tabular trend estimates to raster format in `terra`
[SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
format.

## Usage

``` r
rasterize_trends(
  trends,
  layers = c("abd_ppy", "abd_ppy_lower", "abd_ppy_upper"),
  trim = TRUE
)
```

## Arguments

- trends:

  data frame; trends data for a single species as returned by
  [`load_trends()`](https://ebird.github.io/ebirdst/reference/load_trends.md).

- layers:

  character; column names in the trends data frame to rasterize. These
  columns will become layers in the raster that is created.

- trim:

  logical; flag indicating if the returned raster should be trimmed to
  remove outer rows and columns that are NA. If `trim = FALSE` the
  returned raster will have a global extent, which can be useful if
  rasters will be combined across species with different ranges.

## Value

A
[SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
object.

## Examples

``` r
if (FALSE) { # \dontrun{
# download example trends data if it hasn't already been downloaded
ebirdst_download_trends("yebsap-example")

# load trends
trends <- load_trends("yebsap-example")

# rasterize percent per year trend
rasterize_trends(trends, "abd_ppy")
} # }
```
