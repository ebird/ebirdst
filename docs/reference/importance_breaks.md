# Compute importance breaks for a stewardship connections raster

Derives quantile-based bin breaks for displaying a stewardship
connections raster returned by
[`stewardship_connections()`](https://ebird.github.io/ebirdst/reference/stewardship_connections.md).
The lower tail of the distribution (cells at or below the
`lower_quantile` of non-zero values) is collapsed to a single
low-importance class. The upper portion is divided into `n_bins`
equal-quantile bins.

## Usage

``` r
importance_breaks(x, lower_quantile = 0.5, n_bins = 10, return_raster = FALSE)
```

## Arguments

- x:

  [SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html);
  a single-layer stewardship connection raster, typically the
  `$connection` element returned by
  [`stewardship_connections()`](https://ebird.github.io/ebirdst/reference/stewardship_connections.md).

- lower_quantile:

  numeric; quantile (in `[0, 1)`) of the non-zero cell values used as
  the lower cutoff. Cells below this quantile are pooled into a
  low-importance class (returned as `NA` when `return_raster = TRUE`).
  Defaults to `0.5` (the median of non-zero cells).

- n_bins:

  integer; number of display bins for cells at or above the
  `lower_quantile` cutoff. Defaults to `10`.

- return_raster:

  logical; if `FALSE` (default), return the numeric break values. If
  `TRUE`, return a reclassified
  [SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  with integer cell values in `1:n_bins`; cells below the cutoff are
  `NA`.

## Value

When `return_raster = FALSE` (default), a numeric vector of \`n_bins

- 2`break values:`0`, the `lower_quantile`cutoff, and`n_bins`upper quantile boundaries. When`return_raster
  =
  TRUE`, a single-layer [SpatRaster][terra::SpatRaster] with integer cell values in `1:n_bins` (cells below the cutoff are`NA`), ready to plot with `col
  = ebirdst_palettes(n_bins)\`.

## Details

Colors for the bins can be obtained from
[`ebirdst_palettes()`](https://ebird.github.io/ebirdst/reference/ebirdst_palettes.md)
with `n = n_bins`. The function deliberately does not return a palette
so that colors remain under caller control.

## Examples

``` r
# create a synthetic connection raster
r <- terra::rast(nrows = 10, ncols = 10, vals = runif(100, 0, 0.5))

# default: return breaks (length = n_bins + 2 = 12)
brks <- importance_breaks(r)
length(brks)
#> [1] 12

# return a classified raster for direct plotting
r_class <- importance_breaks(r, return_raster = TRUE)
terra::plot(r_class, col = ebirdst_palettes(10))
```
