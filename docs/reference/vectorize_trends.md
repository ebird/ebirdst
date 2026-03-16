# Convert Trends Data Products to points or circles

The eBird trends data are stored in a tabular format, where each row
gives the trend estimate for a single cell in a 27 km x 27 km equal area
grid. For many applications, an explicitly spatial format is more
useful. This function uses the cell center coordinates to convert the
tabular trend estimates to points or circles in
[sf](https://r-spatial.github.io/sf/reference/sf.html) format. Trends
can be converted to points or to circles with areas roughly proportional
to the relative abundance within that 27 km grid cell. These
abundance-scaled circles are what is used to produce the trends maps on
the eBird Status and Trends website.

## Usage

``` r
vectorize_trends(trends, output = c("circles", "points"), crs = 4326)
```

## Arguments

- trends:

  data frame; trends data for a single species as returned by
  [`load_trends()`](https://ebird.github.io/ebirdst/reference/load_trends.md).

- output:

  character; "points" outputs spatial points while "circles" outputs
  circles with areas roughly proportional to the relative abundance
  within that 27 km grid cell.

- crs:

  character or `sf`
  [crs](https://r-spatial.github.io/sf/reference/st_crs.html) object;
  coordinate reference system to output the results in. For points,
  unprojected latitude-longitude coordinates (the default) are most
  typical, while for circles use whatever equal area CRS you intend to
  use when mapping the data otherwise the "circles" will appear skewed.

## Value

Vetorized trends data as an
[sf](https://r-spatial.github.io/sf/reference/sf.html) object.

## Examples

``` r
if (FALSE) { # \dontrun{
# download example trends data if it hasn't already been downloaded
ebirdst_download_trends("yebsap-example")

# load trends
trends <- load_trends("yebsap-example")

# vectorize as points
vectorize_trends(trends, "points")
# vectorize as circles
vectorize_trends(trends, "circles", crs = "+proj=eqearth")
} # }
```
