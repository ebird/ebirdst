# Assign points to a spacetime grid

Given a set of points in space and (optionally) time, define a regular
grid with given dimensions, and return the grid cell index for each
point.

## Usage

``` r
assign_to_grid(
  points,
  coords = NULL,
  is_lonlat = FALSE,
  res,
  jitter_grid = TRUE,
  grid_definition = NULL
)
```

## Arguments

- points:

  data frame; points with spatial coordinates `x` and `y`, and an
  optional time coordinate `t`.

- coords:

  character; names of the spatial and temporal coordinates in the input
  dataframe. Only provide these names if you want to overwrite the
  default coordinate names: `c("x", "y", "t")` or
  `c("longitude", "latitude", "t")` if `is_lonlat = TRUE`.

- is_lonlat:

  logical; if the points are in unprojected, lon-lat coordinates. In
  this case, the input data frame should have columns `"longitude"` and
  `"latitude"` and the points will be projected to an equal area Eckert
  IV CRS prior to grid assignment.

- res:

  numeric; resolution of the grid in the `x`, `y`, and `t` dimensions,
  respectively. If only 2 dimensions are provided, a space only grid
  will be generated. The units of `res` are the same as the coordinates
  in the input data unless `is_lonlat` is true in which case the `x` and
  `y` resolution should be provided in meters.

- jitter_grid:

  logical; whether to jitter the location of the origin of the grid to
  introduce some randomness.

- grid_definition:

  list; object defining the grid via the `origin` and `resolution`
  components. To assign multiple sets of points to exactly the same
  grid, `assign_to_grid()` returns a data frame with a `grid_definition`
  attribute that can be passed to subsequent calls to
  `assign_to_grid()`. `res` and `jitter` are ignored if
  `grid_definition` is provided.

## Value

Data frame with the indices of the space-only and spacetime grid cells.
This data frame will have a `grid_definition` attribute that can be used
to reconstruct the grid.

## Examples

``` r
set.seed(1)

# generate some example points
points_xyt <- data.frame(x = runif(100), y = runif(100), t = rnorm(100))
# assign to grid
cells <- assign_to_grid(points_xyt, res = c(0.1, 0.1, 0.5))

# assign a second set of points to the same grid
assign_to_grid(points_xyt, grid_definition = attr(cells, "grid_definition"))
#> # A tibble: 100 × 2
#>    cell_xy cell_xyt
#>    <chr>   <chr>   
#>  1 4-7     4-7-4   
#>  2 5-4     5-4-5   
#>  3 7-3     7-3-3   
#>  4 10-10   10-10-6 
#>  5 3-7     3-7-4   
#>  6 10-3    10-3-9  
#>  7 10-2    10-2-7  
#>  8 8-5     8-5-7   
#>  9 7-10    7-10-6  
#> 10 2-7     2-7-9   
#> # ℹ 90 more rows

# assign lon-lat points to a 10km space-only grid
points_ll <- data.frame(longitude = runif(100, min = -180, max = 180),
                        latitude = runif(100, min = -90, max = 90))
assign_to_grid(points_ll, res = c(10000, 10000), is_lonlat = TRUE)
#> # A tibble: 100 × 1
#>    cell_xy  
#>    <chr>    
#>  1 2960-1224
#>  2 3184-781 
#>  3 2110-1687
#>  4 1254-617 
#>  5 2407-1571
#>  6 244-1415 
#>  7 3172-924 
#>  8 2894-1604
#>  9 1203-769 
#> 10 2118-1   
#> # ℹ 90 more rows

# overwrite default coordinate names, 5km by 1 week grid
points_names <- data.frame(lon = runif(100, min = -180, max = 180),
                           lat = runif(100, min = -90, max = 90),
                           day = sample.int(365, size = 100))
assign_to_grid(points_names,
               res = c(5000, 5000, 7),
               coords = c("lon", "lat", "day"),
               is_lonlat = TRUE)
#> # A tibble: 100 × 2
#>    cell_xy   cell_xyt    
#>    <chr>     <chr>       
#>  1 5348-68   5348-68-49  
#>  2 2294-1332 2294-1332-40
#>  3 2577-1839 2577-1839-16
#>  4 5159-3343 5159-3343-26
#>  5 867-2655  867-2655-5  
#>  6 5944-2704 5944-2704-19
#>  7 2254-1551 2254-1551-41
#>  8 3453-166  3453-166-51 
#>  9 3515-2926 3515-2926-9 
#> 10 4736-1401 4736-1401-33
#> # ℹ 90 more rows
```
