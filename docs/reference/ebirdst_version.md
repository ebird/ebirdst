# eBird Status and Trends Data Products version

Identify the version of the eBird Status and Trends Data Products that
this version of the R package works with. Versions are defined by the
year that all model estimates are made for.

## Usage

``` r
ebirdst_version()
```

## Value

A list with two components: for the eBird Status Data Products

- `status_version_year`: version year for the eBird Status Data Products

- `trends_version_year` version year for the eBird Trends Data Products

## Examples

``` r
ebirdst_version()
#> $status_version_year
#> [1] 2025
#> 
#> $trends_version_year
#> [1] 2022
#> 
```
