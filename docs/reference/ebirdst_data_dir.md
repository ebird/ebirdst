# Path to eBird Status and Trends data download directory

Identify and return the path to the default download directory for eBird
Status and Trends data products. This directory can be defined by
setting the environment variable `EBIRDST_DATA_DIR`, otherwise the
directory returned by `tools::R_user_dir("ebirdst", which = "data")`
will be used.

## Usage

``` r
ebirdst_data_dir()
```

## Value

The path to the data download directory.

## Examples

``` r
ebirdst_data_dir()
#> [1] "/Users/mes335/data/ebirdst"
```
