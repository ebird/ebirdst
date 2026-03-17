# Inventory of downloaded eBird Status and Trends data

Returns a summary of all eBird Status and Trends data packages currently
downloaded to disk, with separate rows for the Status and Trends data
products for each species.

## Usage

``` r
ebirdst_data_inventory(path = ebirdst_data_dir())

# S3 method for class 'ebirdst_inventory'
print(x, ...)
```

## Arguments

- path:

  character; directory where data are stored. Defaults to
  [`ebirdst_data_dir()`](https://ebird.github.io/ebirdst/reference/ebirdst_data_dir.md).

- x:

  an `ebirdst_inventory` object as returned by
  `ebirdst_data_inventory()`.

- ...:

  ignored.

## Value

A tibble of class `ebirdst_inventory` with one row per data package
found on disk, with columns `species_code`, `common_name`,
`scientific_name`, `version_year`, `dataset` ("status" or "trends"),
`n_files`, and `size_mb`. The object has a compact print method that
displays the inventory grouped by version year and dataset.

## Examples

``` r
if (FALSE) { # \dontrun{
# inventory of all downloaded data
ebirdst_data_inventory()

# inventory for a specific directory
ebirdst_data_inventory("/path/to/data")
} # }
```
