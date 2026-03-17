# Delete downloaded eBird Status and Trends data

Deletes downloaded eBird Status and Trends data packages for specified
species and/or version years. When called interactively without
`force = TRUE`, prints a summary of the data to be deleted and prompts
for confirmation before proceeding.

## Usage

``` r
ebirdst_delete(
  species = NULL,
  year = NULL,
  path = ebirdst_data_dir(),
  force = FALSE
)
```

## Arguments

- species:

  character; one or more species given as eBird species codes,
  scientific names, or English common names. If NULL (the default), data
  for all species are included.

- year:

  integer; one or more version years. If NULL (the default), data for
  all years are included.

- path:

  character; directory where data are stored. Defaults to
  [`ebirdst_data_dir()`](https://ebird.github.io/ebirdst/reference/ebirdst_data_dir.md).

- force:

  logical; if TRUE, skip the interactive confirmation prompt and delete
  without asking. Required when running in a non-interactive session.

## Value

Invisibly returns a character vector of the paths of the deleted
directories.

## Examples

``` r
if (FALSE) { # \dontrun{
# review and confirm deletion of example data
ebirdst_delete(species = "yebsap-example")

# delete all data for a given version year without prompting
ebirdst_delete(year = 2021, force = TRUE)

# delete a specific species and year
ebirdst_delete(species = "Yellow-bellied Sapsucker", year = 2022,
               force = TRUE)
} # }
```
