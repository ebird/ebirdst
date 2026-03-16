# Load full annual cycle map parameters

Get the map parameters used on the eBird Status and Trends website to
optimally display the full annual cycle data. This includes bins for the
abundance data, a projection, and an extent to map. The extent is the
spatial extent of non-zero data across the full annual cycle and the
projection is optimized for this extent.

## Usage

``` r
load_fac_map_parameters(species, path = ebirdst_data_dir())
```

## Arguments

- species:

  character; the species to load data for, given as a scientific name,
  common name or six-letter species code (e.g. "woothr"). The full list
  of valid species is in the
  [ebirdst_runs](https://ebird.github.io/ebirdst/reference/ebirdst_runs.md)
  data frame included in this package. To download the example dataset,
  use `"yebsap-example"`.

- path:

  character; directory to download the data to. All downloaded files
  will be placed in a sub-directory of this directory named for the data
  version year, e.g. "2020" for the 2020 Status Data Products. Each
  species' data package will then appear in a directory named with the
  eBird species code. Defaults to a persistent data directory, which can
  be found by calling
  [`ebirdst_data_dir()`](https://ebird.github.io/ebirdst/reference/ebirdst_data_dir.md).

## Value

A list containing elements:

- `custom_projection`: a custom projection optimized for the given
  species' full annual cycle

- `fa_extent`: a
  [SpatExtent](https://rspatial.github.io/terra/reference/ext.html)
  object storing the spatial extent of non-zero data for the given
  species in the custom projection

- `res`: a numeric vector with 2 elements giving the target resolution
  of raster in the custom projection

- `fa_extent_projected`: the extent in projected (Equal Earth)
  coordinates

- `weekly_bins`/`weekly_labels`: weekly abundance bins and labels for
  the full annual cycle

- `seasonal_bins`/\`seasonal_labels: seasonal abundance bins and labels
  for the full annual cycle

## Examples

``` r
if (FALSE) { # \dontrun{
# download example data if hasn't already been downloaded
ebirdst_download_status("yebsap-example")

# load configuration parameters
load_fac_map_parameters(path)
} # }
```
