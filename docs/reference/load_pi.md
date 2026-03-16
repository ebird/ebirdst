# Load predictor importance (PI) rasters

The eBird Status models estimate the relative importance of each of the
core environmental predictor used in the model (i.e. the % land and
water cover variables). These predictor importance (PI) data are
converted to ranks (with a rank of 1 being the most important) relative
to the full suite of environmental predictors. The ranks are summarized
to a 27 km resolution raster grid for each predictor, where the cell
values are the average across all models in the ensemble contributing to
that cell. These data are available in raster format provided
`download_pis = TRUE` was used when calling
[`ebirdst_download_status()`](https://ebird.github.io/ebirdst/reference/ebirdst_download_status.md).
PI estimates are available separately for both the occurrence and count
sub-model and only the 30 most important predictors are distributed. Use
`list_available_pis()` to see which predictors have PI data.

## Usage

``` r
load_pi(
  species,
  predictor,
  response = c("occurrence", "count"),
  path = ebirdst_data_dir()
)

list_available_pis(species, path = ebirdst_data_dir())
```

## Arguments

- species:

  character; the species to load data for, given as a scientific name,
  common name or six-letter species code (e.g. "woothr"). The full list
  of valid species is in the
  [ebirdst_runs](https://ebird.github.io/ebirdst/reference/ebirdst_runs.md)
  data frame included in this package. To download the example dataset,
  use `"yebsap-example"`.

- predictor:

  character; the predictor that the PI data should be loaded for. The
  list of predictors that PI data are available for varies by species,
  use `list_available_pis()` to get the list for a given species.

- response:

  character; the model (occurrence or count) that the PI data should be
  loaded for.

- path:

  character; directory to download the data to. All downloaded files
  will be placed in a sub-directory of this directory named for the data
  version year, e.g. "2020" for the 2020 Status Data Products. Each
  species' data package will then appear in a directory named with the
  eBird species code. Defaults to a persistent data directory, which can
  be found by calling
  [`ebirdst_data_dir()`](https://ebird.github.io/ebirdst/reference/ebirdst_data_dir.md).

## Value

A
[SpatRaster](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
object with the PI ranks for the given predictor. For migrants, the
estimates are weekly and the raster will have 52 layers, where the layer
names are the dates (`MM-DD` format) of the midpoint of each week. For
residents, a single year round layer is returned.

`list_available_pis()` returns a data frame listing the top 30
predictors for which PI rasters can be loaded. In addition to the
predictor names, the mean range-wide rank (`rank_mean`) is given as well
as the integer rank (`rank`) relative to the full suite of predictors
(environmental and effort).

## Functions

- `list_available_pis()`: list the predictors that have PI information
  for this species.

## Examples

``` r
if (FALSE) { # \dontrun{
# download example data if hasn't already been downloaded
ebirdst_download_status("yebsap-example", download_pis = TRUE)

# identify the top predictor
top_preds <- list_available_pis("yebsap-example")
print(top_preds[1, ])

# load predictor importance raster of top predictor for occurrence
load_pi("yebsap-example", top_preds$predictor[1])
} # }
```
