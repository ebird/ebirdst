# ebirdst 4.2023.0

- Transition to having all the `load_*()` functions download directly rather than having to call `ebirdst_download_status()`
- Converted vignettes to Quarto and moved them to website-only pkgdown articles; the package no longer ships built-in vignettes to CRAN (documentation lives at <https://ebird.github.io/ebirdst/>)
- Add `ebirdst_regional_stats()` to load regional summary statistics for all species
- Add `ebirdst_data_inventory()` and `ebirdst_delete()` to manage files downloaded by `ebirdst`
- Move to air auto-formatting and jarl linting
- Efficiency improvements for `grid_sample()`
- `grid_sample_stratified()` gains a `cell_quantile_cap` argument to limit how many observations a single chronically over-sampled site (e.g. a bird feeder) can contribute

## Test environments

- local MacOS install, R 4.6
- Windows (github actions), R 4.6
- MacOS (github actions), R 4.6
- ubuntu 22.04.1 (github actions), R release, devel, and oldrel-1
- win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 notes

- NOTE: Version contains large components (4.2023.0). We've aligned our version numbers with the version numbers for the API that this package interacts with. The eBird Status and Trends data products are given a version corresponding to a year, with the current version being 2022, so we've included that year in our version number to indicate that this package only works with the 2023 version of the data.

## revdepcheck results

We checked 1 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
