# ebirdst 4.2023.1

- New function `assign_weeks_to_seasons()` identifies which season each of the 52 weeks of the year falls within for a given species; set `return_df = TRUE` to get a data frame with one row per week instead of a character vector
- Removed all functions previously listed as deprecated or defunct (they have been unavailable or erroring since at least v3.2022.1)
- `load_data_coverage()`'s arguments have been reordered so that the required `weeks` argument comes before `product`, which now has a default; calls relying on positional matching of `product` first must be updated
- Downloads are more robust and secure: an on-demand-first backend approach, better handling of interrupted/failed transfers and flaky connections, and the access key is no longer exposed in download error messages
- `grid_sample()` and `grid_sample_stratified()` now support space-only sampling (no time dimension), which the documentation described but was previously impossible to request, plus several other sampling bug fixes
- Miscellaneous bug fixes and efficiency improvements

## Test environments

- local MacOS install, R 4.6
- Windows (github actions), R 4.6
- MacOS (github actions), R 4.6
- ubuntu 22.04.1 (github actions), R release, devel, and oldrel-1
- win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 notes

- NOTE: Version contains large components (4.2023.1). We've aligned our version numbers with the version numbers for the API that this package interacts with. The eBird Status and Trends data products are given a version corresponding to a year, with the current version being 2023, so we've included that year in our version number to indicate that this package only works with the 2023 version of the data.

## revdepcheck results

We checked 1 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
