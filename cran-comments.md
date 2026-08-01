# ebirdst 4.2023.1

- Backend approach to file download has been refactored to an on-demand first approach
- `list_available_pis()` no longer downloads every predictor importance raster to determine availability, only `pi_rangewide.csv`
- The http fallback for VPNs that block https now also applies to file downloads, not just file listings
- Errors for data that can't be found on-demand now include function-specific guidance, e.g. pointing to `list_available_pis()`

## Test environments

- local MacOS install, R 4.6
- Windows (github actions), R 4.6
- MacOS (github actions), R 4.6
- ubuntu 22.04.1 (github actions), R release, devel, and oldrel-1
- win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 notes

- NOTE: Version contains large components (4.2023.1). We've aligned our version numbers with the version numbers for the API that this package interacts with. The eBird Status and Trends data products are given a version corresponding to a year, with the current version being 2022, so we've included that year in our version number to indicate that this package only works with the 2023 version of the data.

## revdepcheck results

We checked 1 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
