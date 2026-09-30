# ebirdst 4.2025.0

- eBird Status Data Products update for 2025

## Test environments

- local MacOS install, R 4.6
- Windows (github actions), R 4.6
- MacOS (github actions), R 4.6
- ubuntu 22.04.1 (github actions), R release, devel, and oldrel-1
- win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 notes

- NOTE: Version contains large components (4.2025.0). We've aligned our version numbers with the version numbers for the API that this package interacts with. The eBird Status and Trends data products are given a version corresponding to a year, with the current version being 2023, so we've included that year in our version number to indicate that this package only works with the 2023 version of the data.

## revdepcheck results

We checked 1 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
