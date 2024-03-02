# ebirdst 3.2022.3

- arrow is back on CRAN, move from Suggests back to Imports
- add 6 new species for Australia

## Test environments

- local MacOS install, R 4.3
- Windows (github actions), R 4.3
- MacOS (github actions), R 4.3
- ubuntu 22.04.1 (github actions), R release, devel, and oldrel-1
- win-builder (devel and release)
- R-hub

## R CMD check results

0 errors | 0 warnings | 2 notes

- NOTE: Version contains large components (3.2022.3). We've aligned our version numbers with the version numbers for the API that this package interacts with. The eBird Status and Trends data products are given a version corresponding to a year, with the current version being 2022, so we've included that year in our version number to indicate that this package only works with the 2022 version of the data.
