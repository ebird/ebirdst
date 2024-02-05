# ebirdst 3.2022.2

- switch terminology from "trajectory" to "migration chronology"
- ensure `rasterize_trends()` works for older versions of `terra` (issue #7)
- move `arrow` package to Suggests until it is back on CRAN (see https://github.com/apache/arrow/issues/39806)

## Test environments

- local MacOS install, R 4.3
- Windows (github actions), R 4.3
- MacOS (github actions), R 4.3
- ubuntu 22.04.1 (github actions), R release, devel, and oldrel-1
- win-builder (devel and release)
- R-hub

## R CMD check results

0 errors | 0 warnings | 2 notes

- NOTE: Version contains large components (3.2022.2). We've aligned our version numbers with the version numbers for the API that this package interacts with. The eBird Status and Trends data products are given a version corresponding to a year, with the current version being 2022, so we've included that year in our version number to indicate that this package only works with the 2022 version of the data.
