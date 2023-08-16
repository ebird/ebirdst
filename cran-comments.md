# ebirdst 2.2021.4

- trends data released!
- performance improvements to `ebirdst_habitat()`

## Test environments

- local MacOS install, R 4.3
- Windows (github actions), R 4.3
- MacOS (github actions), R 4.3
- ubuntu 22.04.1 (github actions), R release, devel, and oldrel-1
- win-builder (devel and release)
- R-hub

## R CMD check results

0 errors | 0 warnings | 2 notes

- NOTE: Version contains large components (2.2021.4). We've aligned our version numbers with the version numbers for the API that this package interacts with. The eBird Status and Trends data products are given a version corresponding to a year, with the current version being 2021, so we've included that year in our version number to indicate that this package only works with the 2021 version of the data.
- NOTE: Found the following (possibly) invalid URLs: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019EA000658. Not sure why this is flagged, the URL is working.
