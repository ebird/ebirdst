# ebirdst 2.2021.0

- transition from using raster to terra for handling raster data
- move the following packages from Imports to Suggests: gbm, mgcv, precrec, PresenceAbsence
- move package to the eBird GitHub organization https://github.com/ebird/ebirdst

## Test environments

- local MacOS install, R 4.2
- Windows (github actions), R 4.2
- MacOS (github actions), R 4.2
- ubuntu 22.04.1 (github actions), R release, devel, and oldrel-1
- win-builder (devel and release)
- R-hub

## R CMD check results

0 errors | 0 warnings | 2 notes

- NOTE: Version contains large components (2.2021.0). We've aligned our version numbers with the version numbers for the API that this package interacts with. The eBird Status and Trends data products are given a version corresponding to a year, with the current version being 2021, so we've included that year in our version number to indicate that this package only works with the 2021 version of the data.
- NOTE: win-builder highlights the following two links as "Service Unavailable"; however, I've repeatedly checked them and they're both working.
  - https://doi.org/10.1029/2019EA000658
  - https://doi.org/10.1002/eap.2056
