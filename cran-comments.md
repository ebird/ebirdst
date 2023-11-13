# ebirdst 3.2022.0

- new 2022 status data and trends data released for the first time!
- major overhaul to allow more targeting downloading of data
- all stixel-level results (PPMS/PIs/PDs) removed, replaced with spatialized raster versions
- no more restart required after updating API key
- change package-level documentation as per roxygen2 suggestions

## Test environments

- local MacOS install, R 4.3
- Windows (github actions), R 4.3
- MacOS (github actions), R 4.3
- ubuntu 22.04.1 (github actions), R release, devel, and oldrel-1
- win-builder (devel and release)
- R-hub

## R CMD check results

0 errors | 0 warnings | 2 notes

- NOTE: Version contains large components (3.2022.0). We've aligned our version numbers with the version numbers for the API that this package interacts with. The eBird Status and Trends data products are given a version corresponding to a year, with the current version being 2022, so we've included that year in our version number to indicate that this package only works with the 2022 version of the data.
