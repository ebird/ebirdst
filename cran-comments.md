# ebirdst 3.2023.2

- Add `ebirdst_regional_stats()` to download (on first use) and load regional summary statistics for all species
- Add `ebirdst_data_inventory()` and `ebirdst_delete()` to manage files downloaded by `ebirdst`
- Move to air autoformatting and jarl linting
- Efficiency improvements for grid_sample()

## Test environments

- local MacOS install, R 4.6
- Windows (github actions), R 4.6
- MacOS (github actions), R 4.6
- ubuntu 22.04.1 (github actions), R release, devel, and oldrel-1
- win-builder (devel and release)
- R-hub

## R CMD check results

0 errors | 0 warnings | 1 notes

- NOTE: Version contains large components (3.2023.1). We've aligned our version numbers with the version numbers for the API that this package interacts with. The eBird Status and Trends data products are given a version corresponding to a year, with the current version being 2022, so we've included that year in our version number to indicate that this package only works with the 2023 version of the data.
