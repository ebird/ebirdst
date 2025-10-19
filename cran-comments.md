# ebirdst 3.2023.1

- added a function to generate abundance-scaled circles for trends
- fixed bug preventing tibbles from being passed to grid sampling functions
- clarified documentation for sampling function
- fixed bug in `get_species()` for Yellow-bellied Sapsucker (issue #11)

## Test environments

- local MacOS install, R 4.5
- Windows (github actions), R 4.5
- MacOS (github actions), R 4.5
- ubuntu 22.04.1 (github actions), R release, devel, and oldrel-1
- win-builder (devel and release)
- R-hub

## R CMD check results

0 errors | 0 warnings | 1 notes

- NOTE: Version contains large components (3.2023.1). We've aligned our version numbers with the version numbers for the API that this package interacts with. The eBird Status and Trends data products are given a version corresponding to a year, with the current version being 2022, so we've included that year in our version number to indicate that this package only works with the 2023 version of the data.
