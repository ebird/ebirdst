# Changelog

## ebirdst 3.2023.1

CRAN release: 2025-10-19

- added a function to generate abundance-scaled circles for trends
- fixed bug preventing tibbles from being passed to grid sampling
  functions
- clarified documentation for sampling function
- fixed bug in
  [`get_species()`](https://ebird.github.io/ebirdst/reference/get_species.md)
  for Yellow-bellied Sapsucker

## ebirdst 3.2023.0

CRAN release: 2025-05-07

- update for 2023 data release
- add capability to download and load data coverage layers
- Northern Goshawk species code was incorrect
- on some VPNs downloading from https raises an error, switch to http in
  these cases
- update vignettes: add links to YouTube, expand applications, add API
  vignette

## ebirdst 3.2022.3

CRAN release: 2024-03-05

- arrow is back on CRAN, move from Suggests back to Imports
- add 6 new species for Australia

## ebirdst 3.2022.2

CRAN release: 2024-02-23

- switch terminology from “trajectory” to “migration chronology”
- ensure
  [`rasterize_trends()`](https://ebird.github.io/ebirdst/reference/rasterize_trends.md)
  works for older versions of `terra` (issue
  [\#7](https://github.com/ebird/ebirdst/issues/7))
- move `arrow` package to Suggests until it is back on CRAN (see
  <https://github.com/apache/arrow/issues/39806>)

## ebirdst 3.2022.1

CRAN release: 2023-12-08

- Documented which functions are deprecated and defunct relative to
  version 2.2021.3 under topics `ebirdst-defunct` and
  `ebirdst-deprecated` and added them back into the package. This allows
  other packages to conditionally reference them when 2.2021.3 is
  installed while still passing CRAN checks.

## ebirdst 3.2022.0

CRAN release: 2023-11-15

- new 2022 status data and trends data released for the first time!
- major overhaul to allow more targeting downloading of data
- all stixel-level results (PPMS/PIs/PDs) removed, replaced with
  spatialized raster versions
- no more restart required after updating API key
- change package-level documentation as per roxygen2 suggestions

## ebirdst 2.2021.3

CRAN release: 2023-05-09

- fix bug causing stixels with missing bounds to raise an error in
  [`ebirdst_habitat()`](https://ebird.github.io/ebirdst/reference/ebirdst-defunct.md)
- add a function to estimate MCC-F1 for
  [`ebirdst_ppms()`](https://ebird.github.io/ebirdst/reference/ebirdst-defunct.md)

## ebirdst 2.2021.2

CRAN release: 2023-04-27

- add a more robust grid sampling function.

## ebirdst 2.2021.1

CRAN release: 2023-04-06

- release the final batch of 300 species for 2021 bringing the total to
  2,282

## ebirdst 2.2021.0

CRAN release: 2023-01-18

- transition from using raster to terra for handling raster data
- move the following packages from Imports to Suggests: gbm, mgcv,
  precrec, PresenceAbsence
- move package to the eBird GitHub organization
  <https://github.com/ebird/ebirdst>

## ebirdst 1.2021.3

CRAN release: 2023-01-11

- patch to fix a bug introduced in last release causing missing config
  files from data downloads \[issue
  [\#44](https://github.com/ebird/ebirdst/issues/44)\]

## ebirdst 1.2021.2

CRAN release: 2023-01-06

- fix bug causing species with same base code to be downloaded together,
  e.g. leafly also downloads leafly2 \[issue
  [\#43](https://github.com/ebird/ebirdst/issues/43)\]

## ebirdst 1.2021.1

CRAN release: 2022-12-07

- fix bug with extent in
  [`load_fac_map_parameters()`](https://ebird.github.io/ebirdst/reference/load_fac_map_parameters.md),
  GitHub issue [\#40](https://github.com/ebird/ebirdst/issues/40)
- use dynamic PAT cutoff in PPM calculations
- update species list to account for second release of eBird data this
  year

## ebirdst 1.2021.0

CRAN release: 2022-11-09

- update for the v2021 eBird Status and Trends data

## ebirdst 1.2020.1

CRAN release: 2022-07-08

- CRAN checks found files created and left behind in ~/Desktop,
  relocated test files to tempdir() and deleting after test completion
  with withr::defer()

## ebirdst 1.2020.0

CRAN release: 2022-07-07

- major update to align with the new eBird Status and Trends API
- update to align with the 2020 eBird Status Data Products
- transition from rappdirs to tools::R_user_dir() for handling download
  directories
- all new vignettes

## ebirdst 0.3.5

CRAN release: 2022-04-01

- bug fix: API update is causing all data downloads to fail

## ebirdst 0.3.4

CRAN release: 2022-03-16

- rename master branch to main on GitHub requires different download
  path for example data

## ebirdst 0.3.3

CRAN release: 2021-11-12

- move example data to GitHub

## ebirdst 0.3.2

CRAN release: 2021-09-15

- again try to prevent tests and examples from leaving files behind to
  pass CRAN checks

## ebirdst 0.3.1

CRAN release: 2021-08-18

- prevent tests and examples from leaving files behind to pass CRAN
  checks

## ebirdst 0.3.1

CRAN release: 2021-08-18

- prevent tests and examples from leaving files behind to pass CRAN
  checks

## ebirdst 0.3.0

CRAN release: 2021-08-10

- add support for new data structures used for 2020 eBird Status and
  Trends
- functionality to handle partial dependence data added
- overhaul of package API to be more intuitive and streamlined
- all documentation and vignettes updated

## ebirdst 0.2.2

CRAN release: 2021-01-16

- add support for variable ensemble support in `compute_ppms()`

## ebirdst 0.2.1

CRAN release: 2020-03-23

- bug fix: corrected date types in seasonal definitions
- bug fix: fixed possibility that ebirdst_extent could produce invalid
  date (day 366 of 2015)
- added import of pipe operator
- `velox` was archived, removed dependency from Suggests
- `fasterize` was archived, removed dependency from Imports

## ebirdst 0.2.0

CRAN release: 2020-02-26

- change maintainer to Matthew Strimas-Mackey
- update to access 2019 status and trends data
- partial dependence data no longer available, all references to PDs
  removed
- bug fix:
  [`load_raster()`](https://ebird.github.io/ebirdst/reference/load_raster.md)
  gave incorrect names to seasonal rasters
- bug fix: didn’t properly implement quantile binning
- [`date_to_st_week()`](https://ebird.github.io/ebirdst/reference/date_to_st_week.md)
  gets the status and trends week for a give vector of dates

## ebirdst 0.1.0

CRAN release: 2019-04-04

- first CRAN release
