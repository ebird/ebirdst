# Changelog

## ebirdst 4.2023.1

- New function
  [`assign_weeks_to_seasons()`](https://ebird.github.io/ebirdst/reference/assign_weeks_to_seasons.md)
  identifies which season each of the 52 weeks of the year falls within
  for a given species, only considering seasons meeting a minimum expert
  review quality score; set `return_df = TRUE` to get a data frame with
  one row per week instead of a character vector
- Removed all functions previously listed as deprecated or defunct
  (`abundance_palette()`, `ebirdst_download()`, `ebirdst_extent()`,
  `ebirdst_habitat()`, `ebirdst_ppms()`, `ebirdst_ppms_ts()`,
  `ebirdst_subset()`, `load_pds()`, `load_pis()`, `load_predictions()`,
  `load_stixels()`, `parse_raster_dates()`, `plot_pds()`, `plot_pis()`,
  `project_extent()`, `stixelize()`); they have been unavailable or
  erroring since at least v3.2022.1
- Backend approach to file download has been refactored to an on-demand
  first approach
- [`list_available_pis()`](https://ebird.github.io/ebirdst/reference/load_pi.md)
  no longer downloads every predictor importance raster to determine
  availability, only `pi_rangewide.csv`
- The http fallback for VPNs that block https now also applies to file
  downloads, not just file listings. The fallback is only attempted when
  https fails to reach the server at all, never when the server
  responds, so the access key isn’t sent over an unencrypted connection
  unnecessarily
- Errors for data that can’t be found on-demand now include
  function-specific guidance, e.g. pointing to
  [`list_available_pis()`](https://ebird.github.io/ebirdst/reference/load_pi.md)
- Files are now downloaded to a temporary file and only moved into place
  once the transfer completes. Previously a transfer that was cut short
  part way left a partial file behind, which was treated as a completed
  download and never re-downloaded; a forced re-download that failed
  also deleted the existing local copy of the file
- Downloads that fail for a reason other than the data not being
  available, e.g. a dropped connection, now raise an error saying so
  rather than reporting the data as missing
- The access key is no longer included in download error messages. The
  key is passed to the API in the query string of the request URL, and
  errors from failed downloads quoted that URL, so users reporting a
  download problem were inadvertently sharing their private key.
  Download errors now report the reason for the failure with the key
  redacted
- [`vectorize_trends()`](https://ebird.github.io/ebirdst/reference/vectorize_trends.md)
  now drops locations with zero relative abundance from the circles
  output entirely, rather than giving them a missing radius; when all
  remaining locations share a single abundance quantile bin they now all
  receive the maximum circle radius
- [`ebirdst_palettes()`](https://ebird.github.io/ebirdst/reference/ebirdst_palettes.md)
  now requires `n` to be a whole number, rather than accepting a value
  such as `n = 10.5`
- [`ebirdst_regional_stats()`](https://ebird.github.io/ebirdst/reference/ebirdst_regional_stats.md)
  no longer prints a message while downloading
- [`get_species_path()`](https://ebird.github.io/ebirdst/reference/get_species_path.md)
  now raises its documented “No data package found” error when the data
  directory itself hasn’t been created yet, e.g. on a fresh install,
  rather than a cryptic assertion failure
- Minimum required `terra` version bumped to 1.7-3;
  [`rasterize_trends()`](https://ebird.github.io/ebirdst/reference/rasterize_trends.md)
  no longer carries a compatibility branch for older versions that
  couldn’t rasterize multiple fields in a single call
- A connection-level download failure that persists after the http
  fallback is now retried a couple more times with a short backoff
  before the file is given up on, rather than failing immediately on a
  single transient blip
- [`load_data_coverage()`](https://ebird.github.io/ebirdst/reference/load_data_coverage.md)’s
  arguments have been reordered so that the required `weeks` argument
  comes before `product`, which now has a default; calls relying on
  positional matching of `product` first must be updated,
  e.g. `load_data_coverage("05-10", product = "selection-probability")`
  rather than
  `load_data_coverage("selection-probability", weeks = "05-10")`
- [`load_pi()`](https://ebird.github.io/ebirdst/reference/load_pi.md)
  and
  [`load_ppm()`](https://ebird.github.io/ebirdst/reference/load_ppm.md)
  now check for GeoTIFF read support before downloading, consistent with
  [`load_raster()`](https://ebird.github.io/ebirdst/reference/load_raster.md)
  and
  [`load_data_coverage()`](https://ebird.github.io/ebirdst/reference/load_data_coverage.md),
  rather than only discovering the missing GDAL driver after the
  download completes
- [`get_species()`](https://ebird.github.io/ebirdst/reference/get_species.md)
  now documents that unmatched input returns `NA`, and identifies
  example datasets by the generic `-example` suffix rather than a
  hardcoded `"yebsap-example"` check
- The package now formally opts in to the 3rd edition of testthat
  (`Config/testthat/edition: 3`); the test suite’s remaining `context()`
  calls and `expect_is()` usages, both deprecated since testthat 3.0.0,
  have been removed and replaced with
  `expect_type()`/`expect_s3_class()`/`expect_s4_class()` as appropriate
- [`set_ebirdst_access_key()`](https://ebird.github.io/ebirdst/reference/set_ebirdst_access_key.md)
  now warns when a project-level `.Renviron` is found in the working
  directory, since R gives that file precedence over `~/.Renviron` and
  would otherwise silently prevent the key from being found in a fresh
  session
- Various small bug fixes and typos

## ebirdst 4.2023.0

CRAN release: 2026-07-20

- Transition to having all the `load_*()` functions download directly
  rather than having to call
  [`ebirdst_download_status()`](https://ebird.github.io/ebirdst/reference/ebirdst_download_status.md)
- Converted vignettes to Quarto and moved them to website-only pkgdown
  articles; the package no longer ships built-in vignettes to CRAN
  (documentation lives at <https://ebird.github.io/ebirdst/>)
- Add
  [`ebirdst_regional_stats()`](https://ebird.github.io/ebirdst/reference/ebirdst_regional_stats.md)
  to load regional summary statistics for all species
- Add
  [`ebirdst_data_inventory()`](https://ebird.github.io/ebirdst/reference/ebirdst_data_inventory.md)
  and
  [`ebirdst_delete()`](https://ebird.github.io/ebirdst/reference/ebirdst_delete.md)
  to manage files downloaded by `ebirdst`
- Move to air auto-formatting and jarl linting
- Efficiency improvements for grid_sample()
- [`grid_sample_stratified()`](https://ebird.github.io/ebirdst/reference/grid_sample.md)
  gains a `cell_quantile_cap` argument to limit how many observations a
  single chronically over-sampled site (e.g. a bird feeder) can
  contribute

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
  `ebirdst_habitat()`
- add a function to estimate MCC-F1 for `ebirdst_ppms()`

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
