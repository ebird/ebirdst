# ebirdst 4.2023.1

- New function `assign_weeks_to_seasons()` identifies which season each of the 52 weeks of the year falls within for a given species, only considering seasons meeting a minimum expert review quality score
- Removed all functions previously listed as deprecated or defunct
  (`abundance_palette()`, `ebirdst_download()`, `ebirdst_extent()`,
  `ebirdst_habitat()`, `ebirdst_ppms()`, `ebirdst_ppms_ts()`,
  `ebirdst_subset()`, `load_pds()`, `load_pis()`, `load_predictions()`,
  `load_stixels()`, `parse_raster_dates()`, `plot_pds()`, `plot_pis()`,
  `project_extent()`, `stixelize()`); they have been unavailable or erroring
  since at least v3.2022.1
- Backend approach to file download has been refactored to an on-demand first approach
- `list_available_pis()` no longer downloads every predictor importance raster to determine availability, only `pi_rangewide.csv`
- The http fallback for VPNs that block https now also applies to file downloads, not just file listings. The fallback is only attempted when https fails to reach the server at all, never when the server responds, so the access key isn't sent over an unencrypted connection unnecessarily
- Errors for data that can't be found on-demand now include function-specific guidance, e.g. pointing to `list_available_pis()`
- Files are now downloaded to a temporary file and only moved into place once the transfer completes. Previously a transfer that was cut short part way left a partial file behind, which was treated as a completed download and never re-downloaded; a forced re-download that failed also deleted the existing local copy of the file
- Downloads that fail for a reason other than the data not being available, e.g. a dropped connection, now raise an error saying so rather than reporting the data as missing
- The access key is no longer included in download error messages. The key is passed to the API in the query string of the request URL, and errors from failed downloads quoted that URL, so users reporting a download problem were inadvertently sharing their private key. Download errors now report the reason for the failure with the key redacted
- `vectorize_trends()` now drops locations with zero relative abundance from the circles output entirely, rather than giving them a missing radius; when all remaining locations share a single abundance quantile bin they now all receive the maximum circle radius
- `ebirdst_palettes()` now requires `n` to be a whole number, rather than accepting a value such as `n = 10.5`
- `ebirdst_regional_stats()` no longer prints a message while downloading
- `get_species_path()` now raises its documented "No data package found" error when the data directory itself hasn't been created yet, e.g. on a fresh install, rather than a cryptic assertion failure
- Minimum required `terra` version bumped to 1.7-3; `rasterize_trends()` no longer carries a compatibility branch for older versions that couldn't rasterize multiple fields in a single call
- A connection-level download failure that persists after the http fallback is now retried a couple more times with a short backoff before the file is given up on, rather than failing immediately on a single transient blip
- `load_data_coverage()`'s arguments have been reordered so that the required `weeks` argument comes before `product`, which now has a default; calls relying on positional matching of `product` first must be updated, e.g. `load_data_coverage("05-10", product = "selection-probability")` rather than `load_data_coverage("selection-probability", weeks = "05-10")`
- `load_pi()` and `load_ppm()` now check for GeoTIFF read support before downloading, consistent with `load_raster()` and `load_data_coverage()`, rather than only discovering the missing GDAL driver after the download completes
- `get_species()` now documents that unmatched input returns `NA`, and identifies example datasets by the generic `-example` suffix rather than a hardcoded `"yebsap-example"` check
- The package now formally opts in to the 3rd edition of testthat (`Config/testthat/edition: 3`); the test suite's remaining `context()` calls and `expect_is()` usages, both deprecated since testthat 3.0.0, have been removed and replaced with `expect_type()`/`expect_s3_class()`/`expect_s4_class()` as appropriate
- `set_ebirdst_access_key()` now warns when a project-level `.Renviron` is found in the working directory, since R gives that file precedence over `~/.Renviron` and would otherwise silently prevent the key from being found in a fresh session
- `grid_sample_stratified()` now raises an informative error when a stratifying column (`sample_by` or `year`) contains missing values, rather than silently dropping those rows
- `grid_sample_stratified()` now raises an informative error when `obs_column` contains missing values and `case_control = TRUE`, rather than failing with a cryptic error
- `assign_to_grid()` now correctly errors when a point falls below, not just above, the bounds of a reused `grid_definition`
- `grid_sample_stratified()` no longer recomputes the spatiotemporal grid assignment for the detection-oversampling pool on every one of its up to 25 resampling iterations; the assignment is now computed once and reused
- `grid_sample()` and `grid_sample_stratified()` now support the spatial-only sampling that the documentation described but that was impossible to actually request: passing a 2-element `res` drops the time dimension of the grid, and the temporal element of `coords` becomes optional. Previously any attempt to do so failed with a validation error
- `grid_sample_stratified()` now passes `coords` and `is_lonlat` through to `grid_sample()` in the case where no strata are defined (`by_year = FALSE`, `case_control = FALSE`, and no `sample_by`); previously both were silently dropped and the defaults were used instead
- `grid_sample_stratified(unified_grid = TRUE)` now defines the unified grid using the `res` and `jitter_grid` values passed via `...`; previously it always used a jittered 3 km by 7 day grid and silently ignored those arguments
- Various small bug fixes and typos

# ebirdst 4.2023.0

- Transition to having all the `load_*()` functions download directly rather than having to call `ebirdst_download_status()`
- Converted vignettes to Quarto and moved them to website-only pkgdown articles; the package no longer ships built-in vignettes to CRAN (documentation lives at <https://ebird.github.io/ebirdst/>)
- Add `ebirdst_regional_stats()` to load regional summary statistics for all species
- Add `ebirdst_data_inventory()` and `ebirdst_delete()` to manage files downloaded by `ebirdst`
- Move to air auto-formatting and jarl linting
- Efficiency improvements for grid_sample()
- `grid_sample_stratified()` gains a `cell_quantile_cap` argument to limit how many observations a single chronically over-sampled site (e.g. a bird feeder) can contribute

# ebirdst 3.2023.1

- added a function to generate abundance-scaled circles for trends
- fixed bug preventing tibbles from being passed to grid sampling functions
- clarified documentation for sampling function
- fixed bug in `get_species()` for Yellow-bellied Sapsucker

# ebirdst 3.2023.0

- update for 2023 data release
- add capability to download and load data coverage layers
- Northern Goshawk species code was incorrect
- on some VPNs downloading from https raises an error, switch to http in these cases
- update vignettes: add links to YouTube, expand applications, add API vignette

# ebirdst 3.2022.3

- arrow is back on CRAN, move from Suggests back to Imports
- add 6 new species for Australia

# ebirdst 3.2022.2

- switch terminology from "trajectory" to "migration chronology"
- ensure `rasterize_trends()` works for older versions of `terra` (issue #7)
- move `arrow` package to Suggests until it is back on CRAN (see https://github.com/apache/arrow/issues/39806)

# ebirdst 3.2022.1

- Documented which functions are deprecated and defunct relative to version
  2.2021.3 under topics `ebirdst-defunct` and `ebirdst-deprecated` and added
  them back into the package.  This allows other packages to conditionally 
  reference them when 2.2021.3 is installed while still passing CRAN checks. 

# ebirdst 3.2022.0

- new 2022 status data and trends data released for the first time!
- major overhaul to allow more targeting downloading of data
- all stixel-level results (PPMS/PIs/PDs) removed, replaced with spatialized raster versions
- no more restart required after updating API key
- change package-level documentation as per roxygen2 suggestions

# ebirdst 2.2021.3

- fix bug causing stixels with missing bounds to raise an error in `ebirdst_habitat()`
- add a function to estimate MCC-F1 for `ebirdst_ppms()`

# ebirdst 2.2021.2

- add a more robust grid sampling function.

# ebirdst 2.2021.1

- release the final batch of 300 species for 2021 bringing the total to 2,282

# ebirdst 2.2021.0

- transition from using raster to terra for handling raster data
- move the following packages from Imports to Suggests: gbm, mgcv, precrec, PresenceAbsence
- move package to the eBird GitHub organization https://github.com/ebird/ebirdst

# ebirdst 1.2021.3

- patch to fix a bug introduced in last release causing missing config files from data downloads [issue #44]

# ebirdst 1.2021.2

- fix bug causing species with same base code to be downloaded together, e.g. leafly also downloads leafly2 [issue #43]

# ebirdst 1.2021.1

- fix bug with extent in `load_fac_map_parameters()`, GitHub issue #40
- use dynamic PAT cutoff in PPM calculations
- update species list to account for second release of eBird data this year

# ebirdst 1.2021.0

- update for the v2021 eBird Status and Trends data

# ebirdst 1.2020.1

- CRAN checks found files created and left behind in ~/Desktop, relocated test files to tempdir() and deleting after test completion with withr::defer()  

# ebirdst 1.2020.0

- major update to align with the new eBird Status and Trends API
- update to align with the 2020 eBird Status Data Products
- transition from rappdirs to tools::R_user_dir() for handling download directories
- all new vignettes

# ebirdst 0.3.5

- bug fix: API update is causing all data downloads to fail

# ebirdst 0.3.4

- rename master branch to main on GitHub requires different download path for example data

# ebirdst 0.3.3

- move example data to GitHub

# ebirdst 0.3.2

- again try to prevent tests and examples from leaving files behind to pass CRAN checks

# ebirdst 0.3.1

- prevent tests and examples from leaving files behind to pass CRAN checks

# ebirdst 0.3.1

- prevent tests and examples from leaving files behind to pass CRAN checks

# ebirdst 0.3.0

- add support for new data structures used for 2020 eBird Status and Trends
- functionality to handle partial dependence data added
- overhaul of package API to be more intuitive and streamlined
- all documentation and vignettes updated

# ebirdst 0.2.2

- add support for variable ensemble support in `compute_ppms()`

# ebirdst 0.2.1

- bug fix: corrected date types in seasonal definitions
- bug fix: fixed possibility that ebirdst_extent could produce invalid date (day 366 of 2015)
- added import of pipe operator
- `velox` was archived, removed dependency from Suggests
- `fasterize` was archived, removed dependency from Imports

# ebirdst 0.2.0

- change maintainer to Matthew Strimas-Mackey
- update to access 2019 status and trends data
- partial dependence data no longer available, all references to PDs removed
- bug fix: `load_raster()` gave incorrect names to seasonal rasters
- bug fix: didn't properly implement quantile binning
- `date_to_st_week()` gets the status and trends week for a give vector of dates

# ebirdst 0.1.0

- first CRAN release
