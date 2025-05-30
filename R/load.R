#' Load eBird Status Data Products raster data
#'
#' Each of the eBird Status raster products is packaged as a GeoTIFF file
#' representing predictions on a regular grid. The core products are occurrence,
#' count, relative abundance, and proportion of population. This function loads
#' one of the available data products into R as a
#' [SpatRaster][terra::SpatRaster] object. Note that data must be downloaded
#' using [ebirdst_download_status()] prior to loading it using this function.
#'
#' @param species character; the species to load data for, given as a scientific
#'   name, common name or six-letter species code (e.g. "woothr"). The full list
#'   of valid species is in the [ebirdst_runs] data frame included in this
#'   package. To download the example dataset, use `"yebsap-example"`.
#' @param product character; eBird Status raster product to load: occurrence,
#'   count, relative abundance, or proportion of population. See Details for a
#'   detailed explanation of each of these products.
#' @param period character; temporal period of the estimation. The eBird Status
#'   models make predictions for each week of the year; however, as a
#'   convenience, data are also provided summarized at the seasonal or annual
#'   ("full-year") level.
#' @param metric character; by default, the weekly products provide estimates of
#'   the median value (`metric = "median"`) and the summarized products are the
#'   cell-wise mean across the weeks within the season (`metric = "mean"`).
#'   However, additional variants exist for some of the products. For the weekly
#'   relative abundance, confidence intervals are provided: specify `metric =
#'   "lower"` to get the 10th quantile or `metric = "upper"` to get the 90th
#'   quantile. For the seasonal and annual products, the cell-wise maximum
#'   values across weeks can be obtained with `metric = "max"`.
#' @param resolution character; the resolution of the raster data to load. The
#'   default is to load the native 3 km resolution data; however, for some
#'   applications 9 km or 27 km data may be suitable.
#' @inheritParams ebirdst_download_status
#'
#' @details The core eBird Status data products provide weekly estimates across
#'   a regular spatial grid. They are packaged as rasters with 52 layers, each
#'   corresponding to estimates for a week of the year, and we refer to them as
#'   "cubes" (e.g. the "relative abundance cube"). All estimates are the median
#'   expected value for a standard 2 km, 1 hour eBird Traveling Count by an
#'   expert eBird observer at the optimal time of day and for optimal weather
#'   conditions to observe the given species. These products are:
#'
#' - `occurrence`: the expected probability (0-1) of occurrence of a species.
#' - `count`: the expected count of a species, conditional on its occurrence at
#' the given location.
#' - `abundance`: the expected relative abundance of a species, computed as the
#' product of the probability of occurrence and the count conditional on
#' occurrence.
#' - `proportion-population`: the proportion of the total relative abundance
#' within each cell. This is a derived product calculated by dividing each cell
#' value in the relative abundance raster by the total abundance summed across
#' all cells.
#'
#' In addition to these weekly data cubes, this function provides access to data
#' summarized over different periods. Seasonal cubes are produced by taking the
#' cell-wise mean or max across the weeks within each season. The boundary dates
#' for each season are species specific and are available in `ebirdst_runs`, and
#' if a season failed review no associated layer will be included in the cube.
#' In addition, full-year summaries provide the mean or max across all weeks of
#' the year that fall within a season that passed review. Note that this is not
#' necessarily all 52 weeks of the year. For example, if the estimates for the
#' non-breeding season failed expert review for a given species, the full-year
#' summary for that species will not include the weeks that would fall within
#' the non-breeding season.
#'
#' @return For the weekly cubes, a [SpatRaster][terra::SpatRaster] with 52
#'   layers for the given product, where the layer names are the dates
#'   (`YYYY-MM-DD` format) of the midpoint of each week. Seasonal cubes will
#'   have up to four layers named with the corresponding season. The full-year
#'   products will have a single layer.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # download example data if hasn't already been downloaded
#' ebirdst_download_status("yebsap-example")
#'
#' # weekly relative abundance
#' # note that only 27 km data are available for the example data
#' abd_weekly <- load_raster("yebsap-example", "abundance", resolution = "27km")
#'
#' # the weeks for each layer are stored in the layer names
#' names(abd_weekly)
#' # they can be converted to date objects with as.Date
#' as.Date(names(abd_weekly))
#'
#' # max seasonal abundance
#' abd_seasonal <- load_raster("yebsap-example", "abundance",
#'                             period = "seasonal", metric = "max",
#'                             resolution = "27km")
#' # available seasons in stack
#' names(abd_seasonal)
#' # subset to just breeding season abundance
#' abd_seasonal[["breeding"]]
#' }
load_raster <- function(species,
                        product = c("abundance",
                                    "count",
                                    "occurrence",
                                    "proportion-population"),
                        period = c("weekly",
                                   "seasonal",
                                   "full-year"),
                        metric = NULL,
                        resolution = c("3km", "9km", "27km"),
                        path = ebirdst_data_dir()) {

  stopifnot(is.character(species), length(species) == 1, dir.exists(path))
  product <- match.arg(product)
  period <- match.arg(period)
  resolution <- match.arg(resolution)

  species_code <- get_species(species)
  species_path <- get_species_path(species, path = path,
                                   dataset = "status",
                                   check_downloaded = FALSE)
  if (!dir.exists(species_path)) {
    stop("No data found for the requested species. Ensure that the data were ",
         "downloaded using ebirdst_download_status() and that the 'path' ",
         "argument correctly points to the data download directory.")
  }

  # check that the geotiff driver is installed
  drv <- terra::gdal(drivers = TRUE)
  drv <- drv$name[stringr::str_detect(drv$can, "read")]
  if (!"GTiff" %in% drv) {
    stop("GDAL does not have GeoTIFF support. GeoTIFF support is required to ",
         "load Status and Trends raster data.")
  }

  # load config file
  p <- load_config(species = species_code, path = path)
  v <- p$srd_pred_year

  # only low res data available for example
  is_example <- stringr::str_detect(species_code, "-example")
  if (is_example && !resolution == "27km") {
    stop("The example data only contains 27 km estimates.")
  }

  # full year products only available for migrants
  if (p$summarize_as_resident && period == "full-year") {
    stop("Full-year products are not available for residents, use ",
         "period = 'seasonal' instead.")
  }

  # construct file name and path
  if (period == "weekly") {
    # assess which metric is being requested
    if (is.null(metric)) {
      metric <- "median"
    }
    if (product == "abundance") {
      if (!metric %in% c("median", "lower", "upper")) {
        stop("Valid metrics for weekly abundance data are 'median', 'lower', ",
             "or 'upper'")
      }
    } else {
      if (metric != "median") {
        stop("For this product, metric must be 'median'")
      }
    }

    # construct filename
    file <- stringr::str_glue("{species_code}_{product}_{metric}",
                              "_{resolution}_{v}.tif")
    file <- file.path(species_path, "weekly", file)
  } else {
    # assess which metric is being requested
    if (is.null(metric)) {
      metric <- "mean"
    }
    if (!metric %in% c("mean", "max")) {
      stop("Valid metrics for seasonal or full-year data are 'mean' or 'max.'")
    }

    # construct filename
    file <- stringr::str_glue("{species_code}_{product}_{period}_{metric}",
                              "_{resolution}_{v}.tif")
    file <- file.path(species_path, "seasonal", file)
  }

  # check existence of target file
  if (!file.exists(file)) {
    stop("The file for the requested product does not exist: \n  ", file)
  }

  # load and return raster stack
  return(terra::rast(file))
}


#' Load eBird Trends estimates for a set of species
#'
#' Load the relative abundance trend estimates for a single species or a set of
#' species. Trends are estimated on a 27 km by 27 km grid for a single season
#' per species (breeding, non-breeding, or resident).  Note that data must be
#' downloaded using [ebirdst_download_trends()] prior to loading it using this
#' function.
#'
#' The trends in relative abundance are estimated using a double machine
#' learning model. To quantify uncertainty, an ensemble of 100 estimates is made
#' at each location, each based on a random subsample of eBird data. The
#' estimated trend is the median across the ensemble, and the 80% confidence
#' intervals are the lower 10th and upper 90th percentiles across the ensemble.
#' To access estimates from the individual folds making up the ensemble use
#' `fold_estimates = TRUE`. These fold-level estimates can be used to quantify
#' uncertainty, for example, when calculating the trend for a given region. For
#' further details on the methodology used to estimate trends consult Fink et
#' al. 2023.
#'
#' @inheritParams ebirdst_download_trends
#' @param fold_estimates logical; by default, the trends summarized across the
#'   100-fold ensemble are returned; however, by setting `fold_estimates = TRUE`
#'   the individual fold-level estimates are returned.
#'
#' @references
#'   Fink, D., Johnston, A., Strimas-Mackey, M., Auer, T., Hochachka, W. M.,
#'   Ligocki, S., Oldham Jaromczyk, L., Robinson, O., Wood, C., Kelling, S., &
#'   Rodewald, A. D. (2023). A Double machine learning trend model for citizen
#'   science data. Methods in Ecology and Evolution, 00, 1–14.
#'   https://doi.org/10.1111/2041-210X.14186
#'
#' @return A data frame containing the trends estimates for a set of species.
#'   The following columns are included:
#'   - `species_code`: the alphanumeric eBird species code uniquely identifying
#'   the species.
#'   - `season`:  season that the trend was estimated for: breeding,
#'   non-breeding, or resident.
#'   - `start_year/end_year`: the start and end years of the trend time period.
#'   - `start_date/end_date`: the start and end dates (`MM-DD` format) of the
#'   season for which the trend was estimated.
#'   - `srd_id`: unique integer identifier for the grid cell.
#'   - `longitude/latitude`: longitude and latitude of the grid cell center.
#'   - `abd`: relative abundance estimate for the middle of the trend time
#'   period (e.g. 2014 for a 2007-2021 trend).
#'   - `abd_ppy`: the median estimated percent per year change in relative
#'   abundance.
#'   - `abd_ppy_lower/abd_ppy_upper`: the 80% confidence interval for the
#'   estimated percent per year change in relative abundance.
#'   - `abd_ppy_nonzero`: a logical (TRUE/FALSE) value indicating if the 80%
#'   confidence limits overlap zero (FALSE) or don't overlap zero (TRUE)
#'   - `abd_trend`: the median estimated cumulative change in relative
#'   abundance over the trend time period.
#'   - `abd_trend_lower/abd_trend_upper`: the 80% confidence interval for the
#'   estimated cumulative change in relative abundance over the trend time
#'   period.
#'
#'   If `fold_estimates = TRUE`, a data frame of fold-level trend estimates is
#'   returned with the following columns:
#'   - `species_code`: the alphanumeric eBird species code uniquely identifying
#'   the species.
#'   - `season`:  season that the trend was estimated for: breeding,
#'   non-breeding, or resident.
#'   - `srd_id`: unique integer identifier for the grid cell.
#'   - `abd`: relative abundance estimate for the middle of the trend time
#'   period (e.g. 2014 for a 2007-2021 trend).
#'   - `abd_ppy`: the estimated percent per year change in relative abundance.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # download example trends data if it hasn't already been downloaded
#' ebirdst_download_trends("yebsap-example")
#'
#' # load trends
#' trends <- load_trends("yebsap-example")
#'
#' # load fold-level estimates
#' trends_folds <- load_trends("yebsap-example", fold_estimates = TRUE)
#' }
load_trends <- function(species,
                        fold_estimates = FALSE,
                        path = ebirdst_data_dir()) {
  stopifnot(is.character(species), !is.na(species), dir.exists(path))
  stopifnot(is_flag(fold_estimates))

  v <- ebirdst_version()[["trends_version_year"]]

  # trends species and seaons
  species_code <- get_species(species)
  trends_runs <- ebirdst::ebirdst_runs[ebirdst::ebirdst_runs$has_trends, ]
  season <- trends_runs$trends_season[match(species_code,
                                            trends_runs$species_code)]
  if (any(is.na(season))) {
    stop("The following species do not have trends estimates:\n  ",
         paste(species[is.na(season)], collapse = ", "))
  }

  # get paths to trends parquet files
  trends_paths <- character()
  for (i in seq_along(species_code)) {
    p <- get_species_path(species_code[i],
                          path = path,
                          dataset = "trends",
                          check_downloaded = FALSE)
    if (fold_estimates) {
      f <- stringr::str_glue("{species_code[i]}_{season[i]}_ebird-trends_",
                             "folds_{v}.parquet")
    } else {
      f <- stringr::str_glue("{species_code[i]}_{season[i]}_ebird-trends_",
                             "{v}.parquet")
    }
    trends_paths <- c(trends_paths, file.path(p, "trends", f))
  }

  if (!all(file.exists(trends_paths))) {
    stop("No trends data found for the following species. Ensure that the ",
         "data were downloaded using ebirdst_download_trends() and that the ",
         "'path' argument correctly points to the data download directory.\n  ",
         paste(species[file.exists(trends_paths)], collapse = ", "))
  }

  # load data
  trends <- NULL
  for (pq in trends_paths) {
    trends <- dplyr::bind_rows(trends, arrow::read_parquet(pq))
  }
  return(trends)
}


#' Load eBird Status and Trends Data Coverage Products
#'
#' The data coverage products are packaged as individual GeoTIFF files for each
#' product for each week of the year. This function loads one of the available
#' data products for one or more weeks into R as a
#' [SpatRaster][terra::SpatRaster] object. Note that data must be downloaded
#' using [ebirdst_download_data_coverage()] prior to loading it using this
#' function.
#'
#' @param product character; data coverage raster product to load: spatial
#'   coverage or site selection probability.
#' @param weeks character; one or more weeks (expressed in `"MM-DD"` format) to
#'   load the raster layers for. If this argument is not specified, all
#'   downloaded weeks will be loaded. **Note that these rasters are quite large
#'   so it's recommended to only load a small number of weeks of data at the
#'   same time.**
#' @inheritParams ebirdst_download_status
#'
#' @details In addition to the species-specific data products, the eBird Status
#'   data products include two products providing estimates of weekly data
#'   coverage at 3 km spatial resolution:
#'
#' - `spatial-coverage`: a spatially smoothed estimate of the proportion of the
#' area that was covered by eBird checklists for the given week.
#' - `selection-probability`: a modeled estimate of the probability that the
#' given location and habitat was sampled by eBird data in the given week.
#'
#' @return A [SpatRaster][terra::SpatRaster] with between 1 and 52 layers for
#'   the given product for the given weeks, where the layer names are the dates
#'   (`YYYY-MM-DD` format) of the midpoint of each week.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # download example data if hasn't already been downloaded
#' ebirdst_download_data_coverage()
#'
#' # load a single week of site selection probability data
#' load_data_coverage("selection-probability", weeks = "01-04")
#'
#' # load all weeks of spatial coverage data
#' load_data_coverage("spatial-coverage", weeks = c("01-04", "01-11"))
#' }
load_data_coverage <- function(product = c("spatial-coverage",
                                           "selection-probability"),
                               weeks = NULL,
                               path = ebirdst_data_dir()) {
  product <- match.arg(product)
  stopifnot(is.character(weeks))
  stopifnot(dir.exists(path))

  dc_path <- get_species_path("data_coverage", path = path,
                              dataset = "status",
                              check_downloaded = FALSE)
  if (!dir.exists(dc_path)) {
    stop("No data coverage products were found. Ensure that the data were ",
         "downloaded using ebirdst_download_data_coverage() and that the ",
         "'path' argument correctly points to the data download directory.")
  }

  # check that the geotiff driver is installed
  drv <- terra::gdal(drivers = TRUE)
  drv <- drv$name[stringr::str_detect(drv$can, "read")]
  if (!"GTiff" %in% drv) {
    stop("GDAL does not have GeoTIFF support. GeoTIFF support is required to ",
         "load Status and Trends raster data.")
  }

  # generate vector of valid weeks
  valid_weeks <- as.Date(paste(2018, seq(4, 366, 7)), format = "%Y %j")
  valid_weeks <- format(valid_weeks, format = "%m-%d")
  if (!is.null(weeks) && !all(weeks %in% valid_weeks)) {
    stop("The following weeks are invalid: ",
         paste(weeks[!weeks %in% valid_weeks], collapse = ", "), "\n",
         "Valid weeks include: ", paste(valid_weeks, collapse = ", "))
  }
  # subset to selected weeks
  if (!is.null(weeks)) {
    valid_weeks <- intersect(valid_weeks, weeks)
  }
  valid_weeks <- paste(ebirdst_version()[["status_version_year"]],
                       valid_weeks, sep = "-")

  # construct filenames
  product <- paste0(product, "_mean")
  files <- stringr::str_glue("{product}_{valid_weeks}.tif")
  files <- file.path(dc_path, product, files)

  # check existence of target files
  if (!all(file.exists(files))) {
    missing_files <- files[!file.exists(files)]
    stop("The files for the requested product do not exist. You may need to ",
         "download them using ebirdst_download_data_coverage(): \n  ",
         paste(basename(files), sep = "\n"))
  }

  # load and return raster stack
  return(stats::setNames(terra::rast(files), valid_weeks))
}


#' Load seasonal eBird Status and Trends range polygons
#'
#' Range polygons are defined as the boundaries of non-zero seasonal relative
#' abundance estimates, which are then (optionally) smoothed to produce more
#' aesthetically pleasing polygons using the `smoothr` package.
#'
#' @inheritParams load_raster
#' @param resolution character; the raster resolution from which the range
#'   polygons were derived.
#' @param smoothed logical; whether smoothed or unsmoothed ranges should be
#'   loaded.
#'
#' @return An `sf` update containing the seasonal range boundaries, with each
#'   season provided as a different feature.
#' @export
#'
#' @examples
#' \dontrun{
#' # download example data if hasn't already been downloaded
#' ebirdst_download_status("yebsap-example")
#'
#' # load smoothed ranges
#' # note that only 27 km data are provided for the example data
#' ranges <- load_ranges("yebsap-example", resolution = "27km")
#' }
load_ranges <- function(species,
                        resolution = c("9km", "27km"), smoothed = TRUE,
                        path = ebirdst_data_dir()) {
  stopifnot(is.character(species), length(species) == 1, dir.exists(path))
  stopifnot(is.logical(smoothed), length(smoothed) == 1)
  resolution <- match.arg(resolution)

  species_code <- get_species(species)
  species_path <- get_species_path(species, path = path,
                                   dataset = "status",
                                   check_downloaded = FALSE)
  if (!dir.exists(species_path)) {
    stop("No data found for the requested species. Ensure that the data were ",
         "downloaded using ebirdst_download_status() and that the 'path' ",
         "argument correctly points to the data download directory.")
  }

  # load config file
  p <- load_config(species = species_code, path = path)
  v <- p$srd_pred_year

  # only low res data available for example
  is_example <- stringr::str_detect(species_code, "-example")
  if (is_example && !resolution == "27km") {
    stop("The example data only contains 27 km estimates.")
  }

  # define filename
  label <- ifelse(smoothed, "smooth", "raw")
  file <- stringr::str_glue("{species_code}_range_{label}",
                            "_{resolution}_{v}.gpkg")
  file <- file.path(species_path, "ranges", file)

  # check existence of target file
  if (!file.exists(file)) {
    stop("The file for the requested product does not exist: \n  ", file)
  }

  # load polygons
  p <- sf::read_sf(dsn = file, layer = "range")

  return(p)
}


#' Load regional summary statistics
#'
#' Load seasonal summary statistics for regions consisting of countries and
#' states/provinces.
#'
#' @inheritParams load_raster
#'
#' @return A data frame containing regional summary statistics with columns:
#'   - `species_code`: alphanumeric eBird species code.
#'   - `region_type`: `country` for countries or `state` for states, provinces,
#'   or other sub-national regions.
#'   - `region_code`: alphanumeric code for the region.
#'   - `region_name`: English name of the region.
#'   - `continent_code`: alphanumeric code for continent that this region
#'   belongs to.
#'   - `continent_name`: name of the continent that this region belongs to.
#'   - `season`: name of the season that the summary statistics were calculated
#'   for.
#'   - `abundance_mean`: mean relative abundance in the region.
#'   - `total_pop_percent`: proportion of the seasonal modeled population
#'   falling within the region.
#'   - `range_percent_occupied`: the proportion of the region occupied by the
#'   species during the given season.
#'   - `range_total_percent`: the proportion of the species seasonal range
#'   falling within the region.
#'   - `range_days_occupation`: number of days of the season that the region was
#'   occupied by this species.
#' @export
#'
#' @examples
#' \dontrun{
#' # download example data if hasn't already been downloaded
#' ebirdst_download_status("yebsap-example")
#'
#' # load configuration parameters
#' regional <- load_regional_stats("yebsap-example")
#' }
load_regional_stats <- function(species, path = ebirdst_data_dir()) {
  stopifnot(is.character(species), length(species) == 1, dir.exists(path))

  species_code <- get_species(species)
  species_path <- get_species_path(species, path = path,
                                   dataset = "status",
                                   check_downloaded = FALSE)
  if (!dir.exists(species_path)) {
    stop("No data found for the requested species. Ensure that the data were ",
         "downloaded using ebirdst_download_status() and that the 'path' ",
         "argument correctly points to the data download directory.")
  }

  file <- file.path(species_path, "regional_stats.csv")
  if(!file.exists(file)) {
    stop("The regional summary stats file could not be found. To download ",
         "file, use `ebirst_download_status(download_regional = TRUE)`.")
  }
  # load stats
  stats <- dplyr::as_tibble(utils::read.csv(file, na = "", row.names = NULL))
  stats[["region_area_km2"]] <- NULL
  return(stats)
}


#' Load eBird Status Data Products configuration file
#'
#' Load the configuration file for an eBird Status run. This configuration file
#' is mostly for internal use and contains a variety of parameters used in the
#' modeling process.
#'
#' @inheritParams load_raster
#'
#' @return A list with the run configuration parameters.
#' @export
#'
#' @examples
#' \dontrun{
#' # download example data if hasn't already been downloaded
#' ebirdst_download_status("yebsap-example")
#'
#' # load configuration parameters
#' p <- load_config("yebsap-example")
#' }
load_config <- function(species, path = ebirdst_data_dir()) {
  stopifnot(is.character(species), length(species) == 1, dir.exists(path))

  species_code <- get_species(species)
  species_path <- get_species_path(species, path = path,
                                   dataset = "status",
                                   check_downloaded = FALSE)
  if (!dir.exists(species_path)) {
    stop("No data found for the requested species. Ensure that the data were ",
         "downloaded using ebirdst_download_status() and that the 'path' ",
         "argument correctly points to the data download directory.")
  }

  cfg_file <- file.path(species_path, "config.json")
  if(!file.exists(cfg_file)) {
    stop("The file 'config.json' does not exist in: ", species_path)
  }
  # load configuration file
  p <- jsonlite::read_json(cfg_file, simplifyVector = TRUE)
  names(p) <- tolower(names(p))
  return(p)
}


#' Load full annual cycle map parameters
#'
#' Get the map parameters used on the eBird Status and Trends website to
#' optimally display the full annual cycle data. This includes bins for the
#' abundance data, a projection, and an extent to map. The extent is the spatial
#' extent of non-zero data across the full annual cycle and the projection is
#' optimized for this extent.
#'
#' @inheritParams load_raster
#'
#' @return A list containing elements:
#' - `custom_projection`: a custom projection optimized for the given species'
#'    full annual cycle
#' - `fa_extent`: a [SpatExtent][terra::ext()] object storing the spatial
#' extent of non-zero
#' data for the given species in the custom projection
#' - `res`: a numeric vector with 2 elements giving the target resolution of
#'    raster in the custom projection
#' - `fa_extent_projected`: the extent in projected (Equal Earth) coordinates
#' - `weekly_bins`/`weekly_labels`: weekly abundance bins and labels for the
#' full annual cycle
#' - `seasonal_bins`/`seasonal_labels: seasonal abundance bins and labels for
#' the full annual cycle
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # download example data if hasn't already been downloaded
#' ebirdst_download_status("yebsap-example")
#'
#' # load configuration parameters
#' load_fac_map_parameters(path)
#' }
load_fac_map_parameters <- function(species, path = ebirdst_data_dir()) {
  stopifnot(is.character(species), length(species) == 1, dir.exists(path))

  species_code <- get_species(species)
  species_path <- get_species_path(species, path = path,
                                   dataset = "status",
                                   check_downloaded = FALSE)
  if (!dir.exists(species_path)) {
    stop("No data found for the requested species. Ensure that the data were ",
         "downloaded using ebirdst_download_status() and that the 'path' ",
         "argument correctly points to the data download directory.")
  }

  # load config file
  p <- load_config(species = species_code, path = path)
  ext_order <- unlist(p$bbox_native)[c("xmin", "xmax", "ymin", "ymax")]

  list(custom_projection = p$projection$crs,
       fa_extent = terra::ext(p$projection$extent),
       res = p$projection$res,
       fa_extent_projected = terra::ext(ext_order),
       weekly_bins = p$bins[["3km"]]$breaks,
       weekly_labels = p$bins[["3km"]]$labels,
       seasonal_bins = p$bins_seasonal[["3km"]]$breaks,
       seasonal_labels = p$bins_seasonal[["3km"]]$labels)
}


#' Load predictor importance (PI) rasters
#'
#' The eBird Status models estimate the relative importance of each of the core
#' environmental predictor used in the model (i.e. the % land and water cover
#' variables). These predictor importance (PI) data are converted to ranks (with
#' a rank of 1 being the most important) relative to the full suite of
#' environmental predictors. The ranks are summarized to a 27 km resolution
#' raster grid for each predictor, where the cell values are the average across
#' all models in the ensemble contributing to that cell. These data are
#' available in raster format provided `download_pis = TRUE` was used when
#' calling [ebirdst_download_status()]. PI estimates are available separately
#' for both the occurrence and count sub-model and only the 30 most important
#' predictors are distributed. Use [list_available_pis()] to see which
#' predictors have PI data.
#'
#' @inheritParams load_raster
#' @param predictor character; the predictor that the PI data should be loaded
#'   for. The list of predictors that PI data are available for varies by
#'   species, use [list_available_pis()] to get the list for a given species.
#' @param response character; the model (occurrence or count) that the PI
#'   data should be loaded for.
#'
#' @return A [SpatRaster][terra::SpatRaster] object with the PI ranks for the
#'   given predictor. For migrants, the estimates are weekly and the raster will
#'   have 52 layers, where the layer names are the dates (`MM-DD` format) of the
#'   midpoint of each week. For residents, a single year round layer is
#'   returned.
#'
#' [list_available_pis()] returns a data frame listing the top 30 predictors for
#' which PI rasters can be loaded. In addition to the predictor names, the mean
#' range-wide rank (`rank_mean`) is given as well as the integer rank
#' (`rank`) relative to the full suite of predictors (environmental and effort).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # download example data if hasn't already been downloaded
#' ebirdst_download_status("yebsap-example", download_pis = TRUE)
#'
#' # identify the top predictor
#' top_preds <- list_available_pis("yebsap-example")
#' print(top_preds[1, ])
#'
#' # load predictor importance raster of top predictor for occurrence
#' load_pi("yebsap-example", top_preds$predictor[1])
#' }
load_pi <- function(species, predictor, response = c("occurrence", "count"),
                    path = ebirdst_data_dir()) {
  stopifnot(is.character(species), length(species) == 1, dir.exists(path))
  response <- match.arg(response)

  species_code <- get_species(species)
  species_path <- get_species_path(species, path = path,
                                   dataset = "status",
                                   check_downloaded = FALSE)
  if (!dir.exists(species_path)) {
    stop("No data found for the requested species. Ensure that the data were ",
         "downloaded using ebirdst_download_status() and that the 'path' ",
         "argument correctly points to the data download directory.")
  }

  # construct file name
  year <- load_config(species = species, path = path)[["srd_pred_year"]]
  p <- stringr::str_replace_all(predictor, "_", "-")
  tif <- stringr::str_glue("{species_code}_pi_{response}_{p}_27km_{year}.tif")
  tif <- file.path(species_path, "pis", tif)
  if (!file.exists(tif)) {
    stop("GeoTIFF for ", predictor, " PI could not be found. To download ",
         "PI data use ebirst_download_status(download_pis = TRUE). To list ",
         "predictors that have PI data use list_available_pis().")
  }
  return(terra::rast(tif))
}


#' @describeIn load_pi list the predictors that have PI information for this
#'   species.
#' @export
list_available_pis <- function(species, path = ebirdst_data_dir()) {
  stopifnot(is.character(species), length(species) == 1, dir.exists(path))

  species_code <- get_species(species)
  species_path <- get_species_path(species, path = path,
                                   check_downloaded = FALSE)
  if (!dir.exists(species_path)) {
    stop("No data found for the requested species. Ensure that the data were ",
         "downloaded using ebirdst_download_status() and that the 'path' ",
         "argument correctly points to the data download directory.")
  }

  csv_file <- file.path(species_path, "pis", "pi_rangewide.csv")
  if(!file.exists(csv_file)) {
    stop("The PI data could not be found. To download, use ",
         "`ebirst_download_status(download_pis = TRUE)`.")
  }
  # load ranks
  ranks <- utils::read.csv(csv_file, row.names = NULL, na = "")

  # available pis
  tifs <- list.files(file.path(species_path, "pis"), pattern = "*.tif")
  tifs <- tifs[!stringr::str_detect(tifs, "n-folds")]
  preds <- stringr::str_remove(tifs, "^[^_]+_pi_(occurrence|count)_")
  preds <- stringr::str_extract(preds, "[-a-z0-9]+")
  preds <- unique(stringr::str_replace_all(preds, "-", "_"))
  preds <- preds[preds %in% ranks$predictor]

  # return ranks
  top_ranks <- dplyr::as_tibble(ranks[ranks$predictor %in% preds, ])
  top_ranks <- dplyr::arrange(top_ranks, .data$response, .data$rank)

  return(top_ranks)
}


#' Load predictive performance metric (PPM) rasters
#'
#' eBird Status models are evaluated against a test set of eBird data not used
#' during model training and a suite of predictive performance metrics (PPMs)
#' are calculated. The PPMs for each base model are summarized to a 27 km
#' resolution raster grid, where the cell values are the average across all
#' models in the ensemble contributing to that cell. These data are available in
#' raster format provided `download_ppms = TRUE` was used when calling
#' [ebirdst_download_status()].
#'
#' @inheritParams load_raster
#' @param ppm character; the name of a single metric to load data for. See
#'   Details for definitions of each metric.
#'
#' @details
#' Nineteen predictive performance metrics are provided:
#' - `binary_f1`: F1-score comparing the model predictions converted to binary
#' with the observed detection/non-detection for the test checklists.
#' - `binary_mcc`: Matthews Correlation Coefficient (MCC) comparing the model
#' predictions converted to binary with the observed detection/non-detection
#' for the test checklists.
#' - `binary_prevalence`: the observed detection probability after
#' spatiotemporal subsampling.
#' - `occ_bernoulli_dev`: proportion of Bernoulli deviance explained comparing
#' the predicted occurrence with the observed detection/non-detection for the
#' test checklists.
#' - `occ_bin_spearman`: test observations are binned by predicted encounter
#' rate with bin widths of 0.05, then the mean observed prevalence and predicted
#' encounter rate are calculated within bins. This metric is the Spearman's rank
#' correlation coefficient comparing the observed and predicted binned mean
#' values.
#' - `occ_brier`: the Brier score is the mean squared difference between
#' predicted encounter rate and observed detection/non-detection.
#' - `occ_pr_auc`: the area on the precision-recall curve (PR AUC) generated by
#' comparing the predicted encounter rate with the observed
#' detection/non-detection for the test checklists.
#' - `occ_pr_auc_gt_prev`: the proportion of the ensemble for which the PR AUC
#' is greater than observed prevalence, which indicates that the model is
#' performing better than random guessing.
#' - `occ_pr_auc_normalized`: the PR AUC normalized to account for class
#' imbalance so that a value of 0 represents performance equal to random
#' guessing and a value of 1 represents perfect classification.
#' - `count_log_pearson`: Pearson correlation coefficient comparing the
#' logarithm of the predicted count with the logarithm of the observed count for
#' the subset of test checklists on which the species was detected.
#' - `count_mae`: the mean absolute error (MAE) comparing the observed and
#' predicted counts for the subset of test checklists on which the species was
#' detected.
#' - `count_poisson_dev`: proportion of Poisson deviance explained, comparing
#' the observed and predicted counts for the subset of test checklists on which
#' the species was detected.
#' - `count_rmse`: route mean squared error (RMSE) comparing the observed and
#' predicted counts for the subset of test checklists on which the species was
#' detected.
#' - `count_spearman`: Spearman's rank correlation coefficient comparing the
#' observed and predicted counts for the subset of test checklists on which the
#' species was detected.
#' - `abd_log_pearson`: Pearson correlation coefficient comparing the logarithm
#' of the predicted relative abundance with the logarithm of the observed
#' count for the full set of test checklists.
#' - `abd_mae`: the mean absolute error (MAE) comparing the observed counts and
#' predicted relative abundance for the full set of test checklists.
#' - `abd_poisson_dev`: proportion of Poisson deviance explained, comparing the
#' predicted relative abundance with the observed count for the full set of test
#' checklists.
#' - `abd_rmse`: root mean squared error comparing the predicted relative
#' abundance with the observed count for the full set of test checklists.
#' - `abd_spearman`: Spearman's rank correlation coefficient comparing the
#' predicted relative abundance with the observed count for the full set of
#' test checklists.
#'
#' @return A [SpatRaster][terra::SpatRaster] object with the PPM data. For
#'   migrants, rasters are weekly with  52 layers, where the layer names are the
#'   dates (`MM-DD` format) of the midpoint of each week. For residents, a
#'   single year round layer is returned.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # download example data if hasn't already been downloaded
#' ebirdst_download_status("yebsap-example", download_ppms = TRUE)
#'
#' # load area under the precision-recall curve PPM raster
#' load_ppm("yebsap-example", ppm = "binary_pr_auc")
#' }
load_ppm <- function(species,
                     ppm = c("binary_f1",
                             "binary_mcc",
                             "binary_prevalence",
                             "occ_bernoulli_dev",
                             "occ_bin_spearman",
                             "occ_brier",
                             "occ_pr_auc",
                             "occ_pr_auc_gt_prev",
                             "occ_pr_auc_normalized",
                             "count_log_pearson",
                             "count_mae",
                             "count_poisson_dev",
                             "count_rmse",
                             "count_spearman",
                             "abd_log_pearson",
                             "abd_mae",
                             "abd_poisson_dev",
                             "abd_rmse",
                             "abd_spearman"),
                     path = ebirdst_data_dir()) {
  stopifnot(is.character(species), length(species) == 1, dir.exists(path))
  ppm <- match.arg(ppm)

  species_code <- get_species(species)
  species_path <- get_species_path(species, path = path,
                                   dataset = "status",
                                   check_downloaded = FALSE)
  if (!dir.exists(species_path)) {
    stop("No data found for the requested species. Ensure that the data were ",
         "downloaded using ebirdst_download_status() and that the 'path' ",
         "argument correctly points to the data download directory.")
  }

  # construct file name
  year <- load_config(species = species, path = path)[["srd_pred_year"]]
  p <- stringr::str_replace_all(ppm, "_", "-")
  tif <- stringr::str_glue("{species_code}_ppm_{p}_mean_27km_{year}.tif")
  tif <- file.path(species_path, "ppms", tif)
  if (!file.exists(tif)) {
    stop("GeoTIFF for ", ppm, " PPM could not be found. To download ",
         "PPM data use ebirst_download_status(download_ppm = TRUE).")
  }
  return(terra::rast(tif))
}
