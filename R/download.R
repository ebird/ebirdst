#' Download eBird Status Data Products
#'
#' Download eBird Status Data Products for a single species, or for an example
#' species. Downloading Status and Trends data requires an access key, consult
#' [set_ebirdst_access_key()] for instructions on how to obtain and store this
#' key. The example data consist of the results for Yellow-bellied Sapsucker
#' subset to Michigan and are much smaller than the full dataset, making these
#' data quicker to download and process. Only the low resolution (27 km) data
#' are available for the example data. In addition, the example data are
#' accessible without an access key.
#'
#' @param species character; a single species given as a scientific name, common
#'   name or six-letter species code (e.g. "woothr"). The full list of valid
#'   species is in the [ebirdst_runs] data frame included in this package. To
#'   download the example dataset, use `"yebsap-example"`.
#' @param path character; directory to download the data to. All downloaded
#'   files will be placed in a sub-directory of this directory named for the
#'   data version year, e.g. "2020" for the 2020 Status Data Products. Each
#'   species' data package will then appear in a directory named with the eBird
#'   species code. Defaults to a persistent data directory, which can be found
#'   by calling `ebirdst_data_dir()`.
#' @param download_abundance whether to download estimates of abundance and
#'   proportion of population.
#' @param download_occurrence logical; whether to download estimates of
#'   occurrence.
#' @param download_count logical; whether to download estimates of count.
#' @param download_ranges logical; whether to download the range polygons.
#' @param download_regional logical; whether to download the regional summary
#'   stats, e.g. percent of population in regions.
#' @param download_pis logical; whether to download spatial estimates of
#'   predictor importance.
#' @param download_ppms logical; whether to download spatial predictive
#'   performance metrics.
#' @param download_all logical; download all files in the data package.
#'   Equivalent to setting all the `download_` arguments to `TRUE`.
#' @param pattern character; regular expression pattern to supply to
#'   [str_detect()][stringr::str_detect()] to filter files to download. This
#'   filter will be applied in addition to any of the `download_` arguments.
#'   Note that some files are mandatory and will always be downloaded.
#' @param dry_run logical; whether to do a dry run, just listing files that will
#'   be downloaded. This can be useful when testing the use of `pattern` to
#'   filter the files to download.
#' @param force logical; if the data have already been downloaded, should a
#'   fresh copy be downloaded anyway.
#' @param show_progress logical; whether to print download progress information.
#'   Defaults to `interactive()`, so downloads are silent in non-interactive
#'   sessions (e.g. scripts and R Markdown).
#'
#' @details The complete data package for each species contains a large number
#'   of files, all of which are cataloged in the vignettes. Most users will only
#'   require a small subset of these files, so by default this function only
#'   downloads the most commonly used files: GeoTIFFs providing estimate of
#'   relative abundance and proportion of population. For those interested in
#'   additional data products, the arguments starting with `download_` control
#'   the download of these other products. The `pattern` argument provides even
#'   finer grained control over what gets downloaded.
#'
#' @return Path to the folder containing the downloaded data package for the
#'   given species. If `dry_run = TRUE` a list of files to download will be
#'   returned.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # download the example data
#' ebirdst_download_status("yebsap-example")
#'
#' # download the data package for wood thrush
#' ebirdst_download_status("woothr")
#'
#' # use pattern to only download low resolution (27 km) geotiff data
#' # dry_run can be used to see what files will be downloaded
#' ebirdst_download_status("lobcur", pattern = "_27km_", dry_run = TRUE)
#' # use pattern to only download high resolution (3 km) weekly abundance data
#' ebirdst_download_status("lobcur", pattern = "abundance_median_3km",
#'                         dry_run = TRUE)
#' }
ebirdst_download_status <- function(
  species,
  path = ebirdst_data_dir(),
  download_abundance = TRUE,
  download_occurrence = FALSE,
  download_count = FALSE,
  download_ranges = FALSE,
  download_regional = FALSE,
  download_pis = FALSE,
  download_ppms = FALSE,
  download_all = FALSE,
  pattern = NULL,
  dry_run = FALSE,
  force = FALSE,
  show_progress = interactive()
) {
  stopifnot(is.character(species), length(species) == 1)
  stopifnot(is.character(path), length(path) == 1)
  stopifnot(
    is_flag(download_abundance),
    is_flag(download_occurrence),
    is_flag(download_count),
    is_flag(download_ranges),
    is_flag(download_regional),
    is_flag(download_pis),
    is_flag(download_ppms),
    is_flag(download_all)
  )
  stopifnot(is_flag(dry_run))
  stopifnot(is_flag(force))
  stopifnot(is_flag(show_progress))

  # convert to species code
  species <- get_species(species)
  if (is.na(species)) {
    stop(
      "The requested species was not modeled by Status and Trends. ",
      "Consult ebirdst_runs for a complete list of available species."
    )
  }

  # complete list of all available files for this species
  keys <- list_object_keys(species_code = species, dataset = "status")

  # decide which files to download
  keys <- select_status_keys(
    keys,
    download_abundance = download_abundance,
    download_occurrence = download_occurrence,
    download_count = download_count,
    download_ranges = download_ranges,
    download_regional = download_regional,
    download_pis = download_pis,
    download_ppms = download_ppms,
    download_all = download_all,
    pattern = pattern
  )

  # path to data package
  run_path <- file.path(
    path,
    ebirdst_version()[["status_version_year"]],
    species
  )

  # print files to download for dry run
  if (dry_run) {
    message("Downloading Status Data Products for ", species, " to:\n  ", path)
    message(paste(c("File list:", keys), collapse = "\n  "))
    return(invisible(keys))
  }

  if (show_progress) {
    message(stringr::str_glue("Downloading Status Data Products for {species}"))
  }

  fetch_data(
    keys,
    path = path,
    force = force,
    show_progress = show_progress,
    report_existing = TRUE
  )

  return(invisible(normalizePath(run_path)))
}


#' Download eBird Trends Data Products
#'
#' Download eBird Trends Data Products for set of species, or for an example
#' species. Downloading Status and Trends data requires an access key, consult
#' [set_ebirdst_access_key()] for instructions on how to obtain and store this
#' key. The example data consist of the results for Yellow-bellied Sapsucker
#' subset to Michigan and are much smaller than the full dataset, making these
#' data quicker to download and process. The example data are accessible without
#' an access key.
#'
#' @param species character; one or more species given as scientific names,
#'   common names or six-letter species codes (e.g. "woothr"). The full list of
#'   valid species can be viewed in the [ebirdst_runs] data frame included in
#'   this package; species with trends estimates are indicated by the
#'   `has_trends` column. To access the example dataset, use "yebsap-example".
#' @inheritParams ebirdst_download_status
#'
#' @return Character vector of paths to the folders containing the downloaded
#'   data packages for the given species. The trends data will be in the
#'   `trends/` subdirectory.
#' @export
#'
#' @examples
#' \dontrun{
#' # download the example data
#' ebirdst_download_trends("yebsap-example")
#'
#' # download the data package for wood thrush
#' ebirdst_download_trends("woothr")
#'
#' # multiple species can be downloaded at once
#' ebirdst_download_trends(c("Sage Thrasher", "Abert's Towhee"))
#' }
ebirdst_download_trends <- function(
  species,
  path = ebirdst_data_dir(),
  force = FALSE,
  show_progress = interactive()
) {
  stopifnot(is.character(species))
  stopifnot(is.character(path), length(path) == 1)
  stopifnot(is_flag(force))
  stopifnot(is_flag(show_progress))

  # convert to species code
  species_code <- get_species(species)
  if (anyNA(species_code)) {
    stop(
      "The following species were not modeled by Status and Trends. ",
      "Consult ebirdst_runs for a complete list of available species.\n  ",
      paste0(species[is.na(species_code)], collapse = ", ")
    )
  }

  # check that trends are available
  trends_species <- ebirdst::ebirdst_runs[ebirdst::ebirdst_runs$has_trends, ]
  no_trends <- species[!species_code %in% trends_species$species_code]
  if (length(no_trends) > 0) {
    stop(
      "The following species do not have trends estimates. ",
      "Species with trends estimates are identified in ebirdst_runs by the ",
      "has_trends column.\n  ",
      paste0(no_trends, collapse = ", ")
    )
  }

  run_paths <- character()
  for (s in species_code) {
    # complete list of all available files for this species
    keys <- list_object_keys(species_code = s, dataset = "trends")
    # only trends files
    keys <- keys[stringr::str_detect(keys, "/trends/")]

    # path to data package
    run_path <- file.path(path, ebirdst_version()[["trends_version_year"]], s)

    # download
    if (show_progress) {
      message(stringr::str_glue("Downloading Trends Data Products for {s}"))
    }
    fetch_data(
      keys,
      path = path,
      force = force,
      show_progress = show_progress,
      report_existing = TRUE
    )
    run_paths <- c(run_paths, run_path)
  }

  return(invisible(normalizePath(run_paths)))
}


#' Download eBird Status and Trends Data Coverage Products
#'
#' In addition to the species-specific data products, the eBird Status data
#' products include two products providing estimates of weekly data coverage at
#' 3 km spatial resolution: site selection probability and spatial coverage.
#' This function downloads these data products in raster GeoTIFF format.
#'
#' @inheritParams ebirdst_download_status
#'
#' @return Path to the folder containing the downloaded data coverage products.
#' @export
#'
#' @examples
#' \dontrun{
#' # download all data coverage products
#' ebirdst_download_data_coverage()
#'
#' # download just the spatial coverage products
#' ebirdst_download_data_coverage(pattern = "spatial-coverage")
#'
#' # download a single week of data coverage products
#' ebirdst_download_data_coverage(pattern = "01-04")
#'
#' # download all weeks in april
#' ebirdst_download_data_coverage(pattern = "04-")
#' }
ebirdst_download_data_coverage <- function(
  path = ebirdst_data_dir(),
  pattern = NULL,
  dry_run = FALSE,
  force = FALSE,
  show_progress = interactive()
) {
  stopifnot(is.character(path), length(path) == 1)
  stopifnot(is_flag(dry_run))
  stopifnot(is_flag(force))
  stopifnot(is_flag(show_progress))

  # complete list of all available files for this species
  keys <- list_object_keys(species_code = "data_coverage", dataset = "status")
  # path to data package
  run_path <- file.path(
    path,
    ebirdst_version()[["status_version_year"]],
    "data_coverage"
  )

  # apply pattern
  if (!is.null(pattern)) {
    stopifnot(is.character(pattern), length(pattern) == 1, !is.na(pattern))
    pat_match <- stringr::str_detect(basename(keys), pattern = pattern)
    if (!any(pat_match)) {
      stop("No files matched pattern")
    }
    keys <- keys[pat_match]
  }

  # print files to download for dry run
  if (dry_run) {
    message("Downloading Data Coverage Products to:\n  ", path)
    message(paste(c("File list:", keys), collapse = "\n  "))
    return(invisible(keys))
  }

  if (show_progress) {
    message(stringr::str_glue("Downloading Data Coverage Products"))
  }

  fetch_data(
    keys,
    path = path,
    force = force,
    show_progress = show_progress,
    report_existing = TRUE
  )

  return(invisible(normalizePath(run_path)))
}


#' Get the path to the data package for a given species
#'
#' This helper function can be used to get the path to a data package for a
#' given species.
#'
#' @param check_downloaded logical; raise an error if no data have been
#'   downloaded for this species.
#' @param dataset character; whether the path to the Status or Trends data
#'   products should be returned.
#' @inheritParams ebirdst_download_status
#'
#' @return The path to the data package directory.
#' @export
#'
#' @examples
#' \dontrun{
#' # get the path
#' path <- get_species_path("yebsap-example")
#'
#' # get the path to the full data package for yellow-bellied sapsucker
#' # common name, scientific name, or species code can be used
#' path <- get_species_path("Yellow-bellied Sapsucker")
#' path <- get_species_path("Sphyrapicus varius")
#' path <- get_species_path("yebsap")
#' }
get_species_path <- function(
  species,
  path = ebirdst_data_dir(),
  dataset = c("status", "trends"),
  check_downloaded = TRUE
) {
  stopifnot(is.character(species), length(species) == 1)
  stopifnot(is.character(path), length(path) == 1, dir.exists(path))
  stopifnot(is_flag(check_downloaded))
  dataset <- match.arg(dataset)

  if (species == "data_coverage") {
    species_code <- "data_coverage"
  } else {
    species_code <- get_species(species)
  }
  if (is.na(species_code)) {
    stop(species, " does not correspond to a valid Status and Trends species.")
  }
  version_year <- ebirdst_version()[[paste0(dataset, "_version_year")]]
  species_path <- path.expand(file.path(path, version_year, species_code))
  if (check_downloaded && !dir.exists(species_path)) {
    stop("No data package found for species: ", species)
  }
  return(species_path)
}

#' Path to eBird Status and Trends data download directory
#'
#' Identify and return the path to the default download directory for eBird
#' Status and Trends data products. This directory can be defined by setting the
#' environment variable `EBIRDST_DATA_DIR`, otherwise the directory returned by
#' `tools::R_user_dir("ebirdst", which = "data")` will be used.
#'
#' @return The path to the data download directory.
#' @export
#'
#' @examples
#' ebirdst_data_dir()
ebirdst_data_dir <- function() {
  env_var <- Sys.getenv("EBIRDST_DATA_DIR")
  if (is.null(env_var) || env_var == "") {
    return(tools::R_user_dir("ebirdst", which = "data"))
  } else {
    return(env_var)
  }
}


#' eBird Status and Trends Data Products version
#'
#' Identify the version of the eBird Status and Trends Data Products that this
#' version of the R package works with. Versions are defined by the year that
#' all model estimates are made for.
#'
#' @return A list with three components: `status_version_year` is the version year for
#'   the eBird Status Data Products, `trends_version_year` is the version year for the
#'   eBird Trends Data Products, `release_year` is the year this version of the
#'   data were released.
#' @export
#'
#' @examples
#' ebirdst_version()
ebirdst_version <- function() {
  list(
    status_version_year = 2023,
    trends_version_year = 2022,
    release_year = 2025
  )
}
