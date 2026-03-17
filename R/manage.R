#' Inventory of downloaded eBird Status and Trends data
#'
#' Returns a summary of all eBird Status and Trends data packages currently
#' downloaded to disk, with separate rows for the Status and Trends data
#' products for each species.
#'
#' @param path character; directory where data are stored. Defaults to
#'   [ebirdst_data_dir()].
#'
#' @return A tibble of class `ebirdst_inventory` with one row per data package
#'   found on disk, with columns `species_code`, `common_name`,
#'   `scientific_name`, `version_year`, `dataset` ("status" or "trends"),
#'   `n_files`, and `size_mb`. The object has a compact print method that
#'   displays the inventory grouped by version year and dataset.
#'
#' @export
#' @examples
#' \dontrun{
#' # inventory of all downloaded data
#' ebirdst_data_inventory()
#'
#' # inventory for a specific directory
#' ebirdst_data_inventory("/path/to/data")
#' }
ebirdst_data_inventory <- function(path = ebirdst_data_dir()) {
  stopifnot(is.character(path), length(path) == 1)

  empty <- dplyr::tibble(
    species_code = character(),
    common_name = character(),
    scientific_name = character(),
    version_year = integer(),
    dataset = character(),
    n_files = integer(),
    size_mb = numeric()
  )
  class(empty) <- c("ebirdst_inventory", class(empty))

  if (!dir.exists(path)) {
    return(empty)
  }

  # discover version year directories (4-digit numeric names only)
  year_dirs <- list.dirs(path, recursive = FALSE, full.names = FALSE)
  year_dirs <- year_dirs[grepl("^[0-9]{4}$", year_dirs)]

  if (length(year_dirs) == 0) {
    return(empty)
  }

  rows <- list()
  for (yr in year_dirs) {
    yr_path <- file.path(path, yr)
    sp_dirs <- list.dirs(yr_path, recursive = FALSE, full.names = TRUE)

    for (sp_dir in sp_dirs) {
      sp_code <- basename(sp_dir)
      all_files <- list.files(sp_dir, recursive = TRUE, full.names = TRUE)

      # files in the trends/ subdirectory are trends data products; all others
      # are status data products
      trends_dir <- file.path(sp_dir, "trends")
      if (dir.exists(trends_dir)) {
        trends_files <- list.files(trends_dir, recursive = TRUE,
                                   full.names = TRUE)
      } else {
        trends_files <- character(0)
      }
      status_files <- setdiff(all_files, trends_files)

      if (length(status_files) > 0) {
        rows[[length(rows) + 1]] <- dplyr::tibble(
          species_code = sp_code,
          version_year = as.integer(yr),
          dataset = "status",
          n_files = length(status_files),
          size_mb = sum(file.info(status_files)$size, na.rm = TRUE) / 1e6
        )
      }
      if (length(trends_files) > 0) {
        rows[[length(rows) + 1]] <- dplyr::tibble(
          species_code = sp_code,
          version_year = as.integer(yr),
          dataset = "trends",
          n_files = length(trends_files),
          size_mb = sum(file.info(trends_files)$size, na.rm = TRUE) / 1e6
        )
      }
    }
  }

  if (length(rows) == 0) {
    return(empty)
  }

  result <- dplyr::bind_rows(rows)

  # join with ebirdst_runs for species names
  runs_sub <- ebirdst::ebirdst_runs[, c("species_code", "common_name",
                                        "scientific_name")]
  result <- dplyr::left_join(result, runs_sub, by = "species_code")

  # data coverage products are not a species but have a known common name
  result$common_name[result$species_code == "data_coverage"] <- "Data Coverage"

  result <- result[, c("species_code", "common_name", "scientific_name",
                        "version_year", "dataset", "n_files", "size_mb")]
  result <- dplyr::arrange(result, .data$version_year, .data$species_code,
                            .data$dataset)
  class(result) <- c("ebirdst_inventory", class(result))
  return(result)
}


#' Delete downloaded eBird Status and Trends data
#'
#' Deletes downloaded eBird Status and Trends data packages for specified
#' species and/or version years. When called interactively without
#' `force = TRUE`, prints a summary of the data to be deleted and prompts for
#' confirmation before proceeding.
#'
#' @param species character; one or more species given as eBird species codes,
#'   scientific names, or English common names. If NULL (the default), data for
#'   all species are included.
#' @param year integer; one or more version years. If NULL (the default), data
#'   for all years are included.
#' @param path character; directory where data are stored. Defaults to
#'   [ebirdst_data_dir()].
#' @param force logical; if TRUE, skip the interactive confirmation prompt and
#'   delete without asking. Required when running in a non-interactive session.
#'
#' @return Invisibly returns a character vector of the paths of the deleted
#'   directories.
#'
#' @export
#' @examples
#' \dontrun{
#' # review and confirm deletion of example data
#' ebirdst_delete(species = "yebsap-example")
#'
#' # delete all data for a given version year without prompting
#' ebirdst_delete(year = 2021, force = TRUE)
#'
#' # delete a specific species and year
#' ebirdst_delete(species = "Yellow-bellied Sapsucker", year = 2022,
#'                force = TRUE)
#' }
ebirdst_delete <- function(species = NULL, year = NULL,
                            path = ebirdst_data_dir(), force = FALSE) {
  stopifnot(is.character(path), length(path) == 1)
  stopifnot(is_flag(force))
  if (!is.null(species)) stopifnot(is.character(species), length(species) >= 1)
  if (!is.null(year)) {
    stopifnot(is_integer(year), length(year) >= 1, all(year > 0))
    year <- as.integer(year)
  }

  inv <- ebirdst_data_inventory(path)
  if (nrow(inv) == 0) {
    message("No eBird Status and Trends data found in: ", path)
    return(invisible(character(0)))
  }

  # filter by year
  if (!is.null(year)) {
    inv <- inv[inv$version_year %in% year, ]
  }

  # filter by species; data_coverage passes through directly since it is not
  # in ebirdst_runs and will not be resolved by get_species()
  if (!is.null(species)) {
    is_data_cov <- tolower(trimws(species)) == "data_coverage"
    resolved_codes <- character(0)
    if (any(!is_data_cov)) {
      codes <- get_species(species[!is_data_cov])
      unrecognized <- species[!is_data_cov][is.na(codes)]
      if (length(unrecognized) > 0) {
        warning("Unrecognized species, skipping: ",
                paste(unrecognized, collapse = ", "), call. = FALSE)
      }
      resolved_codes <- codes[!is.na(codes)]
    }
    if (any(is_data_cov)) {
      resolved_codes <- c(resolved_codes, "data_coverage")
    }
    inv <- inv[inv$species_code %in% resolved_codes, ]
  }

  if (nrow(inv) == 0) {
    message("No matching data found.")
    return(invisible(character(0)))
  }

  # build unique target directories (one per species-year regardless of dataset,
  # since both status and trends data reside in the same directory)
  target_dirs <- unique(file.path(path, inv$version_year, inv$species_code))

  # safety check: all targets must be within the base path
  norm_base <- normalizePath(path, mustWork = FALSE)
  norm_targets <- normalizePath(target_dirs, mustWork = FALSE)
  safe <- startsWith(norm_targets, paste0(norm_base, .Platform$file.sep))
  if (!all(safe)) {
    stop("Safety check failed: some target directories are outside the base ",
         "path.")
  }

  # non-interactive guard
  if (!force && !interactive()) {
    stop("Cannot prompt for confirmation in a non-interactive session. ",
         "Use force = TRUE to delete without prompting.")
  }

  if (!force) {
    cat(sprintf("The following data packages will be deleted from:\n  %s\n\n",
                path))
    print(inv)
    cat("\n")

    answer <- readline("Are you sure you want to delete these data? [y/N] ")
    if (!tolower(trimws(answer)) %in% c("y", "yes")) {
      message("Deletion cancelled.")
      return(invisible(character(0)))
    }
  }

  # delete
  deleted_paths <- character(0)
  for (d in target_dirs) {
    if (unlink(d, recursive = TRUE) == 0) {
      deleted_paths <- c(deleted_paths, d)
    } else {
      warning("Failed to delete: ", d, call. = FALSE)
    }
  }

  # remove any year directories that are now empty
  affected_years <- unique(file.path(path, inv$version_year))
  for (yr_dir in affected_years) {
    if (dir.exists(yr_dir)) {
      # full.names = FALSE returns only subdir names, not yr_dir itself
      if (length(list.dirs(yr_dir, recursive = FALSE, full.names = FALSE)) == 0 &&
          length(list.files(yr_dir, recursive = TRUE)) == 0) {
        unlink(yr_dir, recursive = TRUE)
      }
    }
  }

  message("Deleted ", length(deleted_paths), " director",
          if (length(deleted_paths) == 1) "y" else "ies",
          " (", format_size(sum(inv$size_mb) * 1e6), ").")
  return(invisible(deleted_paths))
}


#' @param x an `ebirdst_inventory` object as returned by
#'   [ebirdst_data_inventory()].
#' @param ... ignored.
#'
#' @rdname ebirdst_data_inventory
#' @export
print.ebirdst_inventory <- function(x, ...) {
  n_sp <- length(unique(x$species_code))
  n_pkg <- nrow(x)
  total_size <- format_size(sum(x$size_mb) * 1e6)

  cat(sprintf("eBird Status and Trends data: %d %s, %d %s (%s)\n",
              n_sp, "species",
              n_pkg, if (n_pkg == 1L) "package" else "packages",
              total_size))

  if (n_pkg == 0L) {
    return(invisible(x))
  }

  # iterate over unique (version_year, dataset) groups in data order
  groups <- unique(x[, c("version_year", "dataset")])
  for (i in seq_len(nrow(groups))) {
    yr <- groups$version_year[i]
    ds <- groups$dataset[i]
    grp <- x[x$version_year == yr & x$dataset == ds, ]

    ds_label <- if (ds == "status") "Status" else "Trends"
    grp_size <- format_size(sum(grp$size_mb) * 1e6)
    cat(sprintf("\n%d %s Data Products (%s)\n", yr, ds_label, grp_size))

    for (j in seq_len(nrow(grp))) {
      nm <- if (!is.na(grp$common_name[j])) {
        sprintf("%s (%s)", grp$common_name[j], grp$species_code[j])
      } else {
        grp$species_code[j]
      }
      sp_size <- format_size(grp$size_mb[j] * 1e6)
      cat(sprintf("  %s: %d files, %s\n", nm, grp$n_files[j], sp_size))
    }
  }

  return(invisible(x))
}


# internal ----

format_size <- function(bytes) {
  if (bytes >= 1e9) sprintf("%.1f GB", bytes / 1e9)
  else if (bytes >= 1e6) sprintf("%.1f MB", bytes / 1e6)
  else if (bytes >= 1e3) sprintf("%.1f KB", bytes / 1e3)
  else sprintf("%.0f B", bytes)
}
