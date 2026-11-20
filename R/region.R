#' Compute weekly proportion of population within regions
#'
#' Given one or more species and a region (defined as an [sf][sf::sf] or
#' [SpatVector][terra::SpatVector] polygon object), compute the proportion of
#' each species' global modeled population that falls within each region
#' polygon for each of the 52 weeks of the year. The proportion-of-population
#' raster data are downloaded automatically on first use.
#'
#' @param species character; one or more species given as scientific names,
#'   common names, or six-letter species codes. Mixes of name types are
#'   accepted. The full list of valid species is in [ebirdst_runs].
#' @param region [sf][sf::sf] or [SpatVector][terra::SpatVector]; polygon
#'   features defining the region(s) of interest. Each row (feature) is
#'   treated as a separate region.
#' @param resolution character; resolution of the raster data: `"27km"`
#'   (default), `"9km"`, or `"3km"`. Note that the example dataset only
#'   provides 27 km data. If 9 km or 3 km data are requested for more than 25
#'   species that require downloading, a warning is issued; in an interactive
#'   session confirmation is required before the downloads begin.
#' @param bind_regions logical; if `TRUE`, all non-geometry attribute columns
#'   from `region` are appended to the output. Defaults to `FALSE`.
#' @inheritParams load_raster
#'
#' @return A data frame with one row per combination of species, week, and
#'   region polygon. Columns:
#'   - `species_code`: the alphanumeric eBird species code.
#'   - `season`: season the week falls within (`"breeding"`,
#'     `"nonbreeding"`, `"prebreeding_migration"`, `"postbreeding_migration"`,
#'     or `"resident"`), or `NA` if the week does not fall within a named
#'     season that passed expert review.
#'   - `week`: midpoint date (`YYYY-MM-DD`) of the week.
#'   - `proportion_population`: proportion of the species' global modeled
#'     population within the region for that week. Computed as a weighted sum
#'     of proportion-of-population raster cells, where boundary cells are
#'     weighted by the fraction of the cell area covered by the polygon.
#'   - `.region_number`: integer (1-indexed) indicating which polygon in
#'     `region` the row corresponds to. Only present when `nrow(region) > 1`.
#'   - Additional columns from `region` if `bind_regions = TRUE`.
#'
#' @export
#' @examples
#' \dontrun{
#' # download example data
#' ebirdst_download_status("yebsap-example")
#'
#' # build a region polygon from the raster extent
#' rast <- load_raster("yebsap-example", "abundance", resolution = "27km")
#' region <- sf::st_as_sf(
#'   terra::as.polygons(terra::ext(rast), crs = terra::crs(rast))
#' )
#'
#' chron <- migration_chronology("yebsap-example", region, resolution = "27km")
#' head(chron)
#' }
migration_chronology <- function(
  species,
  region,
  resolution = c("27km", "9km", "3km"),
  bind_regions = FALSE,
  path = ebirdst_data_dir(),
  force = FALSE,
  show_progress = interactive()
) {
  stopifnot(is.character(species), length(species) >= 1)
  stopifnot(inherits(region, c("sf", "SpatVector")))
  stopifnot(nrow(region) >= 1)
  stopifnot(is_flag(bind_regions))
  stopifnot(is.character(path), length(path) == 1)
  stopifnot(is_flag(force), is_flag(show_progress))
  resolution <- match.arg(resolution)

  return(prop_pop_region_extract(
    species = species,
    region = region,
    period = "weekly",
    resolution = resolution,
    bind_regions = bind_regions,
    path = path,
    force = force,
    show_progress = show_progress
  ))
}


#' Compute seasonal proportion of population within regions
#'
#' Given one or more species and a region (defined as an [sf][sf::sf] or
#' [SpatVector][terra::SpatVector] polygon object), compute the proportion of
#' each species' global modeled population that falls within each region
#' polygon for each season. The proportion-of-population raster data are
#' downloaded automatically on first use.
#'
#' @inheritParams migration_chronology
#'
#' @return A data frame with one row per combination of species, season, and
#'   region polygon. Columns:
#'   - `species_code`: the alphanumeric eBird species code.
#'   - `season`: the season name (`"breeding"`, `"nonbreeding"`,
#'     `"prebreeding_migration"`, `"postbreeding_migration"`, or
#'     `"resident"`). Only seasons that passed expert review are included.
#'   - `proportion_population`: proportion of the species' global modeled
#'     population within the region during that season. Computed as a weighted
#'     sum of proportion-of-population raster cells, where boundary cells are
#'     weighted by the fraction of the cell area covered by the polygon.
#'   - `.region_number`: integer (1-indexed) indicating which polygon in
#'     `region` the row corresponds to. Only present when `nrow(region) > 1`.
#'   - Additional columns from `region` if `bind_regions = TRUE`.
#'
#' @export
#' @examples
#' \dontrun{
#' # download example data
#' ebirdst_download_status("yebsap-example")
#'
#' # build a region polygon from the raster extent
#' rast <- load_raster("yebsap-example", "abundance", resolution = "27km")
#' region <- sf::st_as_sf(
#'   terra::as.polygons(terra::ext(rast), crs = terra::crs(rast))
#' )
#'
#' seas <- seasonal_proportion_population(
#'   "yebsap-example", region,
#'   resolution = "27km"
#' )
#' seas
#' }
seasonal_proportion_population <- function(
  species,
  region,
  resolution = c("27km", "9km", "3km"),
  bind_regions = FALSE,
  path = ebirdst_data_dir(),
  force = FALSE,
  show_progress = interactive()
) {
  stopifnot(is.character(species), length(species) >= 1)
  stopifnot(inherits(region, c("sf", "SpatVector")))
  stopifnot(nrow(region) >= 1)
  stopifnot(is_flag(bind_regions))
  stopifnot(is.character(path), length(path) == 1)
  stopifnot(is_flag(force), is_flag(show_progress))
  resolution <- match.arg(resolution)

  return(prop_pop_region_extract(
    species = species,
    region = region,
    period = "seasonal",
    resolution = resolution,
    bind_regions = bind_regions,
    path = path,
    force = force,
    show_progress = show_progress
  ))
}


# internal ----

prop_pop_region_extract <- function(
  species,
  region,
  period,
  resolution,
  bind_regions,
  path,
  force,
  show_progress
) {
  # resolve species to codes; error on any unrecognised name
  species_codes <- get_species(species)
  bad <- species[is.na(species_codes)]
  if (length(bad) > 0) {
    stop(
      "The following species were not found in ebirdst_runs. ",
      "Use get_species() to look up valid species codes.\n  ",
      paste(bad, collapse = "\n  ")
    )
  }

  # convert region to SpatVector; keep original for bind_regions join
  if (inherits(region, "sf")) {
    region_vect <- terra::vect(region)
  } else {
    region_vect <- region
  }
  n_regions <- nrow(region_vect)

  # check how many species will need downloading at the requested resolution
  v <- ebirdst_version()[["status_version_year"]]
  if (period == "weekly") {
    file_suffix <- sprintf(
      "proportion-population_median_%s_%s.tif",
      resolution,
      v
    )
    sub_dir <- "weekly"
  } else {
    file_suffix <- sprintf(
      "proportion-population_seasonal_mean_%s_%s.tif",
      resolution,
      v
    )
    sub_dir <- "seasonal"
  }
  n_to_download <- sum(vapply(
    species_codes,
    function(sc) {
      fpath <- file.path(
        path,
        v,
        sc,
        sub_dir,
        sprintf("%s_%s", sc, file_suffix)
      )
      return(!file.exists(fpath))
    },
    FUN.VALUE = logical(1)
  ))

  if (resolution != "27km" && n_to_download > 25) {
    msg <- sprintf(
      paste0(
        "Downloading %d species at %s resolution requires significant ",
        "disk space and time."
      ),
      n_to_download,
      resolution
    )
    if (interactive()) {
      message(msg)
      answer <- readline("Continue? [y/n]: ")
      if (!grepl("^[Yy]", answer)) {
        stop("Download cancelled.")
      }
    } else {
      warning(msg)
    }
  }

  # loop over species; collect per-species data frames
  results <- vector("list", length(species_codes))

  for (i in seq_along(species_codes)) {
    sc <- species_codes[i]

    if (show_progress) {
      message(sprintf("[%d/%d] %s", i, length(species_codes), sc))
    }

    raster <- load_raster(
      sc,
      product = "proportion-population",
      period = period,
      resolution = resolution,
      path = path,
      force = force,
      show_progress = show_progress
    )

    region_proj <- terra::project(region_vect, terra::crs(raster))

    ext <- terra::extract(
      raster,
      region_proj,
      fun = "sum",
      na.rm = TRUE,
      weights = TRUE
    )

    n_layers <- terra::nlyr(raster)
    # ext is n_regions x (1 + n_layers); first column is polygon ID
    # t() + as.numeric() lays out values as: poly1_layer1, poly1_layer2, ...,
    # poly1_layerN, poly2_layer1, ..., matching rep(..., each = n_layers)
    vals <- as.numeric(t(as.matrix(ext[, -1, drop = FALSE])))

    if (period == "weekly") {
      seasons <- assign_weeks_to_seasons(
        sc,
        path = path,
        force = FALSE,
        show_progress = FALSE
      )
      results[[i]] <- data.frame(
        species_code = sc,
        season = rep(seasons, times = n_regions),
        week = rep(as.Date(names(raster)), times = n_regions),
        proportion_population = vals,
        .region_number = rep(seq_len(n_regions), each = n_layers)
      )
    } else {
      results[[i]] <- data.frame(
        species_code = sc,
        season = rep(names(raster), times = n_regions),
        proportion_population = vals,
        .region_number = rep(seq_len(n_regions), each = n_layers)
      )
    }
  }

  combined <- do.call(rbind, results)

  # attach region attribute columns if requested
  if (bind_regions) {
    if (inherits(region, "sf")) {
      region_df <- sf::st_drop_geometry(region)
    } else {
      region_df <- as.data.frame(region)
    }
    region_df$.region_number <- seq_len(nrow(region_df))
    combined <- dplyr::left_join(combined, region_df, by = ".region_number")
  }

  # drop .region_number for single-polygon inputs (not useful as an identifier)
  if (n_regions == 1) {
    combined$.region_number <- NULL
  }

  return(combined)
}
