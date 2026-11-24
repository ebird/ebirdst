#' Compute weekly proportion of population within regions
#'
#' Given one or more species and a region (defined as an [sf][sf::sf],
#' [sfc][sf::sfc], or [SpatVector][terra::SpatVector] polygon object), compute
#' the proportion of each species' global modeled population (or the population
#' within `reference_region`, if supplied) that falls within each region polygon
#' for each of the 52 weeks of the year. The proportion-of-population raster
#' data are downloaded automatically on first use.
#'
#' @param species character; one or more species given as scientific names,
#'   common names, or six-letter species codes. Mixes of name types are
#'   accepted. The full list of valid species is in [ebirdst_runs].
#' @param region [sf][sf::sf], [sfc][sf::sfc], or
#'   [SpatVector][terra::SpatVector]; polygon features defining the region(s) of
#'   interest. Each feature is treated as a separate region.
#' @param reference_region `NULL`, [sf][sf::sf], [sfc][sf::sfc], or
#'   [SpatVector][terra::SpatVector]; optional polygon defining the reference
#'   population relative to which proportions are calculated. If multiple
#'   features are present, they are unioned into a single polygon. When
#'   supplied, proportions are calculated for both `region` and
#'   `reference_region`, and `proportion_population` is computed as the
#'   proportion in `region` divided by the proportion in `reference_region`
#'   (returning `0` when both are `0`, and `NA` when only the reference
#'   proportion is `0`). This can be used, for example, to calculate
#'   continental proportions, such as the proportion of the North American
#'   population within the given `region`. Defaults to `NULL` (proportions are
#'   relative to the global modeled population).
#' @param resolution character; resolution of the raster data: `"27km"`
#'   (default), `"9km"`, or `"3km"`. Note that the example dataset only
#'   provides 27 km data. If 9 km or 3 km data are requested for more than 25
#'   species that require downloading, a warning is issued; in an interactive
#'   session confirmation is required before the downloads begin.
#' @param bind_regions logical; if `TRUE`, all non-geometry attribute columns
#'   from `region` are appended to the output (ignored if `region` is an
#'   [sfc][sf::sfc] object, which has no attributes). Defaults to `FALSE`.
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
#'     population (or the population within `reference_region`, if supplied)
#'     within the region for that week. Computed as a weighted sum of
#'     proportion-of-population raster cells, where boundary cells are weighted
#'     by the fraction of the cell area covered by the polygon.
#'   - `proportion_population_reference`: proportion of the species' global
#'     modeled population within `reference_region` for that week. Only present
#'     when `reference_region` is supplied.
#'   - `.region_number`: integer (1-indexed) indicating which polygon in
#'     `region` the row corresponds to. Only present when `region` has more than
#'     one feature.
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
  reference_region = NULL,
  resolution = c("27km", "9km", "3km"),
  bind_regions = FALSE,
  path = ebirdst_data_dir(),
  force = FALSE,
  show_progress = interactive()
) {
  stopifnot(is.character(species), length(species) >= 1)
  stopifnot(inherits(region, c("sf", "sfc", "SpatVector")))
  stopifnot(n_features(region) >= 1)
  stopifnot(
    is.null(reference_region) ||
      (inherits(reference_region, c("sf", "sfc", "SpatVector")) &&
        n_features(reference_region) >= 1)
  )
  stopifnot(is_flag(bind_regions))
  stopifnot(is.character(path), length(path) == 1)
  stopifnot(is_flag(force), is_flag(show_progress))
  resolution <- match.arg(resolution)

  return(prop_pop_region_extract(
    species = species,
    region = region,
    reference_region = reference_region,
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
#' Given one or more species and a region (defined as an [sf][sf::sf],
#' [sfc][sf::sfc], or [SpatVector][terra::SpatVector] polygon object), compute
#' the proportion of each species' global modeled population (or the population
#' within `reference_region`, if supplied) that falls within each region polygon
#' for each season. The proportion-of-population raster data are downloaded
#' automatically on first use.
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
#'     population (or the population within `reference_region`, if supplied)
#'     within the region during that season. Computed as a weighted sum of
#'     proportion-of-population raster cells, where boundary cells are weighted
#'     by the fraction of the cell area covered by the polygon.
#'   - `proportion_population_reference`: proportion of the species' global
#'     modeled population within `reference_region` during that season. Only
#'     present when `reference_region` is supplied.
#'   - `.region_number`: integer (1-indexed) indicating which polygon in
#'     `region` the row corresponds to. Only present when `region` has more than
#'     one feature.
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
  reference_region = NULL,
  resolution = c("27km", "9km", "3km"),
  bind_regions = FALSE,
  path = ebirdst_data_dir(),
  force = FALSE,
  show_progress = interactive()
) {
  stopifnot(is.character(species), length(species) >= 1)
  stopifnot(inherits(region, c("sf", "sfc", "SpatVector")))
  stopifnot(n_features(region) >= 1)
  stopifnot(
    is.null(reference_region) ||
      (inherits(reference_region, c("sf", "sfc", "SpatVector")) &&
        n_features(reference_region) >= 1)
  )
  stopifnot(is_flag(bind_regions))
  stopifnot(is.character(path), length(path) == 1)
  stopifnot(is_flag(force), is_flag(show_progress))
  resolution <- match.arg(resolution)

  return(prop_pop_region_extract(
    species = species,
    region = region,
    reference_region = reference_region,
    period = "seasonal",
    resolution = resolution,
    bind_regions = bind_regions,
    path = path,
    force = force,
    show_progress = show_progress
  ))
}


# internal ----

n_features <- function(x) {
  if (inherits(x, "sfc")) {
    return(length(x))
  }
  return(nrow(x))
}

prop_pop_region_extract <- function(
  species,
  region,
  reference_region,
  period,
  resolution,
  bind_regions,
  path,
  force,
  show_progress
) {
  # resolve species to codes; error on any unrecognised name
  species_codes <- get_species(species, strict = TRUE)

  # convert region to SpatVector; keep original for bind_regions join
  if (inherits(region, c("sf", "sfc"))) {
    region_vect <- terra::vect(region)
  } else {
    region_vect <- region
  }
  n_regions <- nrow(region_vect)

  # convert reference_region to a single unioned SpatVector if supplied
  if (!is.null(reference_region)) {
    if (inherits(reference_region, c("sf", "sfc"))) {
      ref_vect <- terra::vect(reference_region)
    } else {
      ref_vect <- reference_region
    }
    if (nrow(ref_vect) > 1) {
      ref_vect <- terra::aggregate(ref_vect)
    }
  } else {
    ref_vect <- NULL
  }

  prop_pop_download_check(species_codes, period, resolution, path)

  # loop over species; collect per-species data frames
  results <- vector("list", length(species_codes))
  raster_crs <- NULL
  region_proj <- NULL
  ref_proj <- NULL

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

    cur_crs <- terra::crs(raster)
    if (is.null(raster_crs) || !identical(cur_crs, raster_crs)) {
      raster_crs <- cur_crs
      region_proj <- terra::project(region_vect, raster_crs)
      if (!is.null(ref_vect)) {
        ref_proj <- terra::project(ref_vect, raster_crs)
      }
    }

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

    if (!is.null(ref_proj)) {
      ref_ext <- terra::extract(
        raster,
        ref_proj,
        fun = "sum",
        na.rm = TRUE,
        weights = TRUE,
        ID = FALSE
      )
      ref_denom <- rep(as.numeric(ref_ext[1, ]), times = n_regions)
      zero_ref <- !is.na(ref_denom) & ref_denom == 0
      zero_vals <- !is.na(vals) & vals == 0
      vals <- vals / ref_denom
      vals[zero_ref & zero_vals] <- 0
      vals[zero_ref & !zero_vals] <- NA_real_
    }

    if (period == "weekly") {
      seasons <- assign_weeks_to_seasons(
        sc,
        path = path,
        force = FALSE,
        show_progress = FALSE
      )
      df <- data.frame(
        species_code = sc,
        season = rep(seasons, times = n_regions),
        week = rep(as.Date(names(raster)), times = n_regions),
        proportion_population = vals
      )
    } else {
      df <- data.frame(
        species_code = sc,
        season = rep(names(raster), times = n_regions),
        proportion_population = vals
      )
    }
    if (!is.null(ref_proj)) {
      df$proportion_population_reference <- ref_denom
    }
    df$.region_number <- rep(seq_len(n_regions), each = n_layers)
    results[[i]] <- df
  }

  combined <- do.call(rbind, results)

  # attach region attribute columns if requested (sfc has no attributes)
  if (bind_regions && !inherits(region, "sfc")) {
    if (inherits(region, "sf")) {
      region_df <- sf::st_drop_geometry(region)
    } else {
      region_df <- as.data.frame(region)
    }
    if (ncol(region_df) > 0) {
      region_df$.region_number <- seq_len(nrow(region_df))
      combined <- dplyr::left_join(combined, region_df, by = ".region_number")
    }
  }

  # drop .region_number for single-polygon inputs (not useful as an identifier)
  if (n_regions == 1) {
    combined$.region_number <- NULL
  }

  return(combined)
}
