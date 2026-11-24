#' Shared stewardship connections for migratory birds
#'
#' Generate a shared stewardship connections raster mapping where migratory
#' species that occupy an origin region in one stationary season
#' (breeding/non-breeding) concentrate in the other stationary season
#' (non-breeding/breeding). For each species, the destination-season
#' proportion-of-population raster is multiplied by a species weight (the
#' fraction of the global population in the origin region), and the weighted
#' rasters are summed across all retained species.
#'
#' **Input modes:** The species and weights can be defined in three quite
#' distinct ways depending on whether `region_code`, `species_weights`, or
#' `origin_region` is supplied:
#'
#' - **Region-code mode** (`region_code`): provide a region code
#' (e.g. `"USA"`, `"USA-NY"`) matching one present in
#' [ebirdst_regional_stats()]. The candidate species list is derived
#' automatically from that region; `species` optionally subsets it. Weights and
#' thresholds are derived from `total_pop_percent` in the regional stats
#' table — no proportion-of-population rasters are downloaded for filtered-out
#' species.
#' - **Custom polygon mode** (`origin_region`): provide an `sf` or `SpatVector`
#' polygon defining the focal origin region. In this mode, a set of species must
#' be supplied via the `species` argument. For each species the weight and
#' inclusion thresholds are computed by calculating the seasonal
#' proportion-of-population within the supplied polygon.
#' - **Weights-table mode** (`species_weights`): provide a data frame with
#' columns `species_code` and `weight` (numeric, non-negative). Weights are used
#' as supplied; threshold arguments and `include_migration_seasons` are ignored.
#' `species` optionally subsets the rows.
#'
#' **Species selection note:** There is significant download and compute cost
#' for processing large numbers of species. For the custom polygon mode, the
#' typical workflow is to use [ebirdst_regional_stats()], or other datasets, to
#' coarsely subset the species list before calling this function. Region-code
#' mode does this subsetting automatically.
#'
#' This method was developed by the Cornell Lab of Ornithology and Partners in
#' Flight (PIF). The approach is described in:
#'
#' - Partners in Flight (2023). Making Connections for Bird Conservation, 2nd
#'   Edition. PIF Technical Series No. 10.
#'   <https://partnersinflight.org/resources/pif-technical-series-10-making-connections-appendices-and-maps/>
#' - Lello-Smith et al. (2026). Biological Conservation 313: 111551.
#'   <https://doi.org/10.1016/j.biocon.2025.111551>
#'
#' **Important:** The resulting raster shows multispecies abundance
#' concentrations weighted by origin-region population representation — it is
#' **not** a migratory connectivity metric. The individuals present in the
#' destination season are not necessarily the same individuals that occupy the
#' focal region in the origin season.
#'
#' @param species character or `NULL`; one or more species given as scientific
#'   names, common names, or six-letter species codes. Mixes of name types are
#'   accepted. The full list of valid species is in [ebirdst_runs]. In custom
#'   polygon mode, `species` is required. In region-code mode, `species`
#'   optionally subsets the candidate list derived from the regional stats
#'   table. In weights-table mode, `species` optionally subsets the rows of
#'   `species_weights`.
#' @param origin_region `NULL`, [sf][sf::sf], or
#'   [SpatVector][terra::SpatVector]; polygon defining the focal origin region.
#'   Required for polygon mode; must be `NULL` when `region_code` or
#'   `species_weights` is supplied. If multiple features are present, they are
#'   unioned to a single polygon.
#' @param region_code character or `NULL`; a region code matching one
#'   present in [ebirdst_regional_stats()] (e.g. `"USA"`, `"USA-NY"`,
#'   `"TWN"`). Required for region-code mode; must be `NULL` otherwise. Use
#'   [ebirdst_regional_stats()] to browse valid codes via the `region_code`
#'   column.
#' @param species_weights data frame or `NULL`; a data frame with at least two
#'   columns: `species_code` (character, resolved via [get_species()]) and
#'   `weight` (numeric, finite, non-negative). Required for weights-table mode;
#'   must be `NULL` otherwise. Rows with duplicate resolved `species_code`
#'   values are not permitted. Threshold arguments
#'   (`min_origin_threshold`, `min_migration_threshold`) and
#'   `include_migration_seasons` are ignored in this mode.
#' @param origin_season character; the stationary season in which the focal
#'   region is occupied: `"breeding"` (default) or `"nonbreeding"`. The
#'   destination season is the opposite stationary season.
#' @param include_migration_seasons logical; if `TRUE`, the weight for each
#'   species also considers the fraction of the population within
#'   `origin_region` during the pre- and post-breeding migration seasons (taking
#'   the maximum across all origin-season layers present). The default
#'   `min_migration_threshold` is relaxed to `0.5` in this case. Applies only
#'   to polygon mode and region-code mode; ignored in weights-table mode.
#'   Defaults to `FALSE`.
#' @param min_origin_threshold numeric `[0, 1)`; minimum weight (fraction of
#'   global population in the origin region during origin season(s)) required to
#'   retain a species. Species with `weight < min_origin_threshold` are excluded
#'   before the destination-side computation. Applies to polygon mode and
#'   region-code mode; ignored in weights-table mode. Defaults to `0.0001`
#'   (1%). Set to `0` to disable.
#' @param min_migration_threshold numeric `[0, 1)`; minimum fraction by which a
#'   species' population must leave the origin region between the origin and
#'   destination stationary seasons. Applies to polygon mode and region-code
#'   mode; ignored in weights-table mode. Defaults to `0.85` when
#'   `include_migration_seasons = FALSE` and `0.5` when `TRUE`. Set to `0` to
#'   disable.
#' @param destination_mask `NULL`, [sf][sf::sf], [sf][sf::sfc],
#'   [SpatVector][terra::SpatVector], [SpatExtent][terra::SpatExtent], or a
#'   bounding box from [sf::st_bbox()]; spatial object used to crop and mask the
#'   destination-season raster for each species before accumulation. If `NULL`,
#'   the full global raster is used — a `destination_mask` is strongly
#'   recommended for broadly distributed species to reduce memory and compute.
#' @param resolution character; resolution of the raster data: `"27km"`
#'   (default), `"9km"`, or `"3km"`. Note that the example dataset only provides
#'   27 km data. If more than 25 species are requested at finer-than-27km
#'   resolution, or more than 100 at 27km, a message is issued about expected
#'   download size and compute time.
#' @param filename character or `NULL`; file path to write the connection raster
#'   to disk via [terra::writeRaster()]. If `NULL` (default), the raster is
#'   returned in memory only. When `filename` is supplied the return value is
#'   invisible and additional arguments in `...` are passed to
#'   [terra::writeRaster()].
#' @inheritParams load_raster
#' @param ... additional arguments passed to [terra::writeRaster()] when
#'   `filename` is supplied.
#'
#' @return A named list with two elements:
#'   - `connection`: a single-layer [SpatRaster][terra::SpatRaster] giving the
#'     stewardship connection surface — the weighted sum of destination-season
#'     proportion-of-population rasters, with each species weighted by the
#'     fraction of its global population in the origin region.
#'   - `weights`: a data frame with one row per retained species and columns:
#'     - `species_code`: the alphanumeric eBird species code.
#'     - `common_name`: the English common name.
#'     - `origin_season`: which origin-season produced the proportion of
#'     population `NA` in weights-table mode.
#'     - `weight`: fraction of global population in the origin region during the
#'     origin season(s). In region-code mode this is `total_pop_percent` from
#'     [ebirdst_regional_stats()]; in weights-table mode it is the user-supplied
#'     value.
#'     - `rel_pct_decrease`: fraction of the origin-season proportion of
#'     population that left the region by the destination season. `NA` in
#'     weights-table mode.
#' @export
#' @examples
#' \dontrun{
#' # Appalachian warblers with breeding populations in West Virginia
#' warblers <- c("cerwar", "woewar1", "louwat", "buwwar", "hoowar")
#'
#' # West Virginia: rough bounding box in WGS84
#' origin_region <- sf::st_sf(geometry = sf::st_as_sfc(
#'   sf::st_bbox(
#'     c(xmin = -82.6, xmax = -77.7, ymin = 37.2, ymax = 40.6),
#'     crs = sf::st_crs(4326)
#'   )
#' ))
#'
#' # Central America destination mask (Central Mexico through Panama)
#' destination_mask <- sf::st_bbox(
#'     c(xmin = -99, xmax = -80, ymin = 10, ymax = 22),
#'     crs = sf::st_crs(4326)
#' )
#'
#' # region-code mode: weights from ebirdst_regional_stats()
#' result_rc <- stewardship_connections(
#'   warblers,
#'   region_code = "USA-WV",
#'   origin_season = "breeding",
#'   destination_mask = destination_mask,
#'   resolution = "27km"
#' )
#' result_rc$weights
#'
#' # polygon mode: weights derived from proportion-population rasters
#' result_poly <- stewardship_connections(
#'   warblers,
#'   origin_region = origin_region,
#'   origin_season = "breeding",
#'   destination_mask = destination_mask,
#'   resolution = "27km"
#' )
#' result_poly$weights
#'
#' # weights-table mode: supply weights directly, skips threshold filtering
#' wt <- data.frame(
#'   species_code = warblers,
#'   weight = c(0.12, 0.08, 0.15, 0.20, 0.10)
#' )
#' result_wt <- stewardship_connections(
#'   species_weights = wt,
#'   origin_season = "breeding",
#'   destination_mask = destination_mask,
#'   resolution = "27km"
#' )
#'
#' # display the connection surface
#' brks <- importance_breaks(result_poly$connection)
#' pal <- c("#e6e6e6", ebirdst_palettes(length(brks) - 2))
#' terra::plot(
#'   result_poly$connection,
#'   breaks = brks,
#'   col = pal,
#'   axes = FALSE,
#'   box = TRUE,
#'   legend = FALSE,
#'   mar = c(0.25, 0.25, 0.25, 0.25)
#' )
#' }
stewardship_connections <- function(
  species = NULL,
  origin_region = NULL,
  region_code = NULL,
  species_weights = NULL,
  origin_season = c("breeding", "nonbreeding"),
  include_migration_seasons = FALSE,
  min_origin_threshold = 0.0001,
  min_migration_threshold = if (include_migration_seasons) 0.5 else 0.85,
  destination_mask = NULL,
  resolution = c("27km", "9km", "3km"),
  filename = NULL,
  path = ebirdst_data_dir(),
  force = FALSE,
  show_progress = interactive(),
  ...
) {
  stopifnot(is_flag(include_migration_seasons))
  stopifnot(
    is.numeric(min_origin_threshold),
    length(min_origin_threshold) == 1,
    min_origin_threshold >= 0,
    min_origin_threshold <= 1
  )
  stopifnot(
    is.numeric(min_migration_threshold),
    length(min_migration_threshold) == 1,
    min_migration_threshold >= 0,
    min_migration_threshold <= 1
  )
  stopifnot(is.character(path), length(path) == 1)
  stopifnot(is_flag(force), is_flag(show_progress))
  stopifnot(
    is.null(filename) || (is.character(filename) && length(filename) == 1)
  )
  origin_season <- match.arg(origin_season)
  resolution <- match.arg(resolution)

  # exactly one input mode must be supplied
  n_modes <- (!is.null(origin_region)) +
    (!is.null(region_code)) +
    (!is.null(species_weights))
  if (n_modes != 1L) {
    stop(
      "Exactly one of `origin_region`, `region_code`, or `species_weights` ",
      "must be supplied."
    )
  }

  if (origin_season == "breeding") {
    destination_season <- "nonbreeding"
  } else {
    destination_season <- "breeding"
  }
  origin_seasons <- origin_season
  if (include_migration_seasons) {
    origin_seasons <- c(
      origin_seasons,
      "prebreeding_migration",
      "postbreeding_migration"
    )
  }

  # phase 1: derive retained species and weights (mode-specific)
  if (!is.null(origin_region)) {
    if (is.null(species)) {
      stop(
        "`species` is required when `origin_region` is supplied."
      )
    }
    weights <- stewardship_weights_polygon(
      species = species,
      origin_region = origin_region,
      origin_seasons = origin_seasons,
      destination_season = destination_season,
      min_origin_threshold = min_origin_threshold,
      min_migration_threshold = min_migration_threshold,
      resolution = resolution,
      path = path,
      force = force,
      show_progress = show_progress
    )
  } else if (!is.null(region_code)) {
    weights <- stewardship_weights_region_code(
      region_code = region_code,
      species = species,
      origin_seasons = origin_seasons,
      destination_season = destination_season,
      min_origin_threshold = min_origin_threshold,
      min_migration_threshold = min_migration_threshold,
      path = path,
      force = force,
      show_progress = show_progress
    )
  } else {
    weights <- stewardship_weights_table(
      species_weights = species_weights,
      species = species
    )
  }

  weights_df <- weights$weights_df
  n_dropped_origin <- weights$n_dropped_origin
  n_dropped_migration <- weights$n_dropped_migration

  # selection criteria message
  msg <- stringr::str_glue(
    "{nrow(weights_df)} species retained after applying selection thresholds.",
    "{n_dropped_origin} dropped due to low origin proportion of population.",
    "{n_dropped_migration} dropped due to low migration out of origin region.",
    .sep = "\n  ",
    .trim = FALSE
  )
  if (nrow(weights_df) == 0L) {
    stop(
      msg,
      "\nConsider lowering min_origin_threshold or min_migration_threshold."
    )
  }
  message(msg)

  # phase 2: accumulate destination rasters weighted by species weights
  retained_codes <- weights_df$species_code
  prop_pop_download_check(retained_codes, "seasonal", resolution, path)

  running_total <- NULL
  raster_crs <- NULL
  mask_info <- NULL

  if (show_progress) {
    message("Generating connections raster for ", nrow(weights_df), " species")
  }
  for (i in seq_len(nrow(weights_df))) {
    sc <- weights_df$species_code[i]
    weight <- weights_df$weight[i]

    if (show_progress) {
      message(sprintf("  [%d/%d] %s", i, nrow(weights_df), sc))
    }

    r <- load_raster(
      sc,
      product = "proportion-population",
      period = "seasonal",
      resolution = resolution,
      path = path,
      force = force,
      show_progress = FALSE
    )

    if (is.null(raster_crs)) {
      raster_crs <- terra::crs(r)
      mask_info <- normalize_destination_mask(destination_mask, raster_crs)
    }

    if (!(destination_season %in% names(r))) {
      stop(
        destination_season,
        " season layer missing from seasonal raster for",
        sc
      )
    }

    dest_layer <- r[[destination_season]]
    if (!is.null(mask_info$crop_ext)) {
      dest_layer <- terra::crop(dest_layer, mask_info$crop_ext, snap = "out")
    }
    if (!is.null(mask_info$mask_vect)) {
      dest_layer <- terra::mask(dest_layer, mask_info$mask_vect)
    }

    # accumulate one species at a time rather than stacking then calling
    # app(sum): at 3km each layer is ~1 GB, so a 100-species stack would
    # require ~100 GB RAM. sequential addition holds only two rasters in
    # memory at once and benchmarks 3-4x faster than stack + app() due to
    # terra loading all layers simultaneously during chunk processing.
    # similarly, periodic writeRaster/rast() roundtrips to cap gc pressure
    # were benchmarked and found to give no meaningful memory or speed benefit.
    contrib <- dest_layer * weight
    if (is.null(running_total)) {
      running_total <- contrib
    } else {
      # sum required here instead of + since we need to drop missing values
      running_total <- sum(running_total, contrib, na.rm = TRUE)
    }
  }
  result <- list(connection = running_total, weights = weights_df)

  if (!is.null(filename)) {
    terra::writeRaster(running_total, filename, ...)
    return(invisible(result))
  }
  return(result)
}


#' Compute importance breaks for a stewardship connections raster
#'
#' Derives quantile-based bin breaks for displaying a stewardship connections
#' raster returned by [stewardship_connections()]. The lower tail of the
#' distribution (cells at or below the `lower_quantile` of non-zero values) is
#' collapsed to a single low-importance class. The upper portion is divided
#' into `n_bins` equal-quantile bins.
#'
#' Colors for the bins can be obtained from [ebirdst_palettes()] with
#' `n = n_bins`. The function deliberately does not return a palette so that
#' colors remain under caller control.
#'
#' @param x [SpatRaster][terra::SpatRaster]; a single-layer stewardship
#'   connection raster, typically the `$connection` element returned by
#'   [stewardship_connections()].
#' @param lower_quantile numeric; quantile (in `[0, 1)`) of the non-zero cell
#'   values used as the lower cutoff. Cells below this quantile are pooled into
#'   a low-importance class (returned as `NA` when `return_raster = TRUE`).
#'   Defaults to `0.5` (the median of non-zero cells).
#' @param n_bins integer; number of display bins for cells at or above the
#'   `lower_quantile` cutoff. Defaults to `10`.
#' @param return_raster logical; if `FALSE` (default), return the numeric break
#'   values. If `TRUE`, return a reclassified
#'   [SpatRaster][terra::SpatRaster] with integer cell values in `1:n_bins`;
#'   cells below the cutoff are `NA`.
#'
#' @return When `return_raster = FALSE` (default), a numeric vector of `n_bins
#'   + 2` break values: `0`, the `lower_quantile` cutoff, and `n_bins` upper
#'   quantile boundaries. When `return_raster = TRUE`, a single-layer
#'   [SpatRaster][terra::SpatRaster] with integer cell values in `1:n_bins`
#'   (cells below the cutoff are `NA`), ready to plot with
#'   `col = ebirdst_palettes(n_bins)`.
#' @export
#' @examples
#' # create a synthetic connection raster
#' r <- terra::rast(nrows = 10, ncols = 10, vals = runif(100, 0, 0.5))
#'
#' # default: return breaks (length = n_bins + 2 = 12)
#' brks <- importance_breaks(r)
#' length(brks)
#'
#' # return a classified raster for direct plotting
#' r_class <- importance_breaks(r, return_raster = TRUE)
#' terra::plot(r_class, col = ebirdst_palettes(10))
importance_breaks <- function(
  x,
  lower_quantile = 0.5,
  n_bins = 10,
  return_raster = FALSE
) {
  stopifnot(inherits(x, "SpatRaster"), terra::nlyr(x) == 1)
  stopifnot(
    is.numeric(lower_quantile),
    length(lower_quantile) == 1,
    lower_quantile >= 0,
    lower_quantile < 1
  )
  stopifnot(is_count(n_bins), n_bins >= 1)
  stopifnot(is_flag(return_raster))

  v <- terra::values(x, na.rm = TRUE, mat = FALSE)
  v <- v[v > 0]

  if (length(v) == 0) {
    stop("No positive values found in raster x.")
  }

  cutoff <- stats::quantile(v, lower_quantile)
  v_upper <- v[v >= cutoff]
  breaks <- c(
    0,
    unname(stats::quantile(v_upper, seq(0, 1, length.out = n_bins + 1)))
  )

  if (!return_raster) {
    return(unname(breaks))
  }

  rcl <- cbind(
    breaks[seq_len(n_bins) + 1],
    breaks[seq_len(n_bins) + 2],
    seq_len(n_bins)
  )
  r_class <- terra::classify(x, rcl, include.lowest = TRUE, right = FALSE)
  # terra::classify() retains original values for unmatched cells;
  # explicitly mask cells below the first interval's lower bound
  r_class[x < breaks[2]] <- NA
  return(r_class)
}


# internal ----

check_stewardship_species <- function(species_codes, explicit) {
  runs_sub <- ebirdst::ebirdst_runs[
    ebirdst::ebirdst_runs$species_code %in% species_codes,
  ]
  invalid <- character(0)
  invalid_reason <- character(0)
  for (sc in species_codes) {
    run <- runs_sub[runs_sub$species_code == sc, ]
    b_ok <- isTRUE(
      !is.na(run$breeding_quality) && as.integer(run$breeding_quality) > 0
    )
    nb_ok <- isTRUE(
      !is.na(run$nonbreeding_quality) && as.integer(run$nonbreeding_quality) > 0
    )
    if (isTRUE(run$is_resident)) {
      invalid <- c(invalid, sc)
      invalid_reason <- c(invalid_reason, "resident species")
    } else if (isTRUE(run$ocean_modeled)) {
      invalid <- c(invalid, sc)
      invalid_reason <- c(invalid_reason, "ocean-modeled species")
    } else if (!b_ok || !nb_ok) {
      invalid <- c(invalid, sc)
      invalid_reason <- c(
        invalid_reason,
        "failed expert review in a stationary season"
      )
    }
  }
  if (length(invalid) > 0 && explicit) {
    stop(
      "The following species cannot be used in a stewardship connections ",
      "analysis:\n  ",
      paste(sprintf("%s (%s)", invalid, invalid_reason), collapse = "\n  "),
      "\nRemove these species from the input list and try again."
    )
  }
  return(invisible(invalid))
}


stewardship_weights_region_code <- function(
  region_code,
  species,
  origin_seasons,
  destination_season,
  min_origin_threshold,
  min_migration_threshold,
  path,
  force,
  show_progress
) {
  stopifnot(is.character(region_code), length(region_code) == 1)

  stats <- ebirdst_regional_stats(
    path = path,
    force = force,
    show_progress = show_progress
  )

  region_stats <- stats[stats$region_code == region_code, ]
  if (nrow(region_stats) == 0L) {
    stop(
      "No data found for region_code '",
      region_code,
      "' in the regional stats file. ",
      "Use ebirdst_regional_stats() to browse valid region codes."
    )
  }

  candidate_codes <- unique(region_stats$species_code)

  if (!is.null(species)) {
    stopifnot(is.character(species), length(species) >= 1)
    resolved <- get_species(species, strict = TRUE)
    candidate_codes <- intersect(candidate_codes, resolved)
    if (length(candidate_codes) == 0L) {
      stop(
        "None of the supplied species are present in region_code '",
        region_code,
        "'."
      )
    }
  }

  drop_codes <- check_stewardship_species(candidate_codes, explicit = FALSE)
  candidate_codes <- setdiff(candidate_codes, drop_codes)
  n_dropped_ineligible <- length(drop_codes)

  runs_sub <- ebirdst::ebirdst_runs[
    ebirdst::ebirdst_runs$species_code %in% candidate_codes,
  ]
  weights_list <- vector("list", length(candidate_codes))
  n_dropped_origin <- 0L
  n_dropped_migration <- 0L

  for (i in seq_along(candidate_codes)) {
    sc <- candidate_codes[i]
    spp_stats <- region_stats[region_stats$species_code == sc, ]

    # get weights from regional stats table based on proportion population
    origin_rows <- spp_stats[spp_stats$season %in% origin_seasons, ]
    if (nrow(origin_rows) == 0L) {
      n_dropped_origin <- n_dropped_origin + 1L
      next
    }
    max_idx <- which.max(origin_rows$total_pop_percent)
    weight <- origin_rows$total_pop_percent[max_idx]
    max_origin_season <- origin_rows$season[max_idx]

    # drop species with weights below threshold
    if (is.na(weight) || weight < min_origin_threshold || weight == 0) {
      n_dropped_origin <- n_dropped_origin + 1L
      next
    }

    # drop species that don't sufficiently migrate out of origin region
    dest_rows <- spp_stats[spp_stats$season == destination_season, ]
    # missing row means prop pop = 0 for this region/season (rows with 0 are
    # omitted from the table), so treat as fully departed
    if (nrow(dest_rows) == 0L || is.na(dest_rows$total_pop_percent[1L])) {
      dest_max_week <- 0
    } else {
      dest_max_week <- dest_rows$total_pop_percent[1L]
    }

    rel_pct_decrease <- 1 - dest_max_week / weight
    if (rel_pct_decrease < min_migration_threshold) {
      n_dropped_migration <- n_dropped_migration + 1L
      next
    }

    spp_name_row <- runs_sub[runs_sub$species_code == sc, ]
    if (nrow(spp_name_row) > 0) {
      spp_name <- spp_name_row$common_name[1L]
    } else {
      spp_name <- NA_character_
    }

    weights_list[[i]] <- data.frame(
      species_code = sc,
      common_name = spp_name,
      origin_season = max_origin_season,
      weight = weight,
      rel_pct_decrease = rel_pct_decrease,
      stringsAsFactors = FALSE
    )
  }

  weights_df <- dplyr::bind_rows(weights_list)
  row.names(weights_df) <- NULL

  return(list(
    weights_df = weights_df,
    n_dropped_origin = n_dropped_origin + n_dropped_ineligible,
    n_dropped_migration = n_dropped_migration
  ))
}


stewardship_weights_polygon <- function(
  species,
  origin_region,
  origin_seasons,
  destination_season,
  min_origin_threshold,
  min_migration_threshold,
  resolution,
  path,
  force,
  show_progress
) {
  stopifnot(is.character(species), length(species) >= 1)
  stopifnot(inherits(origin_region, c("sf", "SpatVector")))
  stopifnot(nrow(origin_region) >= 1)

  species_codes <- get_species(species, strict = TRUE)

  check_stewardship_species(species_codes, explicit = TRUE)

  prop_pop_download_check(species_codes, "seasonal", resolution, path)

  if (inherits(origin_region, "sf")) {
    origin_vect <- terra::vect(origin_region)
  } else {
    origin_vect <- origin_region
  }
  if (nrow(origin_vect) > 1) {
    origin_vect <- terra::aggregate(origin_vect)
  }

  runs_sub <- ebirdst::ebirdst_runs[
    ebirdst::ebirdst_runs$species_code %in% species_codes,
  ]
  raster_crs <- NULL
  weights_list <- vector("list", length(species_codes))
  n_dropped_origin <- 0L
  n_dropped_migration <- 0L

  if (show_progress) {
    message("Calculating weights for ", length(species_codes), " species")
  }
  for (i in seq_along(species_codes)) {
    sc <- species_codes[i]

    if (show_progress) {
      message(sprintf("  [%d/%d] %s", i, length(species_codes), sc))
    }

    r <- load_raster(
      sc,
      product = "proportion-population",
      period = "seasonal",
      resolution = resolution,
      path = path,
      force = force,
      show_progress = FALSE
    )

    if (is.null(raster_crs)) {
      raster_crs <- terra::crs(r)
      origin_vect_proj <- terra::project(origin_vect, raster_crs)
    }

    # extract seasonal proportion of population within origin region
    seasonal_prop_pop <- terra::extract(
      r,
      origin_vect_proj,
      fun = "sum",
      na.rm = TRUE,
      weights = TRUE,
      ID = FALSE
    )
    seasonal_prop_pop <- stats::setNames(
      as.numeric(seasonal_prop_pop),
      names(seasonal_prop_pop)
    )

    # calculate weight based on origin proportion of population
    present_origin_seasons <- intersect(origin_seasons, names(r))
    if (length(present_origin_seasons) == 0) {
      n_dropped_origin <- n_dropped_origin + 1L
      next
    }
    season_weights <- seasonal_prop_pop[present_origin_seasons]
    max_idx <- which.max(season_weights)
    weight <- unname(season_weights[max_idx])
    max_origin_season <- present_origin_seasons[max_idx]

    # drop species with weights below threshold
    if (weight < min_origin_threshold || weight == 0) {
      n_dropped_origin <- n_dropped_origin + 1L
      next
    }

    # drop species that don't sufficiently migrate out of origin region
    if (!(destination_season %in% names(seasonal_prop_pop))) {
      n_dropped_migration <- n_dropped_migration + 1L
      next
    }
    dest_prop_in_region <- seasonal_prop_pop[destination_season]
    rel_pct_decrease <- 1 - dest_prop_in_region / weight

    if (rel_pct_decrease < min_migration_threshold) {
      n_dropped_migration <- n_dropped_migration + 1L
      next
    }

    spp_name <- runs_sub$common_name[runs_sub$species_code == sc]
    if (length(spp_name) == 0) {
      spp_name <- NA_character_
    }
    weights_list[[i]] <- data.frame(
      species_code = sc,
      common_name = spp_name,
      origin_season = max_origin_season,
      weight = weight,
      rel_pct_decrease = rel_pct_decrease,
      stringsAsFactors = FALSE
    )
  }

  weights_df <- dplyr::bind_rows(weights_list)
  row.names(weights_df) <- NULL

  return(list(
    weights_df = weights_df,
    n_dropped_origin = n_dropped_origin,
    n_dropped_migration = n_dropped_migration
  ))
}

stewardship_weights_table <- function(species_weights, species) {
  stopifnot(is.data.frame(species_weights))
  if (!("species_code" %in% names(species_weights))) {
    stop("species_weights must have a 'species_code' column.")
  }
  if (!("weight" %in% names(species_weights))) {
    stop("species_weights must have a 'weight' column.")
  }
  stopifnot(is.numeric(species_weights$weight))
  if (any(!is.finite(species_weights$weight) | species_weights$weight < 0)) {
    stop(
      "All values in species_weights$weight must be finite and non-negative."
    )
  }

  resolved <- get_species(species_weights$species_code, strict = TRUE)
  species_weights$species_code <- resolved

  if (anyDuplicated(species_weights$species_code) > 0) {
    stop(
      "species_weights contains duplicate species_code values. ",
      "Remove duplicates and try again."
    )
  }

  if (!is.null(species)) {
    stopifnot(is.character(species), length(species) >= 1)
    sub_resolved <- get_species(species, strict = TRUE)
    species_weights <- species_weights[
      species_weights$species_code %in% sub_resolved,
    ]
  }

  if (nrow(species_weights) == 0L) {
    stop("No species remain after subsetting species_weights.")
  }

  check_stewardship_species(species_weights$species_code, explicit = TRUE)

  runs_sub <- ebirdst::ebirdst_runs[
    ebirdst::ebirdst_runs$species_code %in% species_weights$species_code,
  ]
  common_names <- runs_sub$common_name[
    match(species_weights$species_code, runs_sub$species_code)
  ]

  weights_df <- data.frame(
    species_code = species_weights$species_code,
    common_name = common_names,
    origin_season = NA_character_,
    weight = species_weights$weight,
    rel_pct_decrease = NA_real_,
    stringsAsFactors = FALSE
  )
  row.names(weights_df) <- NULL

  return(list(
    weights_df = weights_df,
    n_dropped_origin = 0L,
    n_dropped_migration = 0L
  ))
}


normalize_destination_mask <- function(destination_mask, crs) {
  if (is.null(destination_mask)) {
    return(list(crop_ext = NULL, mask_vect = NULL))
  }

  if (inherits(destination_mask, "bbox")) {
    destination_mask <- sf::st_as_sfc(destination_mask)
  }

  if (inherits(destination_mask, "sfc")) {
    destination_mask <- sf::st_sf(geometry = destination_mask)
  }

  if (inherits(destination_mask, "sf")) {
    mask_vect <- terra::project(terra::vect(destination_mask), crs)
    return(list(
      crop_ext = terra::ext(mask_vect),
      mask_vect = mask_vect
    ))
  }

  if (inherits(destination_mask, "SpatVector")) {
    mask_vect <- terra::project(destination_mask, crs)
    return(list(
      crop_ext = terra::ext(mask_vect),
      mask_vect = mask_vect
    ))
  }

  if (inherits(destination_mask, "SpatExtent")) {
    return(list(
      crop_ext = destination_mask,
      mask_vect = NULL
    ))
  }

  stop(
    "destination_mask must be NULL, sf, SpatVector, SpatExtent, ",
    "or a bounding box from sf::st_bbox()."
  )
}
