#' Spatiotemporal grid sampling of observation data
#'
#' Sample observation data on a spacetime grid to reduce spatiotemporal bias.
#'
#' @param x data frame; observations to sample, including at least the columns
#'   defining the location in space and time. Additional columns can be included
#'   such as features that will later be used in model training.
#' @param coords character; names of the spatial and temporal coordinates. By
#'   default the spatial spatial coordinates should be `longitude` and
#'   `latitude`, and temporal coordinate should be `day_of_year`.
#' @param is_lonlat logical; if the points are in unprojected, lon-lat
#'   coordinates. In this case, the points will be projected to an equal area
#'   Eckert IV CRS prior to grid assignment.
#' @param res numeric; resolution of the spatiotemporal grid in the x, y, and
#'   time dimensions. Unprojected locations are projected to an equal area
#'   coordinate system prior to sampling, and resolution should therefore be
#'   provided in units of meters. The temporal resolution should be in the
#'   native units of the time coordinate in the input data frame, typically it
#'   will be a number of days.
#' @param jitter_grid logical; whether to jitter the location of the origin of
#'   the grid to introduce some randomness.
#' @param sample_size_per_cell integer; number of observations to sample from
#'   each grid cell.
#' @param cell_sample_prop proportion `(0-1]`; if less than 1, only this
#'   proportion of cells will be randomly selected for sampling.
#' @param keep_cell_id logical; whether to retain a unique cell identifier,
#'   stored in column named `.cell_id`.
#' @param grid_definition list defining the spatiotemporal sampling grid as
#'   returned by [assign_to_grid()] in the form of an attribute of the returned
#'   data frame.
#'
#' @details
#' [grid_sample_stratified()] performs stratified case control sampling,
#' independently sampling from strata defined by, for example, year and
#' detection/non-detection. Within each stratum, [grid_sample()] is used to
#' sample the observations on a spatiotemporal grid. In addition, if case
#' control sampling is turned on, the detections are oversampled to increase the
#' frequency of detections in the dataset.
#'
#' The sampling grid is defined, and assignment of locations to cells occurs, in
#' [assign_to_grid()]. Consult the help for that function for further details on
#' how the grid is generated and locations are assigned. Note that by providing
#' 2-element vectors to both `coords` and `res` the time component of the grid
#' can be ignored and spatial-only subsampling is performed.
#'
#' @return A data frame of the spatiotemporally sampled data.
#' @export
#'
#' @examples
#' set.seed(1)
#'
#' # generate some example observations
#' n_obs <- 10000
#' checklists <- data.frame(longitude = rnorm(n_obs, sd = 0.1),
#'                          latitude = rnorm(n_obs, sd = 0.1),
#'                          day_of_year = sample.int(28, n_obs, replace = TRUE),
#'                          year = NA_integer_,
#'                          obs = rpois(n_obs, lambda = 0.1),
#'                          forest_cover = runif(n_obs),
#'                          island = as.integer(runif(n_obs) > 0.95))
#' # add a year column, giving more data to recent years
#' checklists$year <- sample(seq(2016, 2020), size = n_obs, replace = TRUE,
#'                           prob = seq(0.3, 0.7, length.out = 5))
#' # create several rare islands
#' checklists$island[sample.int(nrow(checklists), 9)] <- 2:10
#'
#' # basic spatiotemporal grid sampling
#' sampled <- grid_sample(checklists)
#'
#' # plot original data and grid sampled data
#' par(mar = c(0, 0, 0, 0))
#' plot(checklists[, c("longitude", "latitude")],
#'      pch = 19, cex = 0.3, col = "#00000033",
#'      axes = FALSE)
#' points(sampled[, c("longitude", "latitude")],
#'        pch = 19, cex = 0.3, col = "red")
#'
#' # case control sampling stratified by year and island
#' # return a maximum of 1000 checklists
#' sampled_cc <- grid_sample_stratified(checklists, sample_by = "island",
#'                                      maximum_ss = 1000)
#'
#' # case control sampling increases the prevalence of detections
#' mean(checklists$obs > 0)
#' mean(sampled$obs > 0)
#' mean(sampled_cc$obs > 0)
#'
#' # stratifying by island ensures all levels are retained, even rare ones
#' table(checklists$island)
#' # normal grid sampling loses rare island levels
#' table(sampled$island)
#' # stratified grid sampling retain at least one observation from each level
#' table(sampled_cc$island)
grid_sample <- function(x,
                        coords = c("longitude", "latitude", "day_of_year"),
                        is_lonlat = TRUE,
                        res = c(3000, 3000, 7), jitter_grid = TRUE,
                        sample_size_per_cell = 1,
                        cell_sample_prop = 0.75,
                        keep_cell_id = FALSE,
                        grid_definition = NULL) {
  # input checks
  stopifnot(is.data.frame(x))
  stopifnot(is.character(coords), length(coords) == 3, all(!is.na(coords)),
            all(coords %in% names(x)))
  stopifnot(is_flag(is_lonlat))
  stopifnot(is.numeric(res), length(res) %in% c(2, 3), all(res > 0))
  stopifnot(is_flag(jitter_grid))
  stopifnot(is_count(sample_size_per_cell), sample_size_per_cell > 0)
  stopifnot(is.numeric(cell_sample_prop), length(cell_sample_prop) == 1,
            cell_sample_prop > 0, cell_sample_prop <=1)
  stopifnot(is_flag(keep_cell_id))

  # handle edge case of no observations
  if (nrow(x) == 0) {
    return(x)
  }

  # assign the observations to a grid
  if (!is.null(grid_definition)) {
    cells <- assign_to_grid(points = x,
                            grid_definition = grid_definition)
  } else {
    cells <- assign_to_grid(points = x, res = res,
                            coords = coords,
                            jitter_grid = jitter_grid,
                            is_lonlat = is_lonlat,
                            grid_definition = NULL)
  }

  # was this a spacetime grid or just a spatial grid
  if ("cell_xyt" %in% names(cells)) {
    cells <- cells[, "cell_xyt", drop = FALSE]
  } else if ("cell_xy" %in% names(cells)) {
    cells <- cells[, "cell_xy", drop = FALSE]
  } else {
    stop("Invalid cells generated by assign_to_grid()")
  }
  names(cells) <- "cell_id"
  cells[["row_number"]] <- seq_len(nrow(cells))

  if (any(is.na(cells[["cell_id"]]))) {
    stop("Sampling grid did not capture all points.")
  }

  # only select a subset of cell, keep at least 1
  if (cell_sample_prop < 1) {
    ids <- unique(cells[["cell_id"]])
    n_cells <- max(round(cell_sample_prop * length(ids)), 1)
    keed_ids <- safe_sample(ids, size = n_cells, replace = FALSE)
    cells <- cells[cells[["cell_id"]] %in% keed_ids, ]
    rm(ids, keed_ids)
  }

  # sample from each cell
  sampled <- tapply(cells[["row_number"]],
                    INDEX = cells[["cell_id"]],
                    FUN = safe_sample, size = sample_size_per_cell,
                    simplify = FALSE)
  sampled <- unique(do.call(c, sampled))

  # extract just the sampled rows of the original dataset
  x_sampled <- dplyr::tibble(x[sampled, ])
  if (keep_cell_id) {
    lookup <- match(sampled, cells[["row_number"]])
    x_sampled[[".cell_id"]] <- cells[["cell_id"]][lookup]
  }
  return(x_sampled)
}

#' @param unified_grid logical; whether a single, unified spatiotemporal
#'   sampling grid should be defined and used for all observations in `x` or a
#'   different grid should be used for each stratum.
#' @param by_year logical; whether the sampling should be done by year, i.e.
#'   sampling N observations per grid cell per year, rather than across years,
#'   i.e. N observations per grid cell regardless of year. If using sampling by
#'   year, the input data frame `x` must have a `year` column.
#' @param case_control logical; whether to apply case control sampling whereby
#'   presence and absence are sampled independently.
#' @param obs_column character; if `case_control = TRUE`, this is the name of
#'   the column in `x` that defines detection (`obs_column > 0`) and
#'   non-detection (`obs_column == 0`).
#' @param sample_by character; additional columns in `x` to stratify sampling
#'   by. For example, if a landscape has many small islands (defined by an
#'   `island` variable) and we wish to sample from each independently, use
#'   `sample_by = "island"`.
#' @param min_detection_probability proportion `[0-1)`; the minimum detection
#'   probability in the final dataset. If `case_control = TRUE`, and the
#'   proportion of detections in the grid sampled dataset is below this level,
#'   then additional detections will be added via grid sampling the detections
#'   from the input dataset until at least this proportion of detections appears
#'   in the final dataset. This will typically result in duplication of some
#'   observations in the final dataset. To turn this off this feature use
#'   `min_detection_probability = 0`.
#' @param maximum_ss integer; the maximum sample size in the final dataset. If
#'   the grid sampling yields more than this number of observations,
#'   `maximum_ss` observations will be selected randomly from the full set. Note
#'   that this subsampling will be performed in such a way that all levels of
#'   each strata will have at least one observation within the final dataset,
#'   and therefore it is not truly randomly sampling.
#' @param jitter_columns character; if detections are oversampled to achieve the
#'   minimum detection probability, some observations will be duplicated, and it
#'   can be desirable to slightly "jitter" the values of model training features
#'   for these duplicated observations. This argument defines the column names
#'   in `x` that will be jittered.
#' @param jitter_sd numeric; strength of the jittering in units of standard
#'   deviations, see `jitter_columns`.
#' @param ... additional arguments defining the spatiotemporal grid; passed to
#'   [grid_sample()].
#'
#' @rdname grid_sample
#' @export
grid_sample_stratified <- function(
    x,
    coords = c("longitude", "latitude", "day_of_year"),
    is_lonlat = TRUE,
    unified_grid = FALSE,
    keep_cell_id = FALSE,
    by_year = TRUE,
    case_control = TRUE, obs_column = "obs",
    sample_by = NULL,
    min_detection_probability = 0,
    maximum_ss = NULL,
    jitter_columns = NULL,
    jitter_sd = 0.1,
    ...) {

  # input checks
  stopifnot(is.data.frame(x))
  stopifnot(is_flag(is_lonlat), is_flag(unified_grid), is_flag(keep_cell_id))
  stopifnot(is_flag(case_control))
  if (case_control) {
    stopifnot(is.character(obs_column), length(obs_column) == 1,
              obs_column %in% names(x))
  }
  if (is.null(sample_by)) {
    sample_by <- character()
  }
  stopifnot(is.character(sample_by), all(sample_by %in% names(x)))
  stopifnot(is_flag(by_year))
  if (by_year) {
    stopifnot("year" %in% names(x))
    sample_by <- unique(c("year", sample_by))
  }
  if (!is.null(maximum_ss)) {
    stopifnot(is_count(maximum_ss), maximum_ss > 0)
  }
  stopifnot(is.numeric(min_detection_probability),
            length(min_detection_probability) == 1,
            min_detection_probability >= 0,
            min_detection_probability < 1)

  if (keep_cell_id && !unified_grid) {
    warning("Cell IDs can only be returned if unified_grid = TRUE since ",
            "different grids are used for each stratum otherwise. Setting ",
            "keep_cell_id = FALSE")
    keep_cell_id <- FALSE
  }

  if (!is.null(jitter_columns)) {
    stopifnot(is.character(jitter_columns),
              all(jitter_columns %in% names(x)))
    stopifnot(is.numeric(jitter_sd), length(jitter_sd) == 1, jitter_sd >= 0)

    # ensure columns to jitter are numeric
    col_types <- vapply(x, class, character(1))
    col_types <- col_types[jitter_columns]
    not_numeric <- names(col_types)[!col_types %in% c("numeric", "integer")]
    if (length(not_numeric) > 0) {
      stop("The following jitter columns must be of type numeric or integer: ",
           paste(not_numeric, collapse = ", "))
    }
  }

  # handle edge case of no observations
  if (nrow(x) == 0) {
    return(x)
  }

  # no strata defined
  if (!by_year && !case_control &&
      (is.null(sample_by) || length(sample_by) == 0)) {
    return(grid_sample(x, keep_cell_id = keep_cell_id, ...))
  }

  # subset to just the location and strata columns
  locs <- x[, coords, drop = FALSE]

  # project once now to avoid having to do it for every stratum
  if (is_lonlat) {
    xy <- project_equal_area(locs, coords = coords[1:2])
    # add time dimension
    if (length(coords) == 3) {
      xy[["t"]] <- locs[[coords[3]]]
    }
    locs <- xy
    coords <- c("x", "y", "t")
    is_lonlat <- FALSE
    rm(xy)
  }
  locs[[".row"]] <- seq_len(nrow(locs))

  # subset to just the location and strata columns
  locs <- dplyr::bind_cols(locs, x[, sample_by, drop = FALSE])

  # add presence/absence columns for case control sampling
  # also stratify by this new column
  if (case_control) {
    locs[[".detected"]] <- x[[obs_column]] > 0

    # turn off case control sampling if there is only one class
    if (all(locs[[".detected"]])) {
      message("No non-detections in dataset, setting case_control = FALSE")
      case_control <- FALSE
    } else if (all(!locs[[".detected"]])) {
      message("No detections in dataset, setting case_control = FALSE")
      case_control <- FALSE
    } else {
      sample_by <- c(sample_by, ".detected")
    }
  }

  # split all locs into strata
  locs_split <- split(locs, f = locs[, sample_by, drop = FALSE], drop = TRUE)

  # define the grid once for all observations
  grid_definition <- NULL
  if (unified_grid) {
    cells <- assign_to_grid(points = locs,
                            coords = coords,
                            is_lonlat = is_lonlat,
                            res = c(3000, 3000, 7),
                            jitter_grid = TRUE)
    grid_definition <- attr(cells, "grid_definition")
    rm(cells)
  }

  # sample from each stratum
  sampled <- lapply(locs_split, FUN = grid_sample,
                    coords = coords, is_lonlat = is_lonlat,
                    keep_cell_id = keep_cell_id,
                    grid_definition = grid_definition,
                    ...)
  rm(locs_split)
  sampled <- dplyr::bind_rows(sampled)

  # subsample to decrease sample size to maximum
  # TODO consider adding && nrow(sampled) > maximum_ss here
  if (!is.null(maximum_ss)) {
    sample_prop <- maximum_ss / nrow(sampled)
    if (case_control) {
      # case control sampling on: sample preserving detection probability
      # calculate detection probability for the grid sampled data
      det_prob <- mean(sampled[[".detected"]], na.rm = TRUE)

      # bias to detections if detection probability is below the minimum
      # otherwise preserve the detection probability
      target_det_prob <- max(det_prob, min_detection_probability)
      target_prop_det <- sample_prop * target_det_prob / det_prob
      target_prop_non <- sample_prop * (1 - target_det_prob) / (1 - det_prob)

      # sample detections
      s <- split(sampled, f = sampled[[".detected"]], drop = TRUE)
      sample_by_nodet <- setdiff(sample_by, ".detected")
      if (target_prop_det < 1 && "TRUE" %in% names(s)) {
        s[["TRUE"]] <- sample_stratify(s[["TRUE"]],
                                       prop = target_prop_det,
                                       sample_by = sample_by_nodet)
      }
      # sample non-detections
      if (target_prop_non < 1 && "FALSE" %in% names(s)) {
        s[["FALSE"]] <- sample_stratify(s[["FALSE"]],
                                        prop = target_prop_non,
                                        sample_by = sample_by_nodet)
      }
      sampled <- dplyr::bind_rows(s)
      rm(s)
    } else if (nrow(sampled) > maximum_ss) {
      sampled <- sample_stratify(sampled,
                                 prop = sample_prop,
                                 sample_by = sample_by)
    }
  }

  # if case control sampling is off, return here
  if (!case_control) {
    x_sampled <- x[sampled[[".row"]], ]
    if (keep_cell_id) {
      x_sampled[[".cell_id"]] <- sampled[[".cell_id"]]
    }
    return(x_sampled)
  }

  # oversample detections if detection probability is below minimum
  det_prob <- mean(sampled[[".detected"]], na.rm = TRUE)
  if (det_prob < min_detection_probability && det_prob > 0) {
    # just the detections
    locs_det <- locs[locs[[".detected"]], ]

    # split detections by year, but ignore other strata at this stage
    if (by_year) {
      locs_split <- split(locs_det, f = locs_det[["year"]], drop = TRUE)
    } else {
      locs_split <- list(locs_det)
    }
    rm(locs_det)

    # target probability
    p <- min_detection_probability
    n_non <- sum(!sampled[[".detected"]], na.rm = TRUE)

    # repeatedly add additional grid sampled positives, up to 25 times
    for (i in seq_len(25)) {
      # number of detections to add to reach the target
      n_det <- sum(sampled[[".detected"]], na.rm = TRUE)
      n_det_add <- round(n_non * p / (1 - p) - n_det)
      if (n_det_add <= 0) {
        break
      }

      # grid sample to add new detections
      add_dets <- lapply(locs_split, FUN = grid_sample,
                         coords = coords, is_lonlat = is_lonlat,
                         keep_cell_id = keep_cell_id,
                         grid_definition = grid_definition,
                         ...)
      add_dets <- dplyr::bind_rows(add_dets)
      # only take as many as required to reach the target
      if (nrow(add_dets) > n_det_add) {
        add_dets <- dplyr::slice_sample(add_dets, n = n_det_add)
      }
      sampled <- dplyr::bind_rows(sampled, add_dets)
    }
  }

  # if covariate jittering doesn't need to be performed, return here
  has_duplicates <- anyDuplicated(sampled[[".row"]]) > 0
  if (is.null(jitter_columns) || jitter_sd == 0 || !has_duplicates) {
    x_sampled <- x[sampled[[".row"]], ]
    if (keep_cell_id) {
      x_sampled[[".cell_id"]] <- sampled[[".cell_id"]]
    }
    return(x_sampled)
  }

  # identify duplicate records
  dupes <- duplicated(sampled[[".row"]])
  x_dupes <- x[sampled[[".row"]][dupes], ]
  x_sampled <- x[sampled[[".row"]][!dupes], ]
  if (keep_cell_id) {
    x_dupes[[".cell_id"]] <- sampled[[".cell_id"]][dupes]
    x_sampled[[".cell_id"]] <- sampled[[".cell_id"]][!dupes]
  }

  # estimate standard deviation from full dataset
  col_sd <- apply(x[, jitter_columns, drop = FALSE],
                  MARGIN = 2,
                  FUN = stats::sd, na.rm = TRUE)
  # scale standard deviation
  col_sd <- col_sd * jitter_sd

  # apply jittering to each column
  for (col in jitter_columns) {
    jitter_amt <- col_sd[[col]] * stats::runif(nrow(x_dupes), min = -1, max = 1)
    x_dupes[[col]] <- x_dupes[[col]] + jitter_amt
  }

  # combine and return jittered and unjittered data
  return(dplyr::tibble(dplyr::bind_rows(x_sampled, x_dupes)))
}


#' Assign points to a spacetime grid
#'
#' Given a set of points in space and (optionally) time, define a regular grid
#' with given dimensions, and return the grid cell index for each point.
#'
#' @param points data frame; points with spatial coordinates `x` and `y`, and an
#'   optional time coordinate `t`.
#' @param coords character; names of the spatial and temporal coordinates in the
#'   input dataframe. Only provide these names if you want to overwrite the
#'   default coordinate names: `c("x", "y", "t")` or `c("longitude", "latitude",
#'   "t")` if `is_lonlat = TRUE`.
#' @param is_lonlat logical; if the points are in unprojected, lon-lat
#'   coordinates. In this case, the input data frame should have columns
#'   `"longitude"` and `"latitude"` and the points will be projected to an equal
#'   area Eckert IV CRS prior to grid assignment.
#' @param res numeric; resolution of the grid in the `x`, `y`, and `t`
#'   dimensions, respectively. If only 2 dimensions are provided, a space only
#'   grid will be generated. The units of `res` are the same as the coordinates
#'   in the input data unless `is_lonlat` is true in which case the `x` and `y`
#'   resolution should be provided in meters.
#' @param jitter_grid logical; whether to jitter the location of the origin of
#'   the grid to introduce some randomness.
#' @param grid_definition list; object defining the grid via the `origin` and
#'   `resolution` components. To assign multiple sets of points to exactly the
#'   same grid, [assign_to_grid()] returns a data frame with a `grid_definition`
#'   attribute that can be passed to subsequent calls to [assign_to_grid()].
#'   `res` and `jitter` are ignored if `grid_definition` is provided.
#'
#' @return Data frame with the indices of the space-only and spacetime grid
#'   cells. This data frame will have a `grid_definition` attribute that can be
#'   used to reconstruct the grid.
#' @export
#'
#' @examples
#' set.seed(1)
#'
#' # generate some example points
#' points_xyt <- data.frame(x = runif(100), y = runif(100), t = rnorm(100))
#' # assign to grid
#' cells <- assign_to_grid(points_xyt, res = c(0.1, 0.1, 0.5))
#'
#' # assign a second set of points to the same grid
#' assign_to_grid(points_xyt, grid_definition = attr(cells, "grid_definition"))
#'
#' # assign lon-lat points to a 10km space-only grid
#' points_ll <- data.frame(longitude = runif(100, min = -180, max = 180),
#'                         latitude = runif(100, min = -90, max = 90))
#' assign_to_grid(points_ll, res = c(10000, 10000), is_lonlat = TRUE)
#'
#' # overwrite default coordinate names, 5km by 1 week grid
#' points_names <- data.frame(lon = runif(100, min = -180, max = 180),
#'                            lat = runif(100, min = -90, max = 90),
#'                            day = sample.int(365, size = 100))
#' assign_to_grid(points_names,
#'                res = c(5000, 5000, 7),
#'                coords = c("lon", "lat", "day"),
#'                is_lonlat = TRUE)
assign_to_grid <- function(points,
                           coords = NULL, is_lonlat = FALSE,
                           res, jitter_grid = TRUE,
                           grid_definition = NULL) {
  if (!is.null(grid_definition)) {
    stopifnot(all(c("res", "origin") %in% names(grid_definition)))
    res <- grid_definition[["res"]]
  } else {
    stopifnot(is_flag(jitter_grid))
  }
  stopifnot(is.numeric(res), length(res) %in% c(2, 3), all(res > 0))
  stopifnot(is_flag(is_lonlat))

  if (!is.null(coords)) {
    stopifnot(is.character(coords), length(coords) == length(res),
              all(!is.na(coords)), all(coords %in% names(points)))
    if (is_lonlat) {
      default_names <- c("longitude", "latitude", "t")
    } else {
      default_names <- c("x", "y", "t")
    }
    default_names <- default_names[seq_len(length(coords))]
    names(points)[match(coords, names(points))] <- default_names
  }

  # project to equal area if unprojected
  if (is_lonlat) {
    stopifnot(all(c("longitude", "latitude") %in% names(points)))
    coords <- project_equal_area(points, c("longitude", "latitude"))
    # add temporal columns
    if ("t" %in% names(points)) {
      coords[["t"]] <- points[["t"]]
    }
    points <- coords
    rm(coords)
  }

  # 2d or 3d grid?
  if (length(res) == 2) {
    stopifnot(all(c("x", "y") %in% names(points)))
    points <- points[, c("x", "y")]
    time_dim <- FALSE
  } else {
    stopifnot(all(c("x", "y", "t") %in% names(points)))
    points <- points[, c("x", "y", "t")]
    time_dim <- TRUE
  }

  # define grid based on lower left corner of grid
  if (is.null(grid_definition)) {
    ll <- apply(points, 2, min, na.rm = TRUE)
    # jitter
    if (jitter_grid) {
      jitter_amount <- stats::runif(length(res)) * res
      ll <- ll - jitter_amount
    }
  } else {
    ll <- grid_definition[["origin"]]
  }

  # assign indexes for space dimensions
  x_cell <- 1 + (points[["x"]] - ll[1]) %/% res[1]
  y_cell <- 1 + (points[["y"]] - ll[2]) %/% res[2]
  if (is.null(grid_definition)) {
    grid_dim <- c(x = max(x_cell), y = max(y_cell))
  } else {
    grid_dim <- grid_definition[["dim"]]
  }
  stopifnot(!is.na(x_cell), !is.na(y_cell))
  cells <- dplyr::tibble(cell_xy = paste0(x_cell, "-", y_cell))

  # assign spacetime index
  if (time_dim) {
    t_cell <- 1 + (points[["t"]] - ll[3]) %/% res[3]
    if (is.null(grid_definition)) {
      grid_dim <- c(grid_dim, t = max(t_cell))
    }
    stopifnot(!is.na(t_cell))
    cells[["cell_xyt"]] <- paste0(cells[["cell_xy"]], "-", t_cell)
  }

  # check that none of the points were outside the grid
  if (!is.null(grid_definition)) {
    if (any(x_cell > grid_dim[["x"]]) || any(y_cell > grid_dim[["y"]])) {
      stop("some points fall outside the provided spatial grid")
    }
    if (time_dim && any(t_cell > grid_dim[["t"]])) {
      stop("some points fall outside the provided temporal grid")
    }
  }

  # attach grid definition to return object
  attr(cells, "grid_definition") <- list(res = res, dim = grid_dim, origin = ll)
  return(cells)
}


# internal ----

# project to equal area
project_equal_area <- function(points, coords = c("longitude", "latitude")) {
  stopifnot(all(coords %in% names(points)))
  locs <- terra::vect(as.data.frame(points[, coords]),
                      geom = coords, crs = "epsg:4326")
  locs <- terra::project(locs, y = "+proj=eck4")
  locs <- as.data.frame(terra::crds(locs))
  return(stats::setNames(locs, c("x", "y")))
}

# safe version of sample() that works when length(x) == 1
safe_sample <- function(x, size, ...) {
  if (length(x) <= size || length(x) == 1) {
    return(x)
  }
  sample(x, size = size, ...)
}

sample_stratify <- function(x, prop, sample_by) {
  if (length(sample_by) == 0) {
    x <- list(x)
  } else {
    x <- split(x, f = x[, sample_by, drop = FALSE], drop = TRUE)
  }
  n <- vapply(x, FUN = nrow, FUN.VALUE = integer(1))
  # ensure at least one row is retained from each stratum
  size <- pmax(1, round(n * prop))
  # sample by group
  sampled <- mapply(FUN = dplyr::slice_sample, .data = x, n = size,
                    SIMPLIFY = FALSE)
  dplyr::bind_rows(sampled)
}
