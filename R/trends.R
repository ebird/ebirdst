#' Convert eBird Trends Data Products to raster format
#'
#' The eBird trends data are stored in a tabular format, where each row gives
#' the trend estimate for a single cell in a 27 km x 27 km equal area grid. For
#' many applications, an explicitly spatial format is more useful. This function
#' uses the cell center coordinates to convert the tabular trend estimates to
#' raster format in `terra` [SpatRaster][terra::SpatRaster] format.
#'
#' @param trends data frame; trends data for a single species as returned by
#'   [load_trends()].
#' @param layers character; column names in the trends data frame to rasterize.
#'   These columns will become layers in the raster that is created.
#' @param trim logical; flag indicating if the returned raster should be trimmed
#'   to remove outer rows and columns that are NA. If `trim = FALSE` the
#'   returned raster will have a global extent, which can be useful if rasters
#'   will be combined across species with different ranges.
#'
#' @return A [SpatRaster][terra::SpatRaster] object.
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
#' # rasterize percent per year trend
#' rasterize_trends(trends, "abd_ppy")
#' }
rasterize_trends <- function(trends,
                             layers = c("abd_ppy",
                                        "abd_ppy_lower",
                                        "abd_ppy_upper"),
                             trim = TRUE) {
  stopifnot(is_flag(trim))
  stopifnot(is.character(layers), !is.na(layers), length(layers) > 0)
  valid_layers <- c("abd",
                    "abd_ppy", "abd_ppy_lower", "abd_ppy_upper",
                    "abd_ppy_nonzero",
                    "abd_trend", "abd_trend_lower", "abd_trend_upper")
  if (any(!layers %in% valid_layers)) {
    stop("All layers to be converted must be one of:\n",
         paste(valid_layers, collapse = ", "))
  }

  stopifnot(is.data.frame(trends))
  required_cols <- c("species_code", "season",
                     "srd_id", "longitude", "latitude",
                     layers)
  if (any(!required_cols %in% names(trends))) {
    missing <- required_cols[!required_cols %in% names(trends)]
    stop("Input trends data frame must have the following columns:\n",
         paste(missing, collapse = ", "))
  }

  # check that we don't have fold-level estimates
  if ("fold" %in% names(trends)) {
    stop("Fold-level trend estimates cannot be rasterized. Try loading the ",
         "trends with `fold_estimates = FALSE`.")
  }

  # check that only one species data is present
  n_runs <- nrow(dplyr::distinct(trends, .data$species_code, .data$season))
  if (n_runs > 1) {
    stop("There appears to be data for multiple species/season combinations ",
         "in the trends data frame. Only one trends run can be rasterized ",
         "at a time.")
  }

  # check for duplicate srd cell estimates
  if (nrow(trends) != dplyr::n_distinct(trends$srd_id)) {
    stop("There are multiple rows for some cell estimates. Check that there ",
         "are not multiple rows for the same srd_id value.")
  }

  # raster template
  r <- trends_raster_template()

  # convert to vector format
  subset_cols <- c("latitude", "longitude", layers)
  trends_ss <- dplyr::select(trends, dplyr::all_of(subset_cols))
  v <- terra::vect(trends_ss,
                   geom = c("longitude", "latitude"),
                   crs = "epsg:4326")
  v <- terra::project(v, terra::crs(r))
  rm(trends, trends_ss)

  # rasterize
  # check terra version
  if (utils::packageVersion("terra") >= "1.7-3") {
    trends_raster <- terra::rasterize(v, r, field = layers)
  } else {
    trends_raster <- list()
    for (l in layers) {
      trends_raster[[l]] <- terra::rasterize(v, r, field = l)
    }
    trends_raster <- terra::rast(trends_raster)
  }
  names(trends_raster) <- layers

  if (isTRUE(trim)) {
    trends_raster <- terra::trim(trends_raster)
  }
  return(trends_raster)
}


#' Convert Trends Data Products to points or circles
#'
#' The eBird trends data are stored in a tabular format, where each row gives
#' the trend estimate for a single cell in a 27 km x 27 km equal area grid. For
#' many applications, an explicitly spatial format is more useful. This function
#' uses the cell center coordinates to convert the tabular trend estimates to
#' points or circles in [sf][sf::sf] format. Trends can be converted to points
#' or to circles with areas roughly proportional to the relative abundance
#' within that 27 km grid cell. These abundance-scaled circles are what is used
#' to produce the trends maps on the eBird Status and Trends website.
#'
#' @param trends data frame; trends data for a single species as returned by
#'   [load_trends()].
#' @param output character; "points" outputs spatial points while "circles"
#'   outputs circles with areas roughly proportional to the relative abundance
#'   within that 27 km grid cell.
#' @param crs character or `sf` [crs][sf::st_crs] object; coordinate reference
#'   system to output the results in. For points, unprojected latitude-longitude
#'   coordinates (the default) are most typical, while for circles use whatever
#'   equal area CRS you intend to use when mapping the data otherwise the
#'   "circles" will appear skewed.
#'
#' @returns Vetorized trends data as an [sf][sf::sf] object.
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
#' # vectorize as points
#' vectorize_trends(trends, "points")
#' # vectorize as circles
#' vectorize_trends(trends, "circles", crs = "+proj=eqearth")
#' }
vectorize_trends <- function(trends,
                             output = c("circles", "points"),
                             crs = 4326) {
  req_cols <- c("species_code", "latitude", "longitude", "abd")
  stopifnot(is.data.frame(trends), req_cols %in% names(trends))
  output <- match.arg(output)
  if (is.character(crs) || is.integer(crs) || is.numeric(crs)) {
    stopifnot(length(crs) == 1, !is.na(crs))
    crs <- sf::st_crs(crs)
  }
  stopifnot(inherits(crs, "crs"))

  # convert trends to spatial points and project
  trends_pts <- sf::st_as_sf(trends,
                             coords = c("longitude", "latitude"),
                             crs = 4326)
  trends_pts <- sf::st_transform(trends_pts, crs = crs)
  if (output == "points") {
    return(trends_pts)
  }

  srd_template <- trends_raster_template()
  radius_range <- c(min(terra::res(srd_template)) / 6,
                    0.99 * min(terra::res(srd_template)) / 2)

  # assign radii based on abundance
  trends_pts <- split(trends_pts, trends_pts$species_code)
  for (i in seq_along(trends_pts)) {
    abd <- trends_pts[[i]][["abd"]]
    abd_bins <- unique(stats::quantile(abd[abd > 0], seq(0, 1, by = 0.05)))

    if (length(abd_bins) == 1) {
      radii <- rep(radius_range[2], length.out = length(abd))
    } else {
      if (length(abd_bins) == 2) {
        abd_bins <- c(abd_bins[1], mean(abd_bins), abd_bins[2])
      }
      midpoint_radius <- sqrt(abd_bins[-length(abd_bins)] + diff(abd_bins) / 2)
      abd_bins[1] <- 0
      circle_rad <- scales::rescale(midpoint_radius, to = radius_range)
      trends_pts[[i]][["radii"]] <- categorize(abd, abd_bins, circle_rad)
    }
  }
  trends_pts <- dplyr::bind_rows(trends_pts)

  # buffer based on radius
  sf::st_buffer(trends_pts, dist = trends_pts$radii)
}


#' Convert percent per year trend to cumulative trend
#'
#' @param x numeric; percent per year trend on the 0-100 scale rather than the
#'   0-1 scale.
#' @param n_years integer; number of years.
#'
#' @return A numeric vector of the same length as `x` that contains the
#'   cumulative trend resulting from `n_years` years of compounding annual
#'   trend.
#' @export
#'
#' @examples
#' ppy_trend <- runif(100, min = -100, 100)
#' cumulative_trend <- convert_ppy_to_cumulative(ppy_trend, n_years = 5)
#' cbind(ppy_trend, cumulative_trend)
convert_ppy_to_cumulative <- function(x, n_years) {
  stopifnot(is.numeric(x), is_count(n_years))
  100 * ((1 + x / 100)^n_years - 1)
}


# internal functions ----

trends_raster_template <- function() {
  e <- terra::ext(c(xmin = -20015109.354, xmax = 20036111.0830769,
                    ymin = -6684911.11603599, ymax = 10007554.677))
  crs <- terra::crs(paste("+proj=sinu +lon_0=0 +x_0=0 +y_0=0",
                          "+R=6371007.181 +units=m +no_defs"))
  terra::rast(e, crs = crs, nrows = 626L, ncols = 1502L)
}

categorize <- function (x, breaks, labels) {
  stopifnot(is.numeric(x), is.numeric(breaks),
            is.numeric(labels) || is.character(labels),
            length(labels) == length(breaks) - 1)
  y <- cut(x, breaks)
  lvl <- levels(y)
  labels[match(y, lvl)]
}
