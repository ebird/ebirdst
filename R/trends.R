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
  trends_raster <- terra::rasterize(v, r, field = layers)
  names(trends_raster) <- layers

  if (isTRUE(trim)) {
    trends_raster <- terra::trim(trends_raster)
  }
  return(trends_raster)
}


# internal functions ----

trends_raster_template <- function() {
  e <- terra::ext(c(xmin = -20015109.354, xmax = 20036111.0830769,
                    ymin = -6684911.11603599, ymax = 10007554.677))
  crs <- terra::crs(paste("+proj=sinu +lon_0=0 +x_0=0 +y_0=0",
                          "+R=6371007.181 +units=m +no_defs"))
  terra::rast(e, crs = crs, nrows = 626L, ncols = 1502L)
}
