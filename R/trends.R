#' Convert eBird Trends Data Products to raster format
#'
#' The eBird trends data are stored in a tabular format, where each row gives
#' the trend estimate for a single cell in a 27km x 27km equal area grid. For
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
#' # download the trends data
#' download_ebird_trends()
#'
#' # load breeding season trends for sage thrasher and sagebrush sparrow
#' trends <- load_trends("sagthr_breeding")
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

  # check that only one species data is present
  n_runs <- nrow(dplyr::distinct(trends, .data$species_code, .data$season))
  if (n_runs > 1) {
    stop("There appears to be data for multiple species/season combinations ",
         "in the trends data frame. Only one trends run can be rasterized ",
         "at a time.")
  }

  # raster tempalte
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

#' eBird Trends color palette for mapping
#'
#' A diverging red to blue color palette for mapping eBird Trends. Reds are
#' typically assigned to declining trends and blues for increasing trends.
#'
#' @inheritParams abundance_palette
#'
#' @return A character vector of hex color codes.
#' @export
#'
#' @examples
#' trends_palette(9)
trends_palette <- function(n) {
  stopifnot(is.numeric(n), length(n) == 1, n >= 1)
  trend_reds <- RColorBrewer::brewer.pal(9, "Reds")[3:7]
  trend_blues <- RColorBrewer::brewer.pal(9, "Blues")[3:7]
  trend_cols <- c(rev(trend_reds), "#cfcfcf", trend_blues)
  pal_fun <- grDevices::colorRampPalette(trend_cols)
  return(pal_fun(n))
}
abundance_palette
# internal functions ----

trends_raster_template <- function() {
  e <- terra::ext(c(xmin = -20015109.354, xmax = 20036111.0830769,
                    ymin = -6684911.11603599, ymax = 10007554.677))
  crs <- terra::crs(paste("+proj=sinu +lon_0=0 +x_0=0 +y_0=0",
                          "+R=6371007.181 +units=m +no_defs"))
  terra::rast(e, crs = crs, nrows = 626L, ncols = 1502L)
}
