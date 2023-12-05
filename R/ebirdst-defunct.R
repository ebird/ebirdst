## ebirdst defunct functions
#' @title Defunct functions in package \pkg{ebirdst}.
#' @description The functions listed below are defunct and no longer supported.
#' Calling them will result in an error.
#'
#' When possible alternative functions are suggested.
#'
#'  Many of them supported stixles which were infrequently used and were
#'  dropped from \pkg{ebirdst} with the 2022 data release.
#' @param ... All arguments are now ignored.
#' @name ebirdst-defunct
#' @keywords internal
NULL

#' @rdname ebirdst-defunct
#' @export
ebirdst_download <- function(species,
                             path = ebirdst_data_dir(),
                             tifs_only = TRUE,
                             force = FALSE,
                             show_progress = TRUE,
                             pattern = NULL,
                             dry_run = FALSE) {
  .Defunct("ebirdst_download_status", package = "ebirdst")
}

#' @rdname ebirdst-defunct
#' @export
ebirdst_extent <- function(x, t, ...) {
  .Defunct(package = "ebirdst")
}

#' @rdname ebirdst-defunct
#' @export
ebirdst_habitat <- function(path, ext, data = NULL,
                            stationary_associations = FALSE) {
  .Defunct(package = "ebirdst")
}

#' @rdname ebirdst-defunct
ebirdst_ppms <- function(path, ext, es_cutoff, pat_cutoff) {
  .Defunct(package = "ebirdst")
}

#' @rdname ebirdst-defunct
#' @export
ebirdst_ppms_ts <- function(ath, ext, summarize_by = c("weeks", "months"), ...){
  .Defunct(package = "ebirdst")
}

#' @rdname ebirdst-defunct
#' @export
ebirdst_subset <- function(x, crs) {
  .Defunct(package = "ebirdst")
}

#' @rdname ebirdst-defunct
#' @export
load_pds <- function(path,
                     ext,
                     model = c("occurrence", "count"),
                     return_sf = FALSE) {
  .Defunct(package = "ebirdst")
}

#' @rdname ebirdst-defunct
#' @export
load_pis <- function(path, ext,
                     model = c("occurrence", "count"), return_sf = FALSE){
  .Defunct(package = "ebirdst")
}

#' @rdname ebirdst-defunct
#' @export
load_predictions <- function(path, return_sf = FALSE) {
  .Defunct(package = "ebirdst")
}

#' @rdname ebirdst-defunct
#' @export
parse_raster_dates <- function(x){
  .Defunct(package = "ebirdst")
}

#' @rdname ebirdst-defunct
#' @export
load_stixels <- function(path, ext, return_sf = FALSE) {
  .Defunct(package = "ebirdst")
}

#' @rdname ebirdst-defunct
#' @export
project_extent <- function(x, crs) {
  .Defunct(package = "ebirdst")
}

#' @rdname ebirdst-defunct
#' @export
plot_pds <- function(path, ext, summarize_by = c("weeks", "months"), ...) {
  .Defunct(package = "ebirdst")
}

#' @rdname ebirdst-defunct
#' @export
plot_pis <- function(pis,
                     ext,
                     by_cover_class = TRUE,
                     n_top_pred = 15,
                     pretty_names = TRUE,
                     plot = TRUE) {
  .Defunct(package = "ebirdst")
}

#' @rdname ebirdst-defunct
#' @export
stixelize <- function(x){
  .Defunct(package = "ebirdst")
}


