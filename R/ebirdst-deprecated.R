
#' @title Deprecated functions in package \pkg{ebirdst}.
#' @description The functions listed below are deprecated and support for them
#' will eventually be dropped.
#' Help pages for deprecated functions are
#'   available at \code{help("<function>-deprecated")}.
#' @name ebirdst-deprecated
#' @keywords internal
NULL


#' eBird Status and Trends color palettes for mapping
#'
#' This deprecated function has been replaced by \code{\link{ebirdst_palettes}}.
#' Both functions generate color palettes used for the eBird Status and Trends
#' relative abundance maps.
#'
#' @param n integer; the number of colors to be in the palette.
#' @param season character; the season to generate colors for or "weekly" to
#'   get the color palette used in the weekly abundance animations.
#'
#' @return A character vector of hex color codes.
#' @usage abundance_palette(n,
#'                         season = c("weekly", "breeding",
#'                                    "nonbreeding",
#'                                    "migration",
#'                                    "prebreeding_migration",
#'                                    "postbreeding_migration",
#'                                    "year_round"))
#' @name abundance_palette-deprecated
#' @seealso \code{\link{ebirdst_palettes}} \code{\link{ebirdst-deprecated}}
#' @keywords internal
NULL

#' @rdname ebirdst-deprecated
#' @section \code{abundance_palette}:
#' For \code{abundance_palette}, use \code{\link{ebirdst_palettes}}
#' @export
abundance_palette <- function(n,
                              season = c("weekly", "breeding",
                                         "nonbreeding",
                                         "migration",
                                         "prebreeding_migration",
                                         "postbreeding_migration",
                                         "year_round")) {
  .Deprecated(new = "ebirdst_palettes", package = "ebirdst")
  ebirdst_palettes(n = n, type = season)
}
