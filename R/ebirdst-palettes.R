#' eBird Status and Trends color palettes for mapping
#'
#' Generate the color palettes used for the eBird Status and Trends relative
#' abundance and trends maps.
#'
#' @param n integer; the number of colors to be in the palette.
#' @param type character; the type of color palette: "weekly" for the weekly
#'   relative abundance, "trends" for trends color palette, and a season name for
#'   the seasonal relative abundance. Note that for trends a diverging palette
#'   is returned, while all other palettes are sequential.
#'
#' @return A character vector of hex color codes.
#' @export
#'
#' @examples
#' # breeding season color palette
#' ebirdst_palettes(10, type = "breeding")
ebirdst_palettes <- function(n, type = c("weekly",
                                         "breeding", "nonbreeding",
                                         "migration",
                                         "prebreeding_migration",
                                         "postbreeding_migration",
                                         "year_round",
                                         "trends")) {
  stopifnot(is.numeric(n), length(n) == 1, n >= 1)
  type <- match.arg(type)

  # set base color by season
  col_zero <- "#e6e6e6"
  if (type == "weekly") {
    plsm <- rev(viridisLite::plasma(n - 1, end = 0.9))
    plsm <- stringr::str_remove(plsm, "FF$")
    gry <- grDevices::colorRampPalette(c(col_zero, plsm[1]))
    return(c(gry(4)[2], plsm))
  } else if (type == "trends") {
    trend_reds <- RColorBrewer::brewer.pal(9, "Reds")[4:7]
    trend_blues <- RColorBrewer::brewer.pal(9, "Blues")[4:7]
    trend_cols <- c(rev(trend_reds), "#ffffff", trend_blues)
    pal_fun <- grDevices::colorRampPalette(trend_cols)
    return(pal_fun(n))
  }else if (type == "breeding") {
    base_col <- "#cc503e"
  } else if (type == "nonbreeding") {
    base_col <- "#1d6996"
  } else if (type %in% c("migration", "postbreeding_migration")) {
    base_col <- "#edad08"
  } else if (type == "prebreeding_migration") {
    base_col <- "#73af48"
  } else if (type == "year_round") {
    base_col <- "#6f4070"
  } else {
    stop("Invalid season.")
  }

  # seasonal palettes
  gry <- grDevices::colorRampPalette(c(col_zero, base_col))
  mid <- grDevices::colorRampPalette(c(gry(5)[2], base_col))
  black <- grDevices::colorRampPalette(c(base_col, "#000000"))
  pal <- grDevices::colorRampPalette(c(gry(5)[2], mid(9)[5], base_col,
                                       black(5)[2]))
  return(pal(n))
}
