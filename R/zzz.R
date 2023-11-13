.onAttach <- function(libname, pkgname) {
  v <- ebirdst_version()
  vy <- v[["version_year"]]
  ry <- v[["release_year"]]
  aed <- format(v[["access_end_date"]], "%B %Y")
  citation <- stringr::str_glue(
    "Please cite the eBird Status & Trends data using: ",
    "Fink, D., T. Auer, A. Johnston, M. Strimas-Mackey, S. Ligocki, O. Robinson, ",
    "W. Hochachka, L. Jaromczyk, C. Crowley, K. Dunham, A. Stillman, I. Davies, ",
    "A. Rodewald, V. Ruiz-Gutierrez, C. Wood. {ry}.",
    "eBird Status and Trends, Data Version: {vy}; Released: {ry}. Cornell Lab of",
    "Ornithology, Ithaca, New York. https://doi.org/10.2173/ebirdst.{vy}",
    .sep = "\n  ", .trim = FALSE)
  m <- stringr::str_glue(
    "{citation}",
    "",
    "This version of the package provides access to the {vy} version of the eBird",
    "Status and Trends Data Products. Access to the {vy} data will be provided ",
    "until {aed} when it will be replaced by the {vy + 1} data. At that ",
    "point, you will be required to update this R package and transition to using ",
    "the new data.",
    .sep = "\n", .trim = FALSE
  )
  packageStartupMessage(as.character(m))
}
