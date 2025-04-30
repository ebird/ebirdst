.onAttach <- function(libname, pkgname) {
  v <- ebirdst_version()
  svy <- v[["status_version_year"]]
  tvy <- v[["trends_version_year"]]

  status_citation <- paste(
    "Cite the eBird Status Data Products using: ",
    "Fink, D., T. Auer, A. Johnston, M. Strimas-Mackey, S. Ligocki, O. Robinson,",
    "W. Hochachka, L. Jaromczyk, C. Crowley, K. Dunham, A. Stillman, C. Davis,",
    "M. Stokowski, P. Sharma, V. Pantoja, D. Burgin, P. Crowe, M. Bell, S. Ray,",
    "I. Davies, V. Ruiz-Gutierrez, C. Wood, A. Rodewald. 2024. eBird Status and",
    "Trends, Data Version: 2023; Released: 2025. Cornell Lab of Ornithology, Ithaca,",
    "New York. https://doi.org/10.2173/WZTW8903", sep = "\n  ")

  trends_citation <- paste(
    "Cite the eBird Trends Data Products using: ",
    "Fink, D., T. Auer, A. Johnston, M. Strimas-Mackey, S. Ligocki, O. Robinson, ",
    "W. Hochachka, L. Jaromczyk, C. Crowley, K. Dunham, A. Stillman, I. Davies, ",
    "A. Rodewald, V. Ruiz-Gutierrez, C. Wood. 20. eBird Status and Trends, Data",
    "Version: 2022; Released: 2023. Cornell Lab of Ornithology, Ithaca, New York. ",
    "https://doi.org/10.2173/ebirdst.2022", sep = "\n  ")

  m <- stringr::str_glue(
    "This version of the ebirdst package provides access to the {svy} version of ",
    "the eBird Status Data Products and the {tvy} version of the eBird Trends Data ",
    "Products.\n\n",
    "{status_citation}", "\n\n", "{trends_citation}",
    .sep = "", .trim = FALSE
  )
  packageStartupMessage(as.character(m))
}
