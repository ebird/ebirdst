# clean up
unlink(list.files("man", full.names = TRUE))

# rebuild docs and install
devtools::document()
devtools::install_local(force = TRUE)

# local tests
devtools::test()
tools:::.check_package_datasets(".")

# vignettes, readme, site
devtools::clean_vignettes()
pkgdown::clean_site()
Sys.setenv(BUILD_VIGNETTES = TRUE)
rmarkdown::render("README.Rmd")
unlink("README.html")
pkgdown::build_site()
Sys.unsetenv("BUILD_VIGNETTES")

# local checks
devtools::check()

# checks
devtools::check_win_devel()
devtools::check_win_release()
