# clean up
unlink(list.files("man", full.names = TRUE))

# rebuild docs and install
devtools::document()
pak::local_install(ask = FALSE, dependencies = TRUE)

# local tests
devtools::test()
tools:::.check_package_datasets(".")

# vignettes, readme, site
devtools::clean_vignettes()
pkgdown::clean_site()
Sys.setenv(BUILD_VIGNETTES = TRUE)
devtools::build_readme()
pkgdown::build_site()
Sys.unsetenv("BUILD_VIGNETTES")

# local checks
devtools::check()

# checks
devtools::check_win_devel()
devtools::check_win_release()
