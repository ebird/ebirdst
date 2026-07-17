# rebuild docs and install
devtools::document()
pak::local_install(ask = FALSE, dependencies = TRUE)

# local tests
devtools::test()
tools:::.check_package_datasets(".")

# vignettes, readme, site
pkgdown::clean_site()
devtools::build_readme()
pkgdown::build_site()
# pkgdown turns every root *.md file into a page, ignoring .Rbuildignore;
# remove the local-only CLAUDE.md page it generates
unlink("docs/CLAUDE.html")

# local checks
devtools::check()

# checks
devtools::check_win_devel()
devtools::check_win_release()
