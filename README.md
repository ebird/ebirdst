
<!-- README.md is generated from README.Rmd. Please edit that file -->

# eBird Status and Trends Data Products

<!-- badges: start -->

[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/ebird/ebirdst/workflows/R-CMD-check/badge.svg)](https://github.com/ebird/ebirdst/actions)
[![CRAN
status](https://www.r-pkg.org/badges/version/ebirdst)](https://cran.r-project.org/package=ebirdst)
<!-- badges: end -->

## Overview

The [eBird Status and Trends
project](https://science.ebird.org/en/status-and-trends) at the [Cornell
Lab of Ornithology](https://www.birds.cornell.edu/home) uses
machine-learning models to estimate distributions, relative abundances,
and population trends at high spatial and temporal resolution across the
full annual cycle of 1,119 bird species globally. These models learn the
relationships between bird observations collected through
[eBird](https://ebird.org/home) and a suite of remotely sensed habitat
variables, while accounting for the noise and bias inherent in community
science datasets, including variation in observer behavior and effort.
Interactive maps and visualizations of these model estimates can be
explored [online](https://science.ebird.org/en/status-and-trends), and
the [Status and Trends Data
Products](https://science.ebird.org/en/status-and-trends/download-data)
provide access to the data behind these maps and visualizations. The
`ebirdst` R package provides a set of tools for downloading these data
products, loading them into R, and using them for visualization and
analysis.

## Installation

Install `ebirdst` from GitHub with:

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("ebird/ebirdst")
```

This version of `ebirdst` is designed to work with the eBird Status and
Trends Data Products estimated for the year 2022 released in 2023.
**Users are strongly discouraged from comparing Status and Trends
results between years due to methodological differences between
versions.** If you have accessed and used previous versions and/or may
need access to previous versions for reasons related to reproducibility,
please contact <ebird@cornell.edu> and your request will be considered.

## Data access

Data access is granted through an Access Request Form at:
<https://ebird.org/st/request>. Access with this form generates a key to
be used with this R package and is provided immediately (as long as
commercial use is not requested). Our terms of use have been designed to
be quite permissive in many cases, particularly academic and research
use. When requesting data access, please be sure to carefully read the
terms of use and ensure that your intended use is not restricted.

After completing the Access Request Form, you will be provided a Status
and Trends Data Products access key, which you will need when
downloading data. To store the key so the package can access it when
downloading data, use the function `set_ebirdst_access_key("XXXXX")`,
where `"XXXXX"` is the access key provided to you.

**For those interested in accessing these data outside of R**, the most
widely used data products are available for direct download through the
[Status and Trends
website](https://science.ebird.org/en/status-and-trends). Spatial data
are accessible in widely adopted GeoTIFF and GeoPackage formats, which
can be opened in QGIS, ArcGIS, or other GIS software.

## Versions

The eBird Status and Trends Data Products provide estimates of relative
abundance, and other variables, for a particular year. This estimation
year is used to identify the version of the data products. Each version
of this R package is associated with a particular version of the data.
For example, the current version of the R package is 3.2022.3 and, as
indicated by the year in the version number, it is designed to work with
the 2022 data products. Every year, typically in November, the Status
and Trends Data Products are updated, and users are encouraged to update
this R package and transition to using the new version of the data
products. After the data products are updated, there will be a brief
period where access to the previous version is also provided, allowing
users to finish any analyses with this previous version. If you intended
to continue using the older data products during this period you must
not update the R package.

## Citation

If you use the the eBird Status and Trends data please cite it with:

<blockquote>
Fink, D., T. Auer, A. Johnston, M. Strimas-Mackey, S. Ligocki, O.
Robinson, W. Hochachka, L. Jaromczyk, C. Crowley, K. Dunham, A.
Stillman, I. Davies, A. Rodewald, V. Ruiz-Gutierrez, C. Wood. 2023.
eBird Status and Trends, Data Version: 2022; Released: 2023. Cornell Lab
of Ornithology, Ithaca, New York.
<a href="https://doi.org/10.2173/ebirdst.2022" class="uri">https://doi.org/10.2173/ebirdst.2022</a>
</blockquote>

[Download
BibTeX](https://raw.githubusercontent.com/ebird/ebirdst/main/ebirdst-citation.bib).

## Vignettes

For full package documentation, including a series of vignettes covering
the full spectrum from introductory to advanced usage, please see the
package [website](https://ebird.github.io/ebirdst/). The available
vignettes are:

- [Introduction to eBird Status Data
  Products](https://ebird.github.io/ebirdst/articles/status.html):
  covers data access, available data products, and structure and format
  of data files.
- [eBird Status Data Products
  Applications](https://ebird.github.io/ebirdst/articles/applications.html):
  demonstrates how to work with the raster data products and use them
  for a variety of common applications.
- [eBird Trends Data
  Products](https://ebird.github.io/ebirdst/articles/trends.html):
  covers downloading and working with the eBird Trends Data Products.

## Quick Start

This quick start guide shows how to download data and plot relative
abundance values similar to how they are plotted for the [eBird Status
and Trends weekly abundance
animations](https://science.ebird.org/en/status-and-trends/species/yebsap/abundance-map-weekly).
In this guide, and throughout all package documentation, a simplified
example dataset is used consisting of Yellow-bellied Sapsucker in
Michigan. For a full list of the species available for download, look at
the data frame `ebirst_runs`, which is included in this package.

**IMPORTANT: eBird Status and Trends Data Products are designed to be
downloaded and accessed using this R package. Downloaded data have a
specific file structure and changing file names or locations will
disrupt the ability of functions in this package to access the data. If
you prefer to access data for use outside of R, consider downloading
data via the [eBird Status and Trends
website](https://science.ebird.org/en/status-and-trends/download-data).**

``` r
library(fields)
library(rnaturalearth)
library(sf)
library(terra)
library(ebirdst)

# download example data, yellow-bellied sapsucker in michigan
ebirdst_download_status(species = "yebsap-example")

# load relative abundance raster stack with 52 layers, one for each week
abd <- load_raster("yebsap-example", resolution = "27km")

# load species specific mapping parameters
pars <- load_fac_map_parameters("yebsap-example")
# custom coordinate reference system
crs <- st_crs(pars$custom_projection)
# legend breaks
breaks <- pars$weekly_bins
# legend labels for top, middle, and bottom
labels <- pars$weekly_labels

# the date that each raster layer corresponds to is stored within the labels
weeks <- as.Date(names(abd))
print(weeks)
#>  [1] "2022-01-04" "2022-01-11" "2022-01-18" "2022-01-25" "2022-02-01" "2022-02-08" "2022-02-15" "2022-02-22"
#>  [9] "2022-03-01" "2022-03-08" "2022-03-15" "2022-03-22" "2022-03-29" "2022-04-05" "2022-04-12" "2022-04-19"
#> [17] "2022-04-26" "2022-05-03" "2022-05-10" "2022-05-17" "2022-05-24" "2022-05-31" "2022-06-07" "2022-06-14"
#> [25] "2022-06-21" "2022-06-28" "2022-07-05" "2022-07-12" "2022-07-19" "2022-07-26" "2022-08-02" "2022-08-09"
#> [33] "2022-08-16" "2022-08-23" "2022-08-30" "2022-09-06" "2022-09-13" "2022-09-20" "2022-09-27" "2022-10-04"
#> [41] "2022-10-11" "2022-10-18" "2022-10-25" "2022-11-01" "2022-11-08" "2022-11-15" "2022-11-22" "2022-11-29"
#> [49] "2022-12-06" "2022-12-13" "2022-12-20" "2022-12-27"

# select a week in the middle of the year
abd <- abd[[26]]

# project to species specific coordinates
# the nearest neighbor method preserves cell values across projections
abd_prj <- project(trim(abd), crs$wkt, method = "near")

# get reference data from the rnaturalearth package
# the example data currently shows only the US state of Michigan
wh_states <- ne_states(country = c("United States of America", "Canada"),
                       returnclass = "sf") %>% 
  st_transform(crs = crs) %>% 
  st_geometry()

# start plotting
par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))

# use raster bounding box to set the spatial extent for the plot
bb <- st_as_sfc(st_bbox(trim(abd_prj)))
plot(bb, col = "white", border = "white")
# add background reference data
plot(wh_states, col = "#cfcfcf", border = NA, add = TRUE)

# plot zeroes as light gray
plot(abd_prj, col = "#e6e6e6", maxpixels = ncell(abd_prj),
     axes = FALSE, legend = FALSE, add = TRUE)

# define color palette
pal <- ebirdst_palettes(length(breaks) - 1, type = "weekly")
# plot abundance
plot(abd_prj, col = pal, breaks = breaks, maxpixels = ncell(abd_prj),
     axes = FALSE, legend = FALSE, add = TRUE)

# state boundaries
plot(wh_states, add = TRUE, col = NA, border = "white", lwd = 1.5)

# legend
label_breaks <- seq(0, 1, length.out = length(breaks))
image.plot(zlim = c(0, 1), breaks = label_breaks, col = pal,
           smallplot = c(0.90, 0.93, 0.15, 0.85),
           legend.only = TRUE,
           axis.args = list(at = c(0, 0.5, 1), 
                            labels = round(labels, 2),
                            cex.axis = 0.9, lwd.ticks = 0))
```

<img src="man/figures/README-quick_start-1.png" width="100%" />
