library(dplyr)
library(ebirdst)
library(ggplot2)
library(rnaturalearth)
library(sf)
library(terra)


# Part I: Introduction ----

# explore available status data products species

# examine seasonal dates and quality scores for chimney swift


# ├ Downloading data ----

# default data download location

# list available files for Golden Eagle

# download 3 km estimates for Golden Eagle, a migrant

# download 3 km estimates for Tui, a resident


# ├ Loading data ----

# load weekly relative abundance estimates

# midpoint of weeks corresponding to each layer

# subset to only weeks within may

# average across the weeks in may


# load seasonal relative abundance estimates

# subset to just the breeding season


# load full-year maximum relative abundance


# load seasonal relative abundance estimates for a resident



# ├ Application 1: regional proportion of population ----

# Goal: find the % of the global Golden Eagle population in Wyoming seasonally

# load seasonal relative abundance

# Wyoming boundary polygon from Natural Earth
wy_boundary <- ne_states(iso_a2 = "US") |>
  filter(name == "Wyoming") |>
  st_transform(st_crs(abd_seasonal)) |>
  vect()
# sum of abundance across all cells in Wyoming

# total global abundance

# proportion of global population


# Goal: find the % of the North American population

# north american boundary
na_boundary <- ne_countries(scale = 50) |>
  filter(continent == "North America") |>
  st_union() |>
  st_transform(st_crs(abd_seasonal)) |>
  vect()
# mask to remove data outside North America

# total North American abundance

# proportion of global population



# ├ Application 2: migration chronology ----

# Goal: plot the change in relative abundance throughout the year in for Hooded
# Warbler in Guatemala.

# download and load weekly relative abundance at 3 km

# country boundary for Guatemala
gt_boundary <- ne_countries(country = "Guatemala", scale = 50) |>
  st_transform(st_crs(abd_weekly)) |>
  vect()
# mean weekly abundance within Guatemala


# plot migration chronology
ggplot(abd_gt) +
  aes(x = week, y = abd) +
  geom_line() +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  labs(x = "Week",
       y = "Mean relative abundance in Guatemala",
       title = "Migration chronology for Hooded Warbler in Guatemala")

# ├ Exercises ----

# Exercises 1: make sure your data access is working
# Install ebirdst, set your access key, then run the following code.
# Email Matt if this gives an error.
ebirdst_download_status("Allen's Hummingbird", dry_run = TRUE)

# Exercise 2: regional proportion of population
# Pick a region and a species of interest to you, then estimate the proportion
# of the population within the region for each season. If it's a broadly
# distributed species, also estimate the proportion of the continental
# population.
