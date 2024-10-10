library(dplyr)
library(ebirdst)
library(ggplot2)
library(rnaturalearth)
library(sf)
library(terra)


# Part I: Introduction ----

# explore available status data products species
View(ebirdst_runs)
# examine seasonal dates and quality scores for chimney swift
glimpse(filter(ebirdst_runs, common_name == "Chimney Swift"))

# ├ Downloading data ----

# default data download location
ebirdst_data_dir()
# list available files for Golden Eagle
ebirdst_download_status("Golden Eagle", dry_run = TRUE)
# download 3 km estimates for Golden Eagle, a migrant
ebirdst_download_status("Golden Eagle", pattern = "27km")
# download 3 km estimates for Tui, a resident
ebirdst_download_status("Tui", pattern = "27km")

# ├ Loading data ----

# load weekly relative abundance estimates at 27 km resolution
abd_weekly <- load_raster("Golden Eagle", resolution = "27km")
# midpoint of weeks corresponding to each layer
weeks <- names(abd_weekly)
# subset to only weeks within may
abd_weekly_may <- abd_weekly[[weeks >= "2022-05-01" & weeks <= "2022-05-31"]]
# average across the weeks in may
abd_may <- mean(abd_weekly_may, na.rm = TRUE)

# load seasonal relative abundance estimates at 27 km resolution
abd_seasonal <- load_raster("Golden Eagle",
                            period = "seasonal",
                            resolution = "27km")
# subset to just the breeding season
abd_breeding <- abd_seasonal[["breeding"]]

# load seasonal relative abundance estimates for a resident
load_raster("Tui", period = "seasonal", resolution = "27km")


# ├ Application 1: regional proportion of population ----

# Goal: find the % of the global Golden Eagle population in Wyoming seasonally

# load seasonal relative abundance
abd_seasonal <- load_raster("goleag", period = "seasonal", resolution = "27km")
# Wyoming boundary polygon from Natural Earth
wy_boundary <- ne_states(iso_a2 = "US") |>
  filter(name == "Wyoming") |>
  vect()
# sum of abundance across all cells in Wyoming
wy_abd <- extract(abd_seasonal, wy_boundary,
                  fun = "sum", na.rm = TRUE,
                  weights = TRUE, ID = FALSE)
# total global abundance
total_abd <- global(abd_seasonal, fun = sum, na.rm = TRUE)
# proportion of global population
prop_global_pop <- as.numeric(wy_abd) / total_abd$sum
names(prop_global_pop) <- names(wy_abd)

# Goal: find the % of the North American population

# north american boundary
na_boundary <- ne_countries(scale = 50) |>
  filter(continent == "North America") |>
  st_union() |>
  vect()
# sum of abundance across all cells in North America
na_abd <- extract(abd_seasonal, na_boundary,
                  fun = "sum", na.rm = TRUE,
                  weights = TRUE, ID = FALSE)
# proportion of global population
prop_na_pop <- as.numeric(wy_abd) / as.numeric(na_abd)
names(prop_na_pop) <- names(wy_abd)


# ├ Application 2: migration chronology ----

# Goal: plot the change in relative abundance throughout the year in for Hooded
# Warbler in Guatemala.

# download and load weekly relative abundance at 27 km
ebirdst_download_status("Hooded Warbler", pattern = "abundance_median_27km")
abd_weekly <- load_raster("Hooded Warbler", resolution = "27km")
# country boundary for Guatemala
gt_boundary <- ne_countries(country = "Guatemala", scale = 50) |>
  vect()
# mean weekly abundance within Guatemala
abd_gt <- extract(abd_weekly, gt_boundary,
                  fun = "mean", na.rm = TRUE, ID = FALSE)
abd_gt <- data.frame(week = as.Date(names(abd_weekly)),
                     abd = as.numeric(abd_gt))
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
ebirdst_download("Allen's Hummingbird", dry_run = TRUE)

# Exercise 2: regional proportion of population
# Pick a region and a species of interest to you and estimate the proportion of
# the global population within the region for each season. If it's a broadly
# distributed species, also estimate the proportion of the continental
# population.
