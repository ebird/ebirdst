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
ebirdst_download_status("Golden Eagle", pattern = "3km")
# download 3 km estimates for Tui, a resident
ebirdst_download_status("Tui", pattern = "3km")

# ├ Loading data ----

# load weekly relative abundance estimates
abd_weekly <- load_raster("Golden Eagle")
# midpoint of weeks corresponding to each layer
weeks <- names(abd_weekly)
# subset to only weeks within may
abd_weekly_may <- abd_weekly[[weeks >= "2022-05-01" & weeks <= "2022-05-31"]]
# average across the weeks in may
abd_may <- mean(abd_weekly_may, na.rm = TRUE)

# load seasonal relative abundance estimates
abd_seasonal <- load_raster("Golden Eagle", period = "seasonal")
# subset to just the breeding season
abd_breeding <- abd_seasonal[["breeding"]]

# load full-year maximum relative abundance
abd_max <- load_raster("Golden Eagle", period = "full-year", metric = "max")

# load seasonal relative abundance estimates for a resident
load_raster("Tui", period = "seasonal")


# ├ Application 1: regional proportion of population ----

# Goal: find the % of the global Golden Eagle population in Wyoming seasonally

# load seasonal relative abundance
abd_seasonal <- load_raster("goleag", period = "seasonal")
# Wyoming boundary polygon from Natural Earth
wy_boundary <- ne_states(iso_a2 = "US") |>
  filter(name == "Wyoming") |>
  st_transform(st_crs(abd_seasonal)) |>
  vect()
# sum of abundance across all cells in Wyoming
wy_abd <- extract(abd_seasonal, wy_boundary,
                  fun = sum, na.rm = TRUE,
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
  st_transform(st_crs(abd_seasonal)) |>
  vect()
# mask to remove data outside North America
abd_seasonal_na <- mask(abd_seasonal, na_boundary)
# total North American abundance
na_abd <- global(abd_seasonal_na, fun = sum, na.rm = TRUE)
# proportion of global population
prop_na_pop <- as.numeric(wy_abd) / na_abd$sum
names(prop_na_pop) <- names(wy_abd)


# ├ Application 2: migration chronology ----

# Goal: plot the change in relative abundance throughout the year in for Hooded
# Warbler in Guatemala.

# download and load weekly relative abundance at 3 km
ebirdst_download_status("Hooded Warbler", pattern = "abundance_median_3km")
abd_weekly <- load_raster("Hooded Warbler")
# country boundary for Guatemala
gt_boundary <- ne_countries(country = "Guatemala", scale = 50) |>
  st_transform(st_crs(abd_weekly)) |>
  vect()
# mean weekly abundance within Guatemala
abd_gt <- extract(abd_weekly, gt_boundary,
                  fun = mean, na.rm = TRUE,
                  ID = FALSE)
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
ebirdst_download_status("Allen's Hummingbird", dry_run = TRUE)

# Exercise 2: regional proportion of population
# Pick a region and a species of interest to you, then estimate the proportion
# of the population within the region for each season. If it's a broadly
# distributed species, also estimate the proportion of the continental
# population.


# Part II: eBird Status Data Products Applications ----

# explore all available eBird Status data products
ebirdst_download_status("Golden Eagle", download_all = TRUE, dry_run = TRUE)
# download everything!
ebirdst_download_status("Golden Eagle", download_all = TRUE)
# load the weekly proportion of population at 27 km, then show that it sums to 1
prop_pop <- load_raster("Golden Eagle",
                        product = "proportion-population",
                        resolution = "27km")
global(prop_pop, fun = sum, na.rm = TRUE)
# load the full-year maximum relative abundance
abd_fy_max <- load_raster("Golden Eagle",
                          product = "abundance",
                          period = "full-year",
                          metric = "max")
# load the regional stats
regional_stats <- load_regional_stats("Golden Eagle")
# load smoothed range polygon at 27 km resolution
ranges <- load_ranges("Golden Eagle", resolution = "27km", smoothed = TRUE)

# ├ Application 1: multi-species migration chronology ----

# Goal: plot migration chronologies with uncertainty estimates for the set of
# six grassland species below, comparing changes in the proportion of the
# population in Montana throughout the year.

# species list
grassland_species <- c("Baird's Sparrow",
                       "Bobolink",
                       "Chestnut-collared Longspur",
                       "Sprague's Pipit",
                       "Upland Sandpiper",
                       "Western Meadowlark")

# Montana boundary polygon from Natural Earth
# note: you could use any region here, e.g. a shapefile, read in using read_sf
mt_boundary <- ne_states(iso_a2 = "US") |>
  filter(name == "Montana") |>
  # transform coordinate system to match the raster data
  st_transform("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs") |>
  vect()

# loop over each species, generate a migration chronology for each
chronology <- NULL
for (species in grassland_species) {
  # download the data products for this species
  # only weekly 27km relative abundance, median and confidence limits
  ebirdst_download_status(species,
                          pattern = "abundance_(median|lower|upper)_27km")

  # load the median weekly relative abundance and lower/upper confidence limits
  abd_median <- load_raster(species, resolution = "27km")
  abd_lower <- load_raster(species, metric = "lower", resolution = "27km")
  abd_upper <- load_raster(species, metric = "upper", resolution = "27km")

  # convert from relative abundance to proportion of population
  abd_total <- global(abd_median, fun = sum, na.rm = TRUE)$sum
  pop_median <- abd_median / abd_total
  pop_lower <- abd_lower / abd_total
  pop_upper <- abd_upper / abd_total

  # estimate the proportion of population in Montana with confidence limits
  # note: you could also mask and sum here
  pop_median_region <- extract(pop_median, mt_boundary,
                               fun = sum, na.rm = TRUE, ID = FALSE)
  pop_lower_region <- extract(pop_lower, mt_boundary,
                              fun = sum, na.rm = TRUE, ID = FALSE)
  pop_upper_region <- extract(pop_upper, mt_boundary,
                              fun = sum, na.rm = TRUE, ID = FALSE)

  # convert to a data frame in long format (one row per week)
  chronology <- data.frame(species = species,
                           week = as.Date(names(pop_median_region)),
                           median = as.numeric(pop_median_region),
                           lower = as.numeric(pop_lower_region),
                           upper = as.numeric(pop_upper_region)) |>
    bind_rows(chronology)
}

# plot the migration chronologies for each species
ggplot(chronology) +
  aes(x = week, y = median, color = species, fill = species) +
  geom_ribbon(aes(ymin = lower, ymax = upper), color = NA, alpha = 0.2) +
  geom_line() +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  scale_y_continuous(labels = scales::label_percent()) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Week",
       y = "% of population in Montana",
       title = "Migration chronologies for grassland birds in Montana",
       color = NULL, fill = NULL) +
  theme(legend.position = "bottom")

# ├ Application 2: areas of importance ----

# Goal: map areas of importance during the breeding season for the set of six
# grassland species in Montana.

# start by producing a map of richness for these species
range_mt <- list()
for (species in grassland_species) {
  # download seasonal abundance at 3km
  ebirdst_download_status(species, pattern = "abundance_seasonal_mean_3km")

  # load breeding season relative abundance
  abd <- load_raster(species, period = "seasonal") |>
    subset("breeding")
  # crop and mask to Montana
  abd_masked <- mask(crop(abd, mt_boundary), mt_boundary)
  # convert to binary, presence-absence
  range_mt[[species]] <- abd_masked > 0
}
# sum across species to calculate richness
richness <- sum(rast(range_mt), na.rm = TRUE)
# make a simple map
plot(richness, axes = FALSE)

# for more granularity, combine proportion of population across species
prop_pop_mt <- list()
for (species in grassland_species) {
  # download seasonal abundance at 3km
  ebirdst_download_status(species,
                          pattern = "proportion-population_seasonal_mean_3km")

  # load breeding season proportion of population
  prop_pop <- load_raster(species,
                          product = "proportion-population",
                          period = "seasonal") |>
    subset("breeding")
  # crop and mask to Montana
  prop_pop_mt[[species]] <- mask(crop(prop_pop, mt_boundary), mt_boundary)
}
# take mean across species
importance <- mean(rast(prop_pop_mt), na.rm = TRUE)
# drop zeros
importance <- ifel(importance == 0, NA, importance)
# drop anything below the median
cutoff <- global(importance, quantile, probs = 0.5, na.rm = TRUE) |>
  as.numeric()
importance <- ifel(importance > cutoff, importance, NA)
# make a simple map
plot(importance, axes = FALSE)
plot(mt_boundary, col = "grey", axes = FALSE, add = TRUE)
plot(importance, axes = FALSE, legend = FALSE, add = TRUE)

# make a slightly nicer map
# reproject
importance_proj <- trim(project(importance, "ESRI:102003"))
mt_boundary_proj <- project(mt_boundary, "ESRI:102003")
# basemap
par(mar = c(0, 0, 0, 0))
plot(mt_boundary_proj, col = "grey", axes = FALSE,
     main = "Areas of importance for grassland birds in Montana")
# add importance raster
plot(importance_proj, legend = FALSE, add = TRUE)
# add legend
fields::image.plot(zlim = c(0, 1), legend.only = TRUE,
                   col = viridis::viridis(100),
                   breaks = seq(0, 1, length.out = 101),
                   smallplot = c(0.15, 0.85, 0.12, 0.15),
                   horizontal = TRUE,
                   axis.args = list(at = c(0, 0.5, 1),
                                    labels = c("Low", "Medium", "High"),
                                    fg = "black", col.axis = "black",
                                    cex.axis = 0.75, lwd.ticks = 0.5,
                                    padj = -1.5),
                   legend.args = list(text = "Relative Importance",
                                      side = 3, col = "black",
                                      cex = 1, line = 0))

# ├ Exercises ----

# Exercises 1: repeat the migration chronology application demonstrated in the
# webinar, but try plotting the proportion of the population within the
# contiguous United States (i.e. all US states except Alaska and Hawaii) rather
# than the proportion of global population. Hint: this will require cropping and
# masking the relative abundance rasters to a boundary of the unites states,
# which is provided below.
us_boundary <- ne_states(iso_a2 = "US") |>
  filter(!name %in% c("Alaska", "Hawaii")) |>
  st_union() |>
  st_transform("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs") |>
  vect()

# Exercise 2: select a group of 5-10 species and a region of interest to you.
# Generate a map identifying areas of importance for these species by finding
# the mean proportion of population across the species in your region of
# interest for either the breeding or non-breeding season. Experiment with
# different quantile cutoffs (e.g. median, 70th quantile, 90th quantile) to see
# how that impacts the areas of importance identified.
