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
abd_seasonal[["breeding"]]

# load full-year maximum relative abundance


# load seasonal relative abundance estimates for a resident



# ├ Application 1: regional proportion of population ----

# Goal: find the % of the global Golden Eagle population in Wyoming seasonally

# load seasonal relative abundance
abd_seasonal <- load_raster("Golden Eagle", period = "seasonal")
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
na_abd <- as.numeric(wy_abd) / na_abd$sum
names(na_abd) <- names(wy_abd)


# ├ Application 2: migration chronology ----

# Goal: plot the change in relative abundance throughout the year in for Hooded
# Warbler in Guatemala.

# download and load weekly relative abundance at 3 km
ebirdst_download_status("hoowar", pattern = "abundance_median_3km")
abd_weekly <- load_raster("hoowar")
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

# download everything!

# load the weekly proportion of population at 27 km, then show that it sums to 1

# load the full-year maximum relative abundance

# load the regional stats

# load smoothed range polygon at 27 km resolution


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
# produce binary range rasters for each species in Montana
range_mt <- list()
for (species in grassland_species) {

}
# sum across species to calculate richness

# make a simple map


# for more granularity, use mean proportion of population
# produce proportion of population layers for each species masked to Montana
prop_pop_mt <- list()
for (species in grassland_species) {

}
# take mean across species

# drop zeros

# drop anything below the median


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
# webinar, but try plotting the proportion of the population relative to the
# contiguous United States (i.e. all US states except Alaska and Hawaii) rather
# than the proportion of the global population. Hint: this will require cropping
# and masking the relative abundance rasters to a boundary of the Unites States,
# which is provided below.

# boundary of the contiguous united states
us_boundary <- ne_states(iso_a2 = "US") |>
  filter(!name %in% c("Alaska", "Hawaii")) |>
  st_union() |>
  st_transform("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs") |>
  vect()

# Exercise 2: select a group of 5-10 species and a region of interest to you.
# Generate a map identifying areas of importance for these species by finding
# the mean proportion of population across the species in your region of
# interest for either the breeding or non-breeding season. Experiment with
# different quantile cutoffs (e.g. median, 70th quantile, 80th quantile) to see
# how that impacts the areas of importance identified.


# Part III: eBird Trends Data Products Applications ----

# ├ Exercise 1 Solution ----

# Goal: repeat the multi-species migration chronology application, this time
# showing the proportion of the contiguous US population rather than the
# proportion of the entire modeled population

# species list
grassland_species <- c("Baird's Sparrow",
                       "Bobolink",
                       "Chestnut-collared Longspur",
                       "Sprague's Pipit",
                       "Upland Sandpiper",
                       "Western Meadowlark")

# Montana boundary polygon from Natural Earth
mt_boundary <- ne_states(iso_a2 = "US") |>
  filter(name == "Montana") |>
  # transform coordinate system to match the raster data
  st_transform("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs") |>
  vect()

# boundary of the contiguous United States from Natural Earth
us_boundary <- ne_states(iso_a2 = "US") |>
  filter(!name %in% c("Alaska", "Hawaii")) |>
  st_union() |>
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
  abd_median <- load_raster(species, resolution = "27km") |>
    # masking to US boundary will ensure we have proportions of US population
    mask(us_boundary)
  abd_lower <- load_raster(species, metric = "lower", resolution = "27km") |>
    # masking to US boundary will ensure we have proportions of US population
    mask(us_boundary)
  abd_upper <- load_raster(species, metric = "upper", resolution = "27km") |>
    # masking to US boundary will ensure we have proportions of US population
    mask(us_boundary)

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
# note: there are now some missing values in chronology because some species
# have no US population at all during some weeks (e.g. Upload Sandpiper
# migrates entirely to southern South American)
ggplot(chronology[complete.cases(chronology), ]) +
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

# ├ Exercise 2 Solution ----

# Goal: repeat the areas of importance exercise for a set of 6 woodpeckers in
# Washington this time using a 70th percentile cutoff rather than a median
# cutoff.

# species list, note some are migrants and others are residents
woodpecker_species <- c("Black-backed Woodpecker",
                        "Downy Woodpecker",
                        "Hairy Woodpecker",
                        "Lewis's Woodpecker",
                        "Northern Flicker",
                        "White-headed Woodpecker")

# Washington boundary polygon from Natural Earth
wa_boundary <- ne_states(iso_a2 = "US") |>
  filter(name == "Washington") |>
  # transform coordinate system to match the raster data
  st_transform("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs") |>
  vect()

# produce proportion of population layers for each species masked to Washington
prop_pop_wa <- list()
for (species in woodpecker_species) {
  # download seasonal abundance at 3km
  ebirdst_download_status(species,
                          pattern = "proportion-population_seasonal_mean_3km")

  # check if this species is a resident
  is_resident <- filter(ebirdst_runs, common_name == species)[["is_resident"]]

  # load breeding season proportion of population for migrants
  # and resident season proportion of population for residents
  if (is_resident) {
    season_name <- "resident"
  } else {
    season_name <- "breeding"
  }
  prop_pop <- load_raster(species,
                          product = "proportion-population",
                          period = "seasonal") |>
    subset(season_name)
  # crop and mask to Washington
  prop_pop_wa[[species]] <- mask(crop(prop_pop, wa_boundary), wa_boundary)
}
# take mean across species
importance <- mean(rast(prop_pop_wa), na.rm = TRUE)
# drop zeros
importance <- ifel(importance == 0, NA, importance)
# drop anything below the 70th percentile
cutoff <- global(importance, quantile, probs = 0.7, na.rm = TRUE) |>
  as.numeric()
importance <- ifel(importance > cutoff, importance, NA)

# make a slightly nicer map
# reproject
importance_proj <- trim(project(importance, "ESRI:102003"))
wa_boundary_proj <- project(wa_boundary, "ESRI:102003")
# basemap
par(mar = c(0, 0, 0, 0))
plot(wa_boundary_proj, col = "grey", axes = FALSE,
     main = "Areas of importance for woodpeckers in Washington")
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

# ├ Introduction to the eBird Trends Data Products ----

# ebirdst runs that have trends results

# download trends for Sage Thrasher

# load trends into R session


# convert percent per year trends from data frame format to raster

# convert cumulative trends from data frame format to raster

# save to geotiff for use in QGIS or ArcGIS


# convert to spatial points for use with sf package

# save to GeoPackage for use in QGIS or ArcGIS, could also save to shapefile


# conversion from annual to cumulative trend
ppy_trend <- 2
start_year <- 2012
end_year <- 2022
cumulative_trend <- 100 * ((1 + ppy_trend / 100)^(end_year - start_year) - 1)

# ├ Application 1: regional trends ----

# Goal: estimate the % per year trend for Sage Thrasher in the Great Basin Bird
# Conservation Region (BCR 9). BCR polygons can be downloaded from
# https://www.birdscanada.org/bird-science/nabci-bird-conservation-regions

# load BCR polygons that you've downloaded and unzipped, then subset to BCR 9


# load trend estimates for Sage Thrasher

# convert to spatial sf format

# subset to just cells within BCR 19

# calculate abundance-weighted regional trend


# ├ Application 2: multi-region trends with uncertainty ----

# Goal: estimate the % per year trend with 80% confidence limits for Sage
# Thrasher for each state in the contiguous United States.

# polygon boundaries of each state in the contiguous US
states <- gadm(country = "USA", level = 1, path = tempdir()) |>
  st_as_sf() |>
  select(state = ISO_1) |>
  filter(!state %in% c("US-AK", "US-HI"))

# load fold-level trend estimates for Sage Thrasher

# convert to spatial sf format

# attach state to the fold-level trends data


# abundance-weighted average trend by region and fold


# summarize across folds for each state

# join trends to state polygons and make a map
trends_states_sf <- left_join(states, trends_states, by = "state")
ggplot(trends_states_sf) +
  geom_sf(aes(fill = abd_ppy_median)) +
  scale_fill_distiller(palette = "Reds",
                       limits = c(NA, 0),
                       na.value = "grey80") +
  guides(fill = guide_colorbar(title.position = "top", barwidth = 15)) +
  labs(title = "Sage Thrasher state-level breeding trends 2012-2022",
       fill = "Relative abundance trend [% change / year]") +
  theme_bw() +
  theme(legend.position = "bottom")

# ├ Application 3: multi-species trends ----

# Goal: estimate the mean trend for a set three species representing the
# sagebrush bird community.

# species list

# ensure that all species are for the same region and season


# download trends for these species and load them into R


# calculate cell-wise mean trend


# convert the points to sf format
# only consider cells where all three species occur


# make a map
ggplot(all_species) +
  geom_sf(aes(color = abd_ppy), size = 2) +
  scale_color_gradient2(low = "#CB181D", high = "#2171B5",
                        limits = c(-4, 4),
                        oob = scales::oob_squish) +
  guides(color = guide_colorbar(title.position = "left", barheight = 15)) +
  labs(title = "Sagebrush species breeding trends (2012-2022)",
       color = "Relative abundance trend [% change / year]") +
  theme_bw() +
  theme(legend.title = element_text(angle = 90))
