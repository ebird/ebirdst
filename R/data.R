#' Data frame of species with eBird Status and Trends Data Products
#'
#' A dataset listing the species for which eBird Status and Trends Data Products
#' are available, with additional information relevant to both the Status and
#' Trends results for each species.
#'
#' For the Status Data Products, the dates defining the boundaries of the
#' seasons are provided in additional to a quality rating from 0-3 for each
#' season. These dates and quality ratings are assigned through a process of
#' [expert review](https://science.ebird.org/status-and-trends/faq#seasons).
#' expert review. Note that missing dates imply that a season failed expert
#' review for that species within that season.
#'
#' Trends Data Products are only available for a subset of species, indicated by
#' the `has_trends` variable, and for each species the trends is estimated for a
#' single season. The two predictive performance metrics (`rsquared` and
#' `beta0`) are based on a comparison of actual and estimated percent per year
#' trends for a suite of simulations (see Fink et al. 2023 for further details).
#' The trends regions are defined as follows:
#'
#' - `aus_nz`: Australia and New Zealand
#' - `iberia`: Spain and Portugal
#' - `india_se_asia`: India, Nepal, Bhutan, Sri Lanka, Thailand, Cambodia,
#' Malaysia, Brunei, Singapore, and Philippines
#' - `japan`: Japan
#' - `north_america`: North America including Mexico, Central America, and the
#' Caribbean, but excluding Nunavut, North West Territories, and Hawaii
#' - `south_africa`: South Africa, Lesotho, and Eswatini
#' - `south_america`: Colombia, Ecuador, Peru, Chile, Argentina, and Uruguay
#' - `taiwan`: Taiwan
#' - `turkey_plus`: Turkey, Cyprus, Israel, Palestine, Greece, Armenia, and
#' Georgia
#'
#' @format A data frame with 29 variables:
#' - `species_code`: alphanumeric eBird species code uniquely identifying the
#' species
#' - `scientific_name`: scientific name.
#' - `common_name`: English common name.
#' - `is_resident`: classifies this species a resident or a migrant.
#' - `breeding_quality`: breeding season quality.
#' - `breeding_start`: breeding season start date.
#' - `breeding_end`: breeding season start date.
#' - `nonbreeding_quality`: non-breeding season quality.
#' - `nonbreeding_start`: non-breeding season start date.
#' - `nonbreeding_end`: non-breeding season start date.
#' - `postbreeding_migration_quality`: post-breeding season quality.
#' - `postbreeding_migration_start`: post-breeding season start date.
#' - `postbreeding_migration_end`: post-breeding season start date.
#' - `prebreeding_migration_quality`: pre-breeding season quality.
#' - `prebreeding_migration_start`: pre-breeding season start date.
#' - `prebreeding_migration_end`: pre-breeding season start date.
#' - `resident_quality`: resident quality.
#' - `resident_start`: for resident species, the year-round start date.
#' - `resident_end`: for resident species, the year-round end date.
#' - `status_version_year`: the release version of the Status data products.
#' - `has_trends`: whether or not this species has trends estimates.
#' - `trends_season`: season that the trend was estimated for: breeding,
#' nonbreeding, or resident.
#' - `trends_region`: the geographic region that the trend model was run for.
#' Note that broadly distributed species (e.g. Barn Swallow) will only have
#' trend estimates for a regional subset of their full range.
#' - `trends_start_year`: start year of the trend time period.
#' - `trends_end_year`: end year of the trend time period.
#' - `trends_start_date`: start date (`MM-DD` format) of the season for which
#' the trend was estimated.
#' - `trends_end_date`: end date (`MM-DD` format) of the season for which the
#' trend was estimated.
#' - `rsquared`: R-squared value comparing the actual and estimated trends from
#' the simulations.
#' - `beta0`: the intercept of a linear model fitting actual vs. estimated
#' trends (`actual ~ estimated`) for the simulations. Positive values of `beta0`
#' indicate that the models are systematically *underestimating* the simulated
#' trend for this species.
#' - `trends_version_year`: the release version of the Trends data products.
"ebirdst_runs"

#' eBird Status and Trends predictor variables
#'
#' A data frame of the predictors used in the eBird Status and Trends models.
#' These include effort variables (e.g. distance traveled, number of observers,
#' etc.) in addition to variables describing the environment (e.g. elevation,
#' land cover, water cover, etc.). The environmental variables are derived by
#' summarizing remotely sensed datasets (described in
#' [ebirdst_predictor_descriptions]) over a 3 km diameter neighborhood around
#' each checklist. For categorical datasets, two variables are generated for
#' each class describing the percent cover (`pland`) and edge density (`ed`).
#'
#' @format A data frame with 150 rows and 4 columns:
#' - `predictor`: predictor name.
#' - `dataset`: dataset name, which can be cross referenced in
#' [ebirdst_predictor_descriptions] for further details.
#' - `class`: class number or name for categorical variables.
#' - `label`: descriptive labels for each predictor variable.
"ebirdst_predictors"

#' eBird Status and Trends predictors descriptions
#'
#' Details on the eBird Status and Trends predictor variables or, for variables
#' all derived from the same dataset, details on the dataset.
#'
#' @format A data frame with 37 rows and 4 columns
#' - `dataset`: dataset name.
#' - `predictor`: predictor name or, if multiple variables are derived from
#' this dataset, the pattern used to generate the names.
#' - `description`: detailed description of the dataset or variable.
#' - `reference`: a reference to consult for further information on the dataset.
"ebirdst_predictor_descriptions"
