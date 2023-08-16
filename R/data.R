#' Data frame of species with eBird Status Data Products
#'
#' A dataset containing the species for which eBird Status Data Products are
#' available In addition, the dates defining the boundaries of the seasons are
#' provided. These seasons are defined on a species-specific basis through
#' expert review. For information on the details of defining seasons, please see
#' the [seasons section of the
#' FAQ](https://ebird.org/science/status-and-trends/faq#seasons). Note that
#' missing dates imply that a season failed expert review for that species
#' within that season.
#'
#' @format A data frame with 15 variables:
#' \describe{
#'   \item{species_code}{Six letter eBird code in eBird Taxonomy v2018}
#'   \item{scientific_name}{Scientific name from eBird Taxonomy v2018}
#'   \item{common_name}{English common name from eBird Taxonomy v2018}
#'   \item{resident}{Classifies this species a resident or a migrant}
#'   \item{breeding_quality}{Breeding season quality}
#'   \item{breeding_range_modeled}{Is the full range modeled?}
#'   \item{breeding_start}{Breeding season start date}
#'   \item{breeding_end}{Breeding season start date}
#'   \item{nonbreeding_quality}{Non-breeding season quality}
#'   \item{nonbreeding_range_modeled}{Is the full range modeled?}
#'   \item{nonbreeding_start}{Non-breeding season start date}
#'   \item{nonbreeding_end}{Non-breeding season start date}
#'   \item{postbreeding_migration_quality}{Post-breeding season quality}
#'   \item{postbreeding_migration_range_modeled}{Is the full range modeled?}
#'   \item{postbreeding_migration_start}{Post-breeding season start date}
#'   \item{postbreeding_migration_end}{Post-breeding season start date}
#'   \item{prebreeding_migration_quality}{Pre-breeding season quality}
#'   \item{prebreeding_migration_range_modeled}{Is the full range modeled?}
#'   \item{prebreeding_migration_start}{Pre-breeding season start date}
#'   \item{prebreeding_migration_end}{Pre-breeding season start date}
#'   \item{resident_quality}{Resident quality}
#'   \item{resident_start}{For resident species, the year-round start date}
#'   \item{resident_end}{For resident species, the year-round end date}
#' }
"ebirdst_runs"

#' Data frame of species with eBird Trends Data Products
#'
#' This data frame contains a list of species for which eBird Trends Data
#' Products are available including model-level information for the trend model.
#' The two predictive performance metrics (`rsquared` and `beta0`) are based
#' on a comparison of actual and estimated percent per year trends for a suite
#' of simulations (see Fink et al. 2023 for further details).
#'
#' @format A data frame with 11 variables:
#' \describe{
#'   \item{run_name}{Combination of species code and run name.}
#'   \item{species_code}{The alphanumeric eBird species code uniquely
#'         identifying the species.}
#'   \item{common_name}{The English common name of the species.}
#'   \item{season}{Season that the trend was estimated for: breeding,
#'         nonbreeding, or resident.}
#'   \item{region}{The geographic region that the trend model was run for. Note
#'         that broadly distributed species (e.g. Barn Swallow) will only have
#'         trend estimates for a regional subset of their full range.}
#'   \item{start_date}{The start date (`MM-DD` format) of the season for which
#'         the trend was estimated.}
#'   \item{end_date}{The end date (`MM-DD` format) of the season for which the
#'         trend was estimated.}
#'   \item{start_year}{The start date of the trend time period.}
#'   \item{end_year}{The end date of the trend time period.}
#'   \item{rsquared}{R-squared value comparing the actual and estimated trends
#'         from the simulations.}
#'   \item{beta0}{The slope of a linear model fitting actual vs. estimated
#'         trends (`actual ~ estimated`) for the simulations. Positive values
#'         of `beta0` indicate that the models are systematically
#'         *underestimating* the simulated trend for this species.}
#' }
"ebirdst_trends_runs"

#' eBird Status and Trends predictors
#'
#' A data frame of the predictors used in the eBird Status and Trends models.
#' These include effort variables (e.g. distance traveled, number of observers,
#' etc.) in addition to land and water cover variables. These landcover
#' variables are derived from the MODIS MCD12Q1 500 m landcover product, and for
#' each land cover class two FRAGSTATS metrics are calculated within a 1.5 km
#' buffer around each checklist: % landcover (PLAND) and edge density (ED).
#'
#' @format A data frame with 74 rows and 5 columns:
#' \describe{
#'   \item{predictor}{Predictor variable name.}
#'   \item{predictor_label}{Descriptive labels for predictors for plotting and
#'         translating the cryptic variables names (e.g. `umd_fs_c1` is
#'         Evergreen Needleleaf Forest.}
#'   \item{lc_class}{For the land and water cover FRAGSTATS variables, this
#'         gives the associated landcover class. It can be used for grouping
#'         and summarizing the four FRAGSTATS metrics to the level of the
#'         landcover class.}
#'   \item{lc_class_label}{Similar to `predictor_label`; however, this variable
#'         gives the FRAGSTATS metrics a single name for the landcover
#'         class.}
#' }
"ebirdst_predictors"

#' eBird Status and Trends weeks
#'
#' eBird Status and Trends predictions are made for each of 52 weeks of the
#' year. This data frame provides the boundaries of the weeks.
#'
#' @format A data frame with 52 rows and 5 columns:
#' \describe{
#'   \item{week_number}{Integer week number from 1-52.}
#'   \item{date}{Date of the midpoint of the week.}
#'   \item{week_midpoint}{Date of the midpoint of the week expressed as a
#'         fraction of the year, i.e. a number from 0-1.}
#'   \item{week_start}{Date of the start of the week expressed as a fraction of
#'         the year, i.e. a number from 0-1.}
#'   \item{week_end}{Date of the end of the week expressed as a fraction of the
#'         year, i.e. a number from 0-1.}
#' }
"ebirdst_weeks"
