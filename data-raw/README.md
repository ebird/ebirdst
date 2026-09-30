# eBird Status and Trends Release Checklist

This document outlines the steps required to update this R package for a new production release. These steps can only be run by members of the Status and Trends team that have access to the various private datasets.

## Define input dataset locations

Set the following environmental variables either permanently in your `~/.Renviron` file or temporarily using `Sys.setenv()`. The three `GS_KEY` variables are Google Sheets keys that can be extracted from the relevant URLs. `EBIRDST_S3_BUCKET` is the `s3://` URL for the S3 bucket containing the eBird Status and Trends data products. To complete the steps outlined in this document you will need read access to the S3 bucket.

```
EBIRDST_S3_BUCKET=
EBIRDST_STATUS_GS_KEY=
EBIRDST_TRENDS_GS_KEY=
EBIRDST_FEATURES_GS_KEY=
```

## Update version years in R function

Update `status_version_year` and `trends_version_year` in the function definition for `ebirdst_version()` in `R/download.R`. For Status, this is the prediction year and for Trends it's the final year of the trend time series. If only one of Status or Trends is being updated in a given year, only update the year for the products that have been updated.

Rebuild the R package with `devtools::install()` prior to proceeding since later steps rely on these version years.

## Generate the internal data frames

Run the following scripts to generate data frames that will be accessible when the package is loaded:

- `data-raw/generate_ebirdst_runs.R`: creates `ebirdst_runs`, which contains one row for each species with Status and/or Trends data products.
- `data-raw/generate_ebirdst_predictors.R`: creates `ebirdst_predictors` and `ebirdst_predictor_descriptions`, which provide information on the predictors used in the modeling process.

## Generate example data

We provide a small example dataset for Yellow-bellied Sapsucker in Michigan that can be accessed without an API key. The data are stored in the GitHub repository: https://github.com/ebird/ebirdst_example-data. Clone that repository locally, run the script `generate_example_data.R` in that repository, then commit and push to GitHub to produce the example dataset.

Then in the main `ebirdst` project, run `generate_example_file_lists.R` to create lists of example files for Status and Trends that are stored as internal package data.

## Update citation

Update `R/zzz.R` to reflect the updated citation for this release.

## Update the change log

Each release Tom Auer compiles a detailed change log for anything that has changed since the last release. This includes everything from updated covariates to base model changes to difference in released data products. Tom will typically share this change log as a Google Document, which should be proofread, converted to markdown, and added to the change log vignette in `vignettes/articles/product-changelog.qmd`.

## Check all vignettes

Run all the code in the vignettes to make sure everything works with the newly released data products. Update as needed.

## Standard R package/CRAN tasks

1. Update the version in `DESCRIPTION` to reflect the release year.
2. Update `NEWS.md` and `cran-comments.md`.
4. Run all the steps in `makefile.R`. 