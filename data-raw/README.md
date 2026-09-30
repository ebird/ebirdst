# eBird Status and Trends Release Checklist

This document outlines the steps required to update this R package for a new production release. These steps can only be run by members of the Status and Trends team that have access to the various private datasets.

## Define input dataset locations

Set the following environmental variables either permanently in your `~/.Renviron` file or temporarily using `Sys.setenv()`. The three `GS_KEY` variables are Google Sheets keys that can be extracted from the relevant URLs. `EBIRDST_S3_BUCKET` is the `s3://` URL for the S3 bucket containing the eBird Status and Trends data products. To complete the steps outlined in this document you will need read access to the S3 bucket.

```
EBIRDST_S3_BUCKET=
EBIRDST_STATUS_GS_KEY=
EBIRDST_TRENDS_GS_KEY=
EBIRdST_FEATUERS_GS_KEY=
```

## Update version years in R function

Update `status_version_year` and `trends_version_year` in the function definition for `ebirdst_version()` in `R/download.R`. For Status, this is the prediction year and for Trends it's the final year of the trend time series. If only one of Status or Trends is being updated in a given year, only update the year for the products that have been updated.

Rebuild the R package with `devtools::install()` prior to proceeding since later steps rely on these version years.

## Update citation

Update `R/zzz.R` to reflect the updated citation for this release.

## Update the change log

Each release Tom Auer compiles a detailed change log for anything that has changed since the last release. This includes everything from updated covariates to base model changes to difference in released data products. Tom will typically share this change log as a Google Document, which should be proofread, converted to markdown, and added to the change log vignette in `vignettes/articles/product-changelog.qmd`.
