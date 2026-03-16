# Defunct functions in package ebirdst.

The functions listed below are defunct and no longer supported. Calling
them will result in an error.

When possible alternative functions are suggested.

Many of them supported stixles which were infrequently used and were
dropped from ebirdst with the 2022 data release.

## Usage

``` r
ebirdst_download(
  species,
  path = ebirdst_data_dir(),
  tifs_only = TRUE,
  force = FALSE,
  show_progress = TRUE,
  pattern = NULL,
  dry_run = FALSE
)

ebirdst_extent(x, t, ...)

ebirdst_habitat(path, ext, data = NULL, stationary_associations = FALSE)

ebirdst_ppms(path, ext, es_cutoff, pat_cutoff)

ebirdst_ppms_ts(ath, ext, summarize_by = c("weeks", "months"), ...)

ebirdst_subset(x, crs)

load_pds(path, ext, model = c("occurrence", "count"), return_sf = FALSE)

load_pis(path, ext, model = c("occurrence", "count"), return_sf = FALSE)

load_predictions(path, return_sf = FALSE)

parse_raster_dates(x)

load_stixels(path, ext, return_sf = FALSE)

project_extent(x, crs)

plot_pds(path, ext, summarize_by = c("weeks", "months"), ...)

plot_pis(
  pis,
  ext,
  by_cover_class = TRUE,
  n_top_pred = 15,
  pretty_names = TRUE,
  plot = TRUE
)

stixelize(x)
```

## Arguments

- ...:

  All arguments are now ignored.
