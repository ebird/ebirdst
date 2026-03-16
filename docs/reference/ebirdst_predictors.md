# eBird Status and Trends predictor variables

A data frame of the predictors used in the eBird Status and Trends
models. These include effort variables (e.g. distance traveled, number
of observers, etc.) in addition to variables describing the environment
(e.g. elevation, land cover, water cover, etc.). The environmental
variables are derived by summarizing remotely sensed datasets (described
in
[ebirdst_predictor_descriptions](https://ebird.github.io/ebirdst/reference/ebirdst_predictor_descriptions.md))
over a 3 km diameter neighborhood around each checklist. For categorical
datasets, two variables are generated for each class describing the
percent cover (`pland`) and edge density (`ed`).

## Usage

``` r
ebirdst_predictors
```

## Format

A data frame with 150 rows and 4 columns:

- `predictor`: predictor name.

- `dataset`: dataset name, which can be cross referenced in
  [ebirdst_predictor_descriptions](https://ebird.github.io/ebirdst/reference/ebirdst_predictor_descriptions.md)
  for further details.

- `class`: class number or name for categorical variables.

- `label`: descriptive labels for each predictor variable.
