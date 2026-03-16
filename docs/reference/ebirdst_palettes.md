# eBird Status and Trends color palettes for mapping

Generate the color palettes used for the eBird Status and Trends
relative abundance and trends maps.

## Usage

``` r
ebirdst_palettes(
  n,
  type = c("weekly", "breeding", "nonbreeding", "migration", "prebreeding_migration",
    "postbreeding_migration", "year_round", "trends")
)
```

## Arguments

- n:

  integer; the number of colors to be in the palette.

- type:

  character; the type of color palette: "weekly" for the weekly relative
  abundance, "trends" for trends color palette, and a season name for
  the seasonal relative abundance. Note that for trends a diverging
  palette is returned, while all other palettes are sequential.

## Value

A character vector of hex color codes.

## Examples

``` r
# breeding season color palette
ebirdst_palettes(10, type = "breeding")
#>  [1] "#DFC0BC" "#DBADA7" "#D89A92" "#D5887D" "#D27568" "#CF6252" "#CC503E"
#>  [8] "#BB4938" "#AA4233" "#993C2E"
```
