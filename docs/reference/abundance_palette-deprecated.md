# eBird Status and Trends color palettes for mapping

This deprecated function has been replaced by
[`ebirdst_palettes`](https://ebird.github.io/ebirdst/reference/ebirdst_palettes.md).
Both functions generate color palettes used for the eBird Status and
Trends relative abundance maps.

## Usage

``` r
abundance_palette(n,
                        season = c("weekly", "breeding",
                                   "nonbreeding",
                                   "migration",
                                   "prebreeding_migration",
                                   "postbreeding_migration",
                                   "year_round"))
```

## Arguments

- n:

  integer; the number of colors to be in the palette.

- season:

  character; the season to generate colors for or "weekly" to get the
  color palette used in the weekly abundance animations.

## Value

A character vector of hex color codes.

## See also

[`ebirdst_palettes`](https://ebird.github.io/ebirdst/reference/ebirdst_palettes.md)
[`ebirdst-deprecated`](https://ebird.github.io/ebirdst/reference/ebirdst-deprecated.md)
