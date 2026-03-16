# Get eBird species code for a set of species

Give a vector of species codes, common names, and/or scientific names,
return a vector of 6-letter eBird species codes. This function will only
look up codes for species for which eBird Status and Trends results
exist.

## Usage

``` r
get_species(x)
```

## Arguments

- x:

  character; vector of species codes, common names, and/or scientific
  names.

## Value

A character vector of eBird species codes.

## Examples

``` r
get_species(c("Black-capped Chickadee", "Poecile gambeli", "carchi"))
#> [1] "bkcchi" "mouchi" "carchi"
```
