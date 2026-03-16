# Get the Status and Trends week that a date falls into

Get the Status and Trends week that a date falls into

## Usage

``` r
date_to_st_week(dates, version = 2022)
```

## Arguments

- dates:

  a vector of dates.

- version:

  One of `2021` for the date scheme used for the 2021 and prior data
  releases or `2022` for the date scheme used in the 2022 and subsequent
  releases. Default is `2022`.

## Value

An integer vector of weeks numbers from 1-52.

## Examples

``` r
d <- as.Date(c("2016-04-08", "2018-12-31", "2014-01-01", "2018-09-04"))
date_to_st_week(d)
#> [1] 15 52  1 36
```
