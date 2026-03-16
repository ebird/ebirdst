# Store the eBird Status and Trends access key

Accessing eBird Status and Trends data requires an access key, which can
be obtained by visiting https://ebird.org/st/request. This key must be
stored as the environment variable `EBIRDST_KEY` in order for
[`ebirdst_download_status()`](https://ebird.github.io/ebirdst/reference/ebirdst_download_status.md)
and
[`ebirdst_download_trends()`](https://ebird.github.io/ebirdst/reference/ebirdst_download_trends.md)
to use it. The easiest approach is to store the key in your `.Renviron`
file so it can always be accessed in your R sessions. Use this function
to set `EBIRDST_KEY` in your `.Renviron` file provided that it is
located in the standard location in your home directory. It is also
possible to manually edit the `.Renviron` file. **The access key is
specific to you and should never be shared or made publicly
accessible.**

## Usage

``` r
set_ebirdst_access_key(key, overwrite = FALSE)
```

## Arguments

- key:

  character; API key obtained by filling out the form at
  https://ebird.org/st/request.

- overwrite:

  logical; should the existing `EBIRDST_KEY` be overwritten if it has
  already been set in .Renviron.

## Value

Edits .Renviron, then returns the path to this file invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
# save the api key, replace XXXXXX with your actual key
set_ebirdst_access_key("XXXXXX")
} # }
```
