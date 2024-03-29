% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trends.R
\name{rasterize_trends}
\alias{rasterize_trends}
\title{Convert eBird Trends Data Products to raster format}
\usage{
rasterize_trends(
  trends,
  layers = c("abd_ppy", "abd_ppy_lower", "abd_ppy_upper"),
  trim = TRUE
)
}
\arguments{
\item{trends}{data frame; trends data for a single species as returned by
\code{\link[=load_trends]{load_trends()}}.}

\item{layers}{character; column names in the trends data frame to rasterize.
These columns will become layers in the raster that is created.}

\item{trim}{logical; flag indicating if the returned raster should be trimmed
to remove outer rows and columns that are NA. If \code{trim = FALSE} the
returned raster will have a global extent, which can be useful if rasters
will be combined across species with different ranges.}
}
\value{
A \link[terra:SpatRaster-class]{SpatRaster} object.
}
\description{
The eBird trends data are stored in a tabular format, where each row gives
the trend estimate for a single cell in a 27 km x 27 km equal area grid. For
many applications, an explicitly spatial format is more useful. This function
uses the cell center coordinates to convert the tabular trend estimates to
raster format in \code{terra} \link[terra:SpatRaster-class]{SpatRaster} format.
}
\examples{
\dontrun{
# download example trends data if it hasn't already been downloaded
ebirdst_download_trends("yebsap-example")

# load trends
trends <- load_trends("yebsap-example")

# rasterize percent per year trend
rasterize_trends(trends, "abd_ppy")
}
}
