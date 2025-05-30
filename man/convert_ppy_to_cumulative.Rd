% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trends.R
\name{convert_ppy_to_cumulative}
\alias{convert_ppy_to_cumulative}
\title{Convert percent per year trend to cumulative trend}
\usage{
convert_ppy_to_cumulative(x, n_years)
}
\arguments{
\item{x}{numeric; percent per year trend on the 0-100 scale rather than the
0-1 scale.}

\item{n_years}{integer; number of years.}
}
\value{
A numeric vector of the same length as \code{x} that contains the
cumulative trend resulting from \code{n_years} years of compounding annual
trend.
}
\description{
Convert percent per year trend to cumulative trend
}
\examples{
ppy_trend <- runif(100, min = -100, 100)
cumulative_trend <- convert_ppy_to_cumulative(ppy_trend, n_years = 5)
cbind(ppy_trend, cumulative_trend)
}
