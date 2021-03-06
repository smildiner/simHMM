% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_jasper.R
\name{credible_interval}
\alias{credible_interval}
\title{Compute upper and lower values of the 95\% credible interval}
\usage{
credible_interval(x, type = c("0.95", "0.99"))
}
\arguments{
\item{x}{posterior values for a parameter}

\item{type}{either one of '0.95' for 95\% CI or '0.99' for 99\% CI}
}
\value{
numeric vector. Element 1 is the lower 95\% CI, element 2 is the upper 95\% CI.
}
\description{
Compute upper and lower values of the 95\% credible interval
}
