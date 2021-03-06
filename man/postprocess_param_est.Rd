% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reshape_results.R
\name{postprocess_param_est}
\alias{postprocess_param_est}
\title{Postprocessing utility function for parameter estimates}
\usage{
postprocess_param_est(z, m)
}
\arguments{
\item{z}{list. contains the MAP values of the parameter estimates.}

\item{m}{integer. Number of hidden states}
}
\value{
data frame in which each parameter estimate has been placed in its own column.
}
\description{
Postprocessing utility function for parameter estimates
}
