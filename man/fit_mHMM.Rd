% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_mHMM.R
\name{fit_mHMM}
\alias{fit_mHMM}
\title{Fit mHMM using mHMMfast}
\usage{
fit_mHMM(
  m,
  n_dep,
  q_emiss,
  gamma,
  emiss,
  iter,
  burnin,
  start_gamma = NULL,
  start_emiss = NULL,
  data_sim,
  light = FALSE,
  save_path = FALSE,
  save_subj_data = TRUE,
  progress_line = FALSE
)
}
\description{
Fit mHMM using mHMMfast
}
