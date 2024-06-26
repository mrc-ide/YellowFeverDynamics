% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc2.R
\name{single_posterior_calc2}
\alias{single_posterior_calc2}
\title{single_posterior_calc2}
\usage{
single_posterior_calc2(
  log_params_prop = c(),
  input_data = list(),
  obs_sero_data = NULL,
  obs_case_data = NULL,
  consts = list()
)
}
\arguments{
\item{log_params_prop}{Proposed values of parameters (natural logarithm of actual parameters)}

\item{input_data}{List of population and vaccination data for multiple regions (created using data input
creation code and usually loaded from RDS file), with cross-reference tables added using input_data_process
in MCMC}

\item{obs_sero_data}{Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
positives}

\item{obs_case_data}{Annual reported case/death data for comparison, by region and year, in format no. cases/no.
deaths}

\item{consts}{= List of constant parameters/flags/etc. loaded to mcmc() (type,
mode_start,prior_settings,dt,n_reps,enviro_data,R0_fixed_values,p_severe_inf,
p_death_severe_inf,add_values list,extra_estimated_params,deterministic, mode_parallel, cluster)}
}
\description{
Function which calculates and outputs posterior likelihood of observing simulated data
}
\details{
This function calculates the posterior likelihood of observing a set of observations (across multiple
regions and data types) for a given proposed parameter set. [TBA]
}
