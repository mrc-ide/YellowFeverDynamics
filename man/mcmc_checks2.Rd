% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc2.R
\name{mcmc_checks2}
\alias{mcmc_checks2}
\title{mcmc_checks2}
\usage{
mcmc_checks2(
  log_params_ini = c(),
  n_regions = 1,
  type = NULL,
  prior_settings = list(type = "zero"),
  enviro_data = NULL,
  R0_fixed_values = NULL,
  add_values = list(vaccine_efficacy = 1, p_rep_severe = 1, p_rep_death = 1, m_FOI_Brazil
    = 1),
  extra_estimated_params = list()
)
}
\arguments{
\item{log_params_ini}{Initial values of parameters (natural logarithm of actual parameters)}

\item{n_regions}{Number of regions}

\item{type}{Type of parameter set (FOI only, FOI+R0, FOI and/or R0 coefficients associated with environmental
covariates); choose from "FOI","FOI+R0","FOI enviro","FOI+R0 enviro"}

\item{prior_settings}{TBA}

\item{enviro_data}{Values of environmental covariates (if in use)}

\item{R0_fixed_values}{Values of R0 to use if type set to "FOI" or "FOI enviro"}

\item{add_values}{TBA}

\item{extra_estimated_params}{TBA}
}
\description{
Perform checks on MCMC inputs
}
\details{
This function, which is called by MCMC(), performs a number of checks on data to be used in fitting to
ensure proper functionality. It verifies that the number of parameters being estimated is consistent with other
settings and that certain values are not outwith sensible boundaries (e.g. probabilities must be between 0 and 1).
}
