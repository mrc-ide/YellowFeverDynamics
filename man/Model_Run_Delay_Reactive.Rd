% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{Model_Run_Delay_Reactive}
\alias{Model_Run_Delay_Reactive}
\title{Model_Run_Delay_Reactive}
\usage{
Model_Run_Delay_Reactive(
  FOI_spillover = c(),
  R0 = c(),
  vacc_data = list(),
  pop_data = list(),
  years_data = c(1940:1941),
  year0 = 1940,
  vaccine_efficacy = 1,
  time_inc = 1,
  mode_start = 1,
  start_SEIRV = list(),
  mode_time = 0,
  n_particles = 1,
  n_threads = 1,
  deterministic = FALSE,
  response_delay = 56,
  p_rep = c(0, 0),
  case_threshold = Inf,
  cluster_threshold = Inf,
  vacc_cov_cam = c(),
  t_cam = 0
)
}
\arguments{
\item{FOI_spillover}{Vector of values of force of infection due to spillover from sylvatic reservoir
(size depends on mode_time)}

\item{R0}{Vector of values of basic reproduction number for urban spread of infection (size depends on mode_time)}

\item{vacc_data}{Projected vaccination-based immunity (assuming vaccine_efficacy = 1) by age group and year}

\item{pop_data}{Population by age group and year}

\item{years_data}{Incremental vector of years denoting years for which to save data}

\item{year0}{First year in population/vaccination data}

\item{vaccine_efficacy}{Proportional vaccine efficacy}

\item{time_inc}{Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)}

\item{mode_start}{Flag indicating how to set initial population immunity level in addition to vaccination \cr
If mode_start = 0, only vaccinated individuals \cr
If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age) \cr
If mode_start = 2, use SEIRV input in list from previous run(s)}

\item{start_SEIRV}{SEIRV data (including E_delay and I_delay) from end of a previous run to use as input}

\item{mode_time}{Type of time dependence of FOI_spillover and R0 to be used: \cr
If mode_time = 0, no time variation (constant values)\cr
If mode_time = 1, FOI/R0 vary annually without seasonality (number of values = number of years to consider) \cr
If mode_time = 2, FOI/R0 vary with monthly seasonality without inter-annual variation (number of values = 12) \cr
If mode_time = 3, FOI/R0 vary with daily seasonality without inter-annual variation (number of values = 365/dt) \cr
If mode_time = 4, FOI/R0 vary annually with monthly seasonality (number of values = 12 * number of years to consider) \cr
If mode_time = 5, FOI/R0 vary annually with daily seasonality (number of values = (365/dt) * number of years to consider)}

\item{n_particles}{number of particles to use}

\item{n_threads}{number of threads to use}

\item{deterministic}{TRUE/FALSE - set model to run in deterministic mode if TRUE}

\item{response_delay}{Delay time in days between a threshold being reached and emergency conditions coming into effect}

\item{p_rep}{Probabilities of an infection being reported as a case before emergency conditions triggered (1st value) or
after emergency conditions triggered (2nd value)}

\item{case_threshold}{Threshold total no. reported cases to trigger emergency conditions}

\item{cluster_threshold}{Threshold current infectious fraction to trigger emergency conditions}

\item{vacc_cov_cam}{Target vaccination coverage by age group during emergency campaign}

\item{t_cam}{Duration in days of emergency vaccination campaign
'}
}
\description{
Runs delay+reactive version of SEIRV model
}
\details{
Accepts epidemiological + population parameters and model settings; runs delay/reactive SEIRV model
for one region over a specified time period for a number of particles/threads and outputs time-dependent SEIRV
values, infection numbers and total force of infection values. This version of the model differs from the standard
one in simulating an emergency vaccination campaign applied when an outbreak is declared (as well as using delay
instead of rate for incubation, infectious period etc.). Case reporting is governed by an additional parameter
p_rep which can also change after a reported outbreak is triggered in order to reflect changes in surveillance.
An outbreak is declared when the number of reported cases or the infected fraction of the population exceed
supplied thresholds.
}
