% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outputs_generate.R
\name{Generate_Dataset_VarFR}
\alias{Generate_Dataset_VarFR}
\title{Generate_Dataset_VarFR}
\usage{
Generate_Dataset_VarFR(
  input_data = list(),
  FOI_values = c(),
  R0_values = c(),
  sero_template = NULL,
  case_template = NULL,
  vaccine_efficacy = 1,
  p_severe_inf = 0.12,
  p_death_severe_inf = 0.39,
  p_rep_severe = 1,
  p_rep_death = 1,
  mode_start = 1,
  start_SEIRV = NULL,
  dt = 1,
  n_reps = 1,
  deterministic = FALSE,
  mode_time = 1,
  output_frame = FALSE
)
}
\arguments{
\item{input_data}{List of population and vaccination data for multiple regions in standard format [TBA]}

\item{FOI_values}{Array of values of the force of infection due to spillover from sylvatic reservoir by region and year}

\item{R0_values}{Array of values of the basic reproduction number for human-human transmission by region and year}

\item{sero_template}{Seroprevalence data template - data frame with region, year, minimum/maximum age, vc_factor [TBA]
and number of samples}

\item{case_template}{Annual reported case/death data template - data frame with region and year}

\item{vaccine_efficacy}{Fractional vaccine efficacy}

\item{p_severe_inf}{Probability of an infection being severe}

\item{p_death_severe_inf}{Probability of a severe infection resulting in death}

\item{p_rep_severe}{Probability of reporting of a severe but non-fatal infection}

\item{p_rep_death}{Probability of reporting of a fatal infection}

\item{mode_start}{Flag indicating how to set initial population immunity level in addition to vaccination
If mode_start=0, only vaccinated individuals
If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity (uniform by age, R0 based only)
If mode_start=2, use SEIRV input in list from previous run(s)
If mode_start=3, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age)}

\item{start_SEIRV}{SEIRV data from end of a previous run to use as input (list of datasets, one per region)}

\item{dt}{Time increment in days to use in model (should be either 1.0, 2.5 or 5.0 days)}

\item{n_reps}{number of stochastic repetitions}

\item{deterministic}{TRUE/FALSE - set model to run in deterministic mode if TRUE}

\item{mode_time}{Type of time dependence of FOI_spillover and R0 to be used: If mode_time=0, no time
  variation (identical to Model_Run())
  If mode_time=1, FOI/R0 vary annually without seasonality (number of values = number of years to consider)
  If mode_time=2, FOI/R0 vary with monthly seasonality without inter-annual variation (number of values = 12)
  If mode_time=3, FOI/R0 vary with daily seasonality without inter-annual variation (number of values = 365/dt)
  If mode_time=4, FOI/R0 vary annually with monthly seasonality (number of values = 12.number of years to consider)
  If mode_time=5, FOI/R0 vary annually with daily seasonality (number of values = (365/dt)*number of years to consider)
'}

\item{output_frame}{Flag indicating whether to output a complete data frame of results in template format (if TRUE)
or calculated values only (if FALSE)}
}
\description{
Generate annual serological and/or case/death data with annually varying FOI/R0
}
\details{
This function is used to generate annual serological and/or case/death data based on templates;
it is a variant of the Generate_Dataset function in the YEP package, identical except for using the
Model_Run_VarFR function instead of the YEP::Model_Run function. Inputs are the same as YEP::Generate_Dataset
except that FOI_values and R0_values are 2-dimensional arrays (dimensions no. regions x no. years). (Parallel
processing functionality has also been removed)
}
