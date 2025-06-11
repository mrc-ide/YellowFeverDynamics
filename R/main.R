# R file for general functions in YellowFeverDynamics package
#-------------------------------------------------------------------------------
#Global variables
t_incubation <- 5 #Time for cases to incubate in mosquito
t_latent <- 5 #Latent period before cases become infectious
t_infectious <- 5 #Time cases remain infectious
#-------------------------------------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.
#' @useDynLib YellowFeverDynamics, .registration = TRUE
#' @importFrom assertthat assert_that
#' @import dde
#' @import dust2
#' @importFrom graphics axis matplot par
#' @importFrom mvtnorm rmvnorm
#' @import odin2
#' @import parallel
#' @importFrom R.utils fileAccess
#' @importFrom stats cov dexp dnbinom prop.test rbinom runif dnorm
#' @importFrom tgp lhs
#' @importFrom truncdist dtrunc
#' @importFrom utils write.csv
#' @import YEP
#-------------------------------------------------------------------------------
# unload DLL when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("YellowFeverDynamics", libpath)
}
#-------------------------------------------------------------------------------
#' @title Model_Run_Delay
#'
#' @description Run SEIRV model version using time delay instead of rate to move individuals from E-I and I-R
#'
#' @details Accepts epidemiological + population parameters and model settings; runs version of SEIRV model which
#' uses a fixed time delay instead of a rate to move individuals from exposed (E) to infectious (I) and from
#' infectious to recovered (R). The model is run for one region over a specified time period for a number of
#' particles/threads and outputs time-dependent SEIRV values, infection numbers and total force of infection values.
#'
#' @param FOI_spillover Vector of values of force of infection due to spillover from sylvatic reservoir
#'   (size depends on mode_time)
#' @param R0 Vector of values of basic reproduction number for urban spread of infection (size depends on mode_time)
#' @param vacc_data Projected vaccination-based immunity (assuming vaccine_efficacy = 1) by age group and year
#' @param pop_data Population by age group and year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param year0 First year in population/vaccination data
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param time_inc Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination \cr
#'  If mode_start = 0, only vaccinated individuals \cr
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age) \cr
#'  If mode_start = 2, use SEIRV input in list from previous run(s)
#' @param start_SEIRV SEIRV data (including E_delay and I_delay) from end of a previous run to use as input
#' @param mode_time Type of time dependence of FOI_spillover and R0 to be used: \cr
#'  If mode_time = 0, no time variation (constant values)\cr
#'  If mode_time = 1, FOI/R0 vary annually without seasonality (number of values = number of years to consider) \cr
#'  If mode_time = 2, FOI/R0 vary with monthly seasonality without inter-annual variation (number of values = 12) \cr
#'  If mode_time = 3, FOI/R0 vary with daily seasonality without inter-annual variation (number of values = 365/dt) \cr
#'  If mode_time = 4, FOI/R0 vary annually with monthly seasonality (number of values = 12 * number of years to consider) \cr
#'  If mode_time = 5, FOI/R0 vary annually with daily seasonality (number of values = (365/dt) * number of years to consider)
#' @param n_particles number of particles to use
#' @param n_threads number of threads to use
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' '
#' @export
#'
Model_Run_Delay <- function(FOI_spillover = c(), R0 = c(), vacc_data = list(), pop_data = list(), years_data = c(1940:1941),
                            year0 = 1940, vaccine_efficacy = 1.0, time_inc = 1.0, mode_start = 1, start_SEIRV = list(),
                            mode_time = 0, n_particles = 1, n_threads = 1, deterministic = FALSE) {

  #TODO Add assert_that functions (NB - Some checks carried out in parameter_setup)
  assert_that(n_particles <=  20, msg = "Number of particles must be 20 or less")

  N_age = length(pop_data[1, ]) #Number of age groups
  nd1 <- (t_incubation + t_latent)/time_inc
  nd2 <- t_infectious/time_inc
  step_begin = ((years_data[1] - year0) * (365/time_inc)) #Step at which data starts being saved for final output
  step_end = ((max(years_data) + 1 - year0) * (365/time_inc)) - 1 #Step at which to end
  t_pts_out = step_end - step_begin + 1 #Number of time points in final output data

  pars = parameter_setup(FOI_spillover, R0, vacc_data, pop_data, years_data, year0,
                         vaccine_efficacy, time_inc, mode_start, start_SEIRV, mode_time)

  #Carrying forward delay from previous run may cause errors
  #if(mode_start == 2){pars$E_delay0 = start_SEIRV$E_delay} else {pars$E_delay0 = rep(0, nd1 * N_age)}
  pars$np_E_delay = nd1 * N_age
  pars$np_I_delay = nd2 * N_age
  pars$E_delay0 = rep(0, pars$np_E_delay)
  pars$I_delay0 = rep(0, pars$np_I_delay)

  x <- dust_system_create(SEIRVModelDelay, pars = pars,
                          n_particles = n_particles, n_threads = n_threads, time = 0, dt = 1,
                          deterministic = deterministic, preserve_particle_dimension = TRUE)
  index = dust_unpack_index(x)
  dust_system_set_state_initial(x)
  x_res <- dust_system_simulate(x, c(step_begin:step_end))

  dimensions = c(N_age, n_particles, t_pts_out)
  output_data = list(day = x_res[1, 1, ], year = x_res[2, 1, ])
  output_data$FOI_total = array(x_res[3, , ]/time_inc, dim = c(n_particles, t_pts_out))
  output_data$S = array(x_res[index$S, , ], dim = dimensions)
  output_data$E = array(x_res[index$E, , ], dim = dimensions)
  output_data$E_delay = array(x_res[index$E_delay, , ], dim = c(N_age, nd1, n_particles, t_pts_out))
  output_data$I = array(x_res[index$I, , ], dim = dimensions)
  output_data$I_delay = array(x_res[index$I_delay, , ], dim = c(N_age, nd2, n_particles, t_pts_out))
  output_data$R = array(x_res[index$R, , ], dim = dimensions)
  output_data$V = array(x_res[index$V, , ], dim = dimensions)
  output_data$C = array(x_res[index$C, , ], dim = dimensions)

  return(output_data)
}
#-------------------------------------------------------------------------------
#' @title Model_Run_Delay_Reactive
#'
#' @description Runs delay+reactive version of SEIRV model
#'
#' @details Accepts epidemiological + population parameters and model settings; runs delay/reactive SEIRV model
#' for one region over a specified time period for a number of particles/threads and outputs time-dependent SEIRV
#' values, infection numbers and total force of infection values. This version of the model differs from the standard
#' one in simulating an emergency vaccination campaign applied when an outbreak is declared (as well as using delay
#' instead of rate for incubation, infectious period etc.). Case reporting is governed by an additional parameter
#' p_rep which can also change after a reported outbreak is triggered in order to reflect changes in surveillance.
#' An outbreak is declared when the number of reported cases or the infected fraction of the population exceed
#' supplied thresholds.
#'
#' @param FOI_spillover Vector of values of force of infection due to spillover from sylvatic reservoir
#'   (size depends on mode_time)
#' @param R0 Vector of values of basic reproduction number for urban spread of infection (size depends on mode_time)
#' @param vacc_data Projected vaccination-based immunity (assuming vaccine_efficacy = 1) by age group and year
#' @param pop_data Population by age group and year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param year0 First year in population/vaccination data
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param time_inc Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination \cr
#'  If mode_start = 0, only vaccinated individuals \cr
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age) \cr
#'  If mode_start = 2, use SEIRV input in list from previous run(s)
#' @param start_SEIRV SEIRV data (including E_delay and I_delay) from end of a previous run to use as input
#' @param mode_time Type of time dependence of FOI_spillover and R0 to be used: \cr
#'  If mode_time = 0, no time variation (constant values)\cr
#'  If mode_time = 1, FOI/R0 vary annually without seasonality (number of values = number of years to consider) \cr
#'  If mode_time = 2, FOI/R0 vary with monthly seasonality without inter-annual variation (number of values = 12) \cr
#'  If mode_time = 3, FOI/R0 vary with daily seasonality without inter-annual variation (number of values = 365/dt) \cr
#'  If mode_time = 4, FOI/R0 vary annually with monthly seasonality (number of values = 12 * number of years to consider) \cr
#'  If mode_time = 5, FOI/R0 vary annually with daily seasonality (number of values = (365/dt) * number of years to consider)
#' @param n_particles number of particles to use
#' @param n_threads number of threads to use
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' @param response_delay Delay time in days between a threshold being reached and emergency conditions coming into effect
#' @param p_rep Probabilities of an infection being reported as a case before emergency conditions triggered (1st value) or
#'   after emergency conditions triggered (2nd value)
#' @param case_threshold Threshold total no. reported cases to trigger emergency conditions
#' @param cluster_threshold Threshold current infectious fraction to trigger emergency conditions
#' @param vacc_cov_cam Target vaccination coverage by age group during emergency campaign
#' @param t_cam Duration in days of emergency vaccination campaign
#' '
#' @export
#'
Model_Run_Delay_Reactive <- function(FOI_spillover = c(), R0 = c(), vacc_data = list(), pop_data = list(), years_data = c(1940:1941),
                                     year0 = 1940, vaccine_efficacy = 1.0, time_inc = 1.0, mode_start = 1, start_SEIRV = list(),
                                     mode_time = 0, n_particles = 1, n_threads = 1, deterministic = FALSE,
                                     response_delay = 56.0, p_rep = c(0.0, 0.0), case_threshold = Inf,
                                     cluster_threshold = Inf, vacc_cov_cam = c(), t_cam = 0) {

  #TODO Add assert_that functions (NB - Some checks carried out in parameter_setup)
  assert_that(n_particles <=  20, msg = "Number of particles must be 20 or less")

  N_age = length(pop_data[1, ]) #Number of age groups
  assert_that(length(vacc_cov_cam) == N_age)
  nd1 <- (t_incubation + t_latent)/time_inc
  nd2 <- t_infectious/time_inc
  nd <- nd1 + nd2
  step_begin = ((years_data[1] - year0) * (365/time_inc)) #Step at which data starts being saved for final output
  step_end = ((max(years_data) + 1 - year0) * (365/time_inc)) - 1 #Step at which to end
  t_pts_out = step_end - step_begin + 1 #Number of time points in final output data

  pars1 = parameter_setup(FOI_spillover, R0, vacc_data, pop_data, years_data, year0,
                          vaccine_efficacy, time_inc, mode_start, start_SEIRV, mode_time)
  n_years = length(pop_data[, 1]) - 1
  inv_365 = 1.0/365.0
  pars2 = list(FOI_spillover = pars1$FOI_spillover, R0 = pars1$R0, vacc_rate_daily = pars1$vacc_rate_daily,
               vacc_cov_cam = vacc_cov_cam, t_cam = t_cam, N_age = pars1$N_age,
               S_0 = pars1$S_0, E_0 = pars1$E_0, I_0 = pars1$I_0, R_0 = pars1$R_0, V_0 = pars1$V_0,
               dP1_all = pars1$dP1_all, dP2_all = pars1$dP2_all, n_years = pars1$n_years, year0 = pars1$year0,
               vaccine_efficacy = pars1$vaccine_efficacy, time_inc = pars1$time_inc, t_incubation = pars1$t_incubation,
               t_latent = pars1$t_latent, t_infectious = pars1$t_infectious, n_t_pts = pars1$n_t_pts,
               response_delay = response_delay, p_rep = p_rep, case_threshold = case_threshold,
               cluster_threshold = cluster_threshold)

  #Carrying forward delay from previous run may cause errors
  #if(mode_start == 2){pars2$E_delay0 = start_SEIRV$E_delay} else {pars2$E_delay0 = rep(0, nd1 * N_age)}
  pars2$np_E_delay = nd1 * N_age
  pars2$np_I_delay = nd2 * N_age
  pars2$E_delay0 = rep(0, pars2$np_E_delay)
  pars2$I_delay0 = rep(0, pars2$np_I_delay)

  #Check that there is no overlap between emergency campaign and other vaccination
  #(Not yet possible to adjust vaccine rates on the fly to deal with overlap)
  for(i in 1:N_age){
    if(vacc_cov_cam[i]>0){assert_that(all(pars2$vacc_rate_daily[i, ] == 0))}
  }

  x <- dust_system_create(SEIRVModelDelayReactive, pars = pars2,
                          n_particles = n_particles, n_threads = n_threads, time = 0, dt = 1,
                          deterministic = deterministic, preserve_particle_dimension = TRUE)
  index = dust_unpack_index(x)
  dust_system_set_state_initial(x)
  x_res <- dust_system_simulate(x, c(step_begin:step_end))

  dims1 = c(n_particles, t_pts_out)
  dims2 = c(N_age, n_particles, t_pts_out)
  output_data = list(day = x_res[1, 1, ], year = x_res[2, 1, ], FOI_total = array(x_res[3, , ]/time_inc, dim = dims1),
                     C_rep_total = array(x_res[4, , ], dim = dims1),
                     flag1 = array(x_res[5, , ], dim = dims1), flag2 = array(x_res[6, , ], dim = dims1),
                     flag3 = array(x_res[7, , ], dim = dims1), flag4 = array(x_res[8, , ], dim = dims1),
                     report_rate = array(x_res[9, , ], dim = dims1))
  output_data$S = array(x_res[index$S, , ], dim = dims2)
  output_data$E = array(x_res[index$E, , ], dim = dims2)
  output_data$E_delay = array(x_res[index$E_delay, , ], dim = c(N_age, nd1, n_particles, t_pts_out))
  output_data$I = array(x_res[index$I, , ], dim = dims2)
  output_data$I_delay = array(x_res[index$I_delay, , ], dim = c(N_age, nd2, n_particles, t_pts_out))
  output_data$R = array(x_res[index$R, , ], dim = dims2)
  output_data$V = array(x_res[index$V, , ], dim = dims2)
  output_data$C = array(x_res[index$C, , ], dim = dims2)
  #output_data$C_rep = array(x_res[index$C_rep, , ], dim = dims2)

  return(output_data)
}
#-------------------------------------------------------------------------------
#' @title Model_Run_Split
#'
#' @description Run full SEIRV model with daily infection output split into sylvatic and urban
#'
#' @details Accepts epidemiological + population parameters and model settings; runs split-infection SEIRV model
#' for one region over a specified time period for a number of particles/threads and outputs time-dependent SEIRV
#' values, infection numbers and total force of infection values. This version of the model differs from the standard
#' one in that infections are split into sylvatic and urban, allowing the relative importance of sylvatic and urban
#' infections to be assessed.
#'
#' @param FOI_spillover Vector of values of force of infection due to spillover from sylvatic reservoir
#'   (size depends on mode_time)
#' @param R0 Vector of values of basic reproduction number for urban spread of infection (size depends on mode_time)
#' @param vacc_data Projected vaccination-based immunity (assuming vaccine_efficacy = 1) by age group and year
#' @param pop_data Population by age group and year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param year0 First year in population/vaccination data
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param time_inc Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination \cr
#'  If mode_start = 0, only vaccinated individuals \cr
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age) \cr
#'  If mode_start = 2, use SEIRV input in list from previous run(s)
#' @param start_SEIRV SEIRV data (including E_delay and I_delay) from end of a previous run to use as input
#' @param mode_time Type of time dependence of FOI_spillover and R0 to be used: \cr
#'  If mode_time = 0, no time variation (constant values)\cr
#'  If mode_time = 1, FOI/R0 vary annually without seasonality (number of values = number of years to consider) \cr
#'  If mode_time = 2, FOI/R0 vary with monthly seasonality without inter-annual variation (number of values = 12) \cr
#'  If mode_time = 3, FOI/R0 vary with daily seasonality without inter-annual variation (number of values = 365/dt) \cr
#'  If mode_time = 4, FOI/R0 vary annually with monthly seasonality (number of values = 12 * number of years to consider) \cr
#'  If mode_time = 5, FOI/R0 vary annually with daily seasonality (number of values = (365/dt) * number of years to consider)
#' @param n_particles number of particles to use
#' @param n_threads number of threads to use
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' '
#' @export
#'
Model_Run_Split <- function(FOI_spillover = c(), R0 = c(), vacc_data = list(), pop_data = list(), years_data = c(1940:1941),
                            year0 = 1940, vaccine_efficacy = 1.0, time_inc = 1.0, mode_start = 1, start_SEIRV = list(),
                            mode_time = 0, n_particles = 1, n_threads = 1, deterministic = FALSE) {

  #TODO Add assert_that functions (NB - Some checks carried out in parameter_setup)
  assert_that(n_particles <=  20, msg = "Number of particles must be 20 or less")

  N_age = length(pop_data[1, ]) #Number of age groups
  step_begin = ((years_data[1] - year0) * (365/time_inc)) #Step at which data starts being saved for final output
  step_end = ((max(years_data) + 1 - year0) * (365/time_inc)) - 1 #Step at which to end
  t_pts_out = step_end - step_begin + 1 #Number of time points in final output data

  x <- dust_system_create(SEIRVModelSplitInfection,
                          pars = parameter_setup(FOI_spillover, R0, vacc_data, pop_data, years_data, year0, vaccine_efficacy,
                                                 time_inc, mode_start, start_SEIRV, mode_time),
                          n_particles = n_particles, n_threads = n_threads, time = 0, dt = 1,
                          deterministic = deterministic, preserve_particle_dimension = TRUE)
  index = dust_unpack_index(x)
  dust_system_set_state_initial(x)
  x_res <- dust_system_simulate(x, c(step_begin:step_end))

  dimensions = c(N_age, n_particles, t_pts_out)
  output_data = list(day = x_res[1, 1, ], year = x_res[2, 1, ],
                     FOI_sylvatic = array(x_res[3, , ]/time_inc, dim = c(n_particles, t_pts_out)),
                     FOI_urban = array(x_res[4, , ]/time_inc, dim = c(n_particles, t_pts_out)))
  output_data$S = array(x_res[index$S, , ], dim = dimensions)
  output_data$E_sylvatic = array(x_res[index$E_sylv, , ], dim = dimensions)
  output_data$E_urban = array(x_res[index$E_urb, , ], dim = dimensions)
  output_data$I_sylvatic = array(x_res[index$I_sylv, , ], dim = dimensions)
  output_data$I_urban = array(x_res[index$I_urb, , ], dim = dimensions)
  output_data$R = array(x_res[index$R, , ], dim = dimensions)
  output_data$V = array(x_res[index$V, , ], dim = dimensions)
  output_data$C_sylvatic = array(x_res[index$C_sylv, , ], dim = dimensions)
  output_data$C_urban = array(x_res[index$C_urb, , ], dim = dimensions)

  return(output_data)
}
#-------------------------------------------------------------------------------
# TODO - Adjust possible output format
#' @title Model_Run_Delay_Many_Reps
#'
#' @description Run delay SEIRV model for single region for large number of repetitions
#'
#' @details Accepts epidemiological + population parameters and model settings; runs delay version of SEIRV model
#' for one region over a specified time period for a number of repetitions and outputs time-dependent SEIRV
#' values, infection numbers and/or total force of infection values. Variation of Model_Run_Delay() used for
#' running a large number of repetitions (>20).
#'
#' @param FOI_spillover Vector of values of force of infection due to spillover from sylvatic reservoir
#'   (size depends on mode_time)
#' @param R0 Vector of values of basic reproduction number for urban spread of infection (size depends on mode_time)
#' @param vacc_data Projected vaccination-based immunity (assuming vaccine_efficacy = 1) by age group and year
#' @param pop_data Population by age group and year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param year0 First year in population/vaccination data
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param time_inc Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination \cr
#'  If mode_start = 0, only vaccinated individuals \cr
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age) \cr
#'  If mode_start = 2, use SEIRV input in list from previous run(s)
#' @param start_SEIRV SEIRV data (including E_delay and I_delay) from end of a previous run to use as input
#' @param mode_time Type of time dependence of FOI_spillover and R0 to be used: \cr
#'  If mode_time = 0, no time variation (constant values)\cr
#'  If mode_time = 1, FOI/R0 vary annually without seasonality (number of values = number of years to consider) \cr
#'  If mode_time = 2, FOI/R0 vary with monthly seasonality without inter-annual variation (number of values = 12) \cr
#'  If mode_time = 3, FOI/R0 vary with daily seasonality without inter-annual variation (number of values = 365/dt) \cr
#'  If mode_time = 4, FOI/R0 vary annually with monthly seasonality (number of values = 12 * number of years to consider) \cr
#'  If mode_time = 5, FOI/R0 vary annually with daily seasonality (number of values = (365/dt) * number of years to consider)
#' @param n_reps Number of repetitions (used to set number of particles and threads)
#' @param division Number of particles/threads to run in one go (up to 20)
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' '
#' @export
#'
Model_Run_Delay_Many_Reps <- function(FOI_spillover = c(), R0 = c(), vacc_data = list(), pop_data = list(),
                                      years_data = c(1940:1941), year0 = 1940, vaccine_efficacy = 1.0, time_inc = 1.0,
                                      mode_start = 1, start_SEIRV = list(), mode_time = 0, n_reps = 1, division = 10,
                                      deterministic = FALSE) {

  assert_that(division <=   20, msg = "Number of particles run at once must be 20 or less")
  n_particles0 = min(division, n_reps)
  n_threads = min(division, n_particles0)
  n_divs = ceiling(n_reps/division)
  if(n_divs ==  1){
    n_particles_list = n_particles0
  } else {
    n_particles_list = c(rep(n_particles0, n_divs - 1), n_reps - (division * (n_divs - 1)))
  }

  N_age = length(pop_data[1, ]) #Number of age groups
  nd1 <- (t_incubation + t_latent)/time_inc
  nd2 <- t_infectious/time_inc
  step_begin = ((years_data[1] - year0) * (365/time_inc)) #Step at which data starts being saved for final output
  step_end = ((max(years_data) + 1 - year0) * (365/time_inc)) - 1 #Step at which to end
  t_pts_out = step_end - step_begin + 1 #Number of time points in final output data

  dim = c(N_age, n_reps, t_pts_out)
  output_data = list(day = rep(NA, t_pts_out), year = rep(NA, t_pts_out), FOI_total = array(NA, c(n_reps, t_pts_out)),
                     S = array(NA, dim), E = array(NA, dim), I = array(NA, dim),
                     R = array(NA, dim), V = array(NA, dim), C = array(NA, dim))

  pars = parameter_setup(FOI_spillover, R0, vacc_data, pop_data, years_data, year0,
                         vaccine_efficacy, time_inc, mode_start, start_SEIRV, mode_time)
  pars$np_E_delay = nd1 * N_age
  pars$np_I_delay = nd2 * N_age
  pars$E_delay0 = rep(0, pars$np_E_delay)
  pars$I_delay0 = rep(0, pars$np_I_delay)
  for(div in 1:n_divs){
    n_particles = n_particles_list[div]

    x <- dust_system_create(SEIRVModelDelay, pars = pars, n_particles = n_particles, n_threads = n_threads, time = 0, dt = 1,
                            deterministic = FALSE, preserve_particle_dimension = TRUE)
    dust_system_set_state_initial(x)
    x_res <- dust_system_simulate(x, times = c(step_begin:step_end))

    if(div ==  1){
      n_p0 = 0
      index = dust_unpack_index(x)
    } else{
      n_p0 = sum(n_particles_list[c(1:(div - 1))])
    }

    n_p_values = c(1:n_particles) + n_p0
    dim = c(N_age, n_particles, t_pts_out)
    output_data$day = x_res[1, 1, ]
    output_data$year = x_res[2, 1, ]
    output_data$FOI_total[n_p_values, ] = array(x_res[3, , ]/time_inc, dim = c(n_particles, t_pts_out))
    output_data$S[, n_p_values, ] = array(x_res[index$S, , ], dim)
    output_data$E[, n_p_values, ] = array(x_res[index$E, , ], dim)
    output_data$I[, n_p_values, ] = array(x_res[index$I, , ], dim)
    output_data$R[, n_p_values, ] = array(x_res[index$R, , ], dim)
    output_data$V[, n_p_values, ] = array(x_res[index$V, , ], dim)
    output_data$C[, n_p_values, ] = array(x_res[index$C, , ], dim)
    x_res = NULL
    gc()
  }

  return(output_data)
}
#-------------------------------------------------------------------------------
#' @title Model_Run_VTrack
#'
#' @description Run SEIR model with separate vaccination track for single region
#'
#' @details Accepts epidemiological + population parameters and model settings; runs SEIR model with
#' separate vaccination track for one region over a specified time period for a number of particles/threads and outputs time-dependent SEIRV
#' values, infection numbers and/or total force of infection values.
#'
#' @param FOI_spillover Vector of values of force of infection due to spillover from sylvatic reservoir
#'   (size depends on mode_time)
#' @param R0 Vector of values of basic reproduction number for urban spread of infection (size depends on mode_time)
#' @param vacc_data Projected vaccination-based immunity (assuming vaccine_efficacy = 1) by age group and year
#' @param pop_data Population by age group and year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param year0 First year in population/vaccination data
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param time_inc Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination \cr
#'  If mode_start = 0, only vaccinated individuals \cr
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age) \cr
#'  If mode_start = 2, use SEIRV input in list from previous run(s) \cr
#' @param start_SEIRV SEIRV data from end of a previous run to use as input (if mode_start = 2)
#'   #TODO - Switch to SEIR multitrack data
#' @param mode_time Type of time dependence of FOI_spillover and R0 to be used: \cr
#'  If mode_time = 0, no time variation (constant values)\cr
#'  If mode_time = 1, FOI/R0 vary annually without seasonality (number of values = number of years to consider) \cr
#'  If mode_time = 2, FOI/R0 vary with monthly seasonality without inter - annual variation (number of values = 12) \cr
#'  If mode_time = 3, FOI/R0 vary with daily seasonality without inter - annual variation (number of values = 365/time_inc) \cr
#'  If mode_time = 4, FOI/R0 vary annually with monthly seasonality (number of values = 12*number of years to consider) \cr
#'  If mode_time = 5, FOI/R0 vary annually with daily seasonality (number of values = (365/time_inc)*number of years to consider)
#' @param n_particles number of particles to use
#' @param n_threads number of threads to use
#' @param deterministic TRUE/FALSE  -  set model to run in deterministic mode if TRUE
#' '
#' @export
#'
Model_Run_VTrack <- function(FOI_spillover = 0.0, R0 = 1.0, vacc_data = list(), pop_data = list(), years_data = c(1940:1941),
                      year0 = 1940, vaccine_efficacy = 1.0, time_inc = 1.0, mode_start = 0,
                      start_SEIRV = list(), mode_time = 0, n_particles = 1, n_threads = 1, deterministic = FALSE) {

  #TODO Add assert_that functions (NB  -  Some checks carried out in parameter_setup)
  assert_that(n_particles <= 20, msg = "Number of particles must be 20 or less")

  N_age = length(pop_data[1, ]) #Number of age groups
  step_begin = ((years_data[1] - year0)*(365/time_inc)) #Step at which data starts being saved for final output
  step_end = ((max(years_data) + 1 - year0)*(365/time_inc)) - 1 #Step at which to end
  t_pts_out = step_end - step_begin + 1 #Number of time points in final output data

  pars = parameter_setup(FOI_spillover, R0, vacc_data, pop_data, years_data, year0,
                         vaccine_efficacy, time_inc, mode_start, start_SEIRV, mode_time)
  pars2=pars
  pars2$S_0=pars2$E_0=pars2$I_0=pars2$R_0=array(0,dim=c(N_age,2))
  pars2$V_0=NULL
  pars2$S_0[,1]=pars$S_0
  pars2$S_0[,2]=pars$V_0 #TODO - split vaccinated between S and R
  pars2$E_0[,1]=pars$E_0
  pars2$I_0[,1]=pars$I_0
  pars2$R_0[,1]=pars$R_0

  x <- dust_system_create(SEIRModelVtrack, pars = pars2,
                          n_particles = n_particles, n_threads = n_threads, time = 0, dt = 1,
                          deterministic = deterministic, preserve_particle_dimension = TRUE)
  index = dust_unpack_index(x)
  dust_system_set_state_initial(x)
  x_res <- dust_system_simulate(x, times = c(step_begin:step_end))

  dim1 = c(N_age,n_particles, t_pts_out)
  dim2 = c(N_age, 2, n_particles, t_pts_out)
  output_data = list(day = x_res[1, 1, ], year = x_res[2, 1, ], FOI_total = x_res[3, , ]/time_inc,
                     S = array(x_res[index$S, , ], dim2), E = array(x_res[index$E, , ], dim2),
                     I = array(x_res[index$I, , ], dim2), R = array(x_res[index$R, , ], dim2),
                     C = array(x_res[index$C, , ], dim1))
  x_res = NULL
  gc()

  return(output_data)
}
#-------------------------------------------------------------------------------
#' @title Model_Versions_Run
#'
#' @description TBA
#'
#' @details TBA
#'
#' @param version Model version to use; choose from "Basic", "Delay", "Reactive" or "Split"
#' @param FOI_spillover Vector of values of force of infection due to spillover from sylvatic reservoir
#'   (size depends on mode_time)
#' @param R0 Vector of values of basic reproduction number for urban spread of infection (size depends on mode_time)
#' @param vacc_data Projected vaccination-based immunity (assuming vaccine_efficacy = 1) by age group and year
#' @param pop_data Population by age group and year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param year0 First year in population/vaccination data
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param time_inc Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination \cr
#'  If mode_start = 0, only vaccinated individuals \cr
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age) \cr
#'  If mode_start = 2, use SEIRV input in list from previous run(s) \cr
#' @param start_SEIRV SEIRV data from end of a previous run to use as input (if mode_start = 2)
#' @param mode_time Type of time dependence of FOI_spillover and R0 to be used: \cr
#'  If mode_time = 0, no time variation (constant values)\cr
#'  If mode_time = 1, FOI/R0 vary annually without seasonality (number of values = number of years to consider) \cr
#'  If mode_time = 2, FOI/R0 vary with monthly seasonality without inter - annual variation (number of values = 12) \cr
#'  If mode_time = 3, FOI/R0 vary with daily seasonality without inter - annual variation (number of values = 365/dt) \cr
#'  If mode_time = 4, FOI/R0 vary annually with monthly seasonality (number of values = 12 * number of years to consider) \cr
#'  If mode_time = 5, FOI/R0 vary annually with daily seasonality (number of values = (365/dt) * number of years to consider)
#' @param n_particles number of particles to use
#' @param n_threads number of threads to use
#' @param deterministic TRUE/FALSE  -  set model to run in deterministic mode if TRUE
#' @param response_delay (Reactive version only) Delay time in days between a threshold being reached and emergency conditions
#'   coming into effect
#' @param p_rep (Reactive version only) Probabilities of an infection being reported as a case before emergency conditions triggered
#'   (1st value) or after emergency conditions triggered (2nd value)
#' @param case_threshold (Reactive version only) Threshold total no. reported cases to trigger emergency conditions
#' @param cluster_threshold (Reactive version only) Threshold current infectious fraction to trigger emergency conditions
#' @param vacc_cov_cam (Reactive version only) Target vaccination coverage by age group during emergency campaign
#' @param t_cam (Reactive version only) Duration in days of emergency vaccination campaign
#' '
#' @export
#'
Model_Versions_Run <- function(version = "Basic", FOI_spillover = 0.0, R0 = 1.0, vacc_data = list(), pop_data = list(),
                               years_data = c(1940:1941), year0 = 1940, vaccine_efficacy = 1.0, time_inc = 1.0,
                               mode_start = 0, start_SEIRV = list(), mode_time = 0, n_particles = 1, n_threads = 1,
                               deterministic = FALSE, response_delay = 56.0, p_rep = c(0.0, 0.0), case_threshold = Inf,
                               cluster_threshold = Inf, vacc_cov_cam = c(), t_cam = 0){

  assert_that(version %in% c("Basic", "Delay", "Reactive", "Split"))

  if(version=="Basic"){
    output <- YEP::Model_Run(FOI_spillover, R0, vacc_data, pop_data, years_data, year0, vaccine_efficacy, time_inc,
                             output_type = "full", mode_start, start_SEIRV, mode_time, n_particles, n_threads, deterministic)
  }
  if(version=="Delay"){
    output <- Model_Run_Delay(FOI_spillover, R0, vacc_data, pop_data, years_data, year0, vaccine_efficacy, time_inc,
                              mode_start, start_SEIRV, mode_time, n_particles, n_threads, deterministic)
  }
  if(version=="Reactive"){
    output <- Model_Run_Delay_Reactive(FOI_spillover, R0, vacc_data, pop_data, years_data, year0, vaccine_efficacy, time_inc,
                                       mode_start, start_SEIRV, mode_time, n_particles, n_threads, deterministic,
                                       response_delay, p_rep, case_threshold, cluster_threshold, vacc_cov_cam, t_cam)
  }
  if(version=="Split"){
    output <- Model_Run_Split(FOI_spillover, R0, vacc_data, pop_data, years_data, year0, vaccine_efficacy, time_inc,
                              mode_start, start_SEIRV, mode_time, n_particles, n_threads, deterministic)
  }

  return(output)
}
