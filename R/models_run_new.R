#Test versions of model run functions for updated Odin files
#-------------------------------------------------------------------------------
#' @title Model_Run_Delay_New
#'
#' @description Run SEIRV model version using time delay instead of rate to move individuals from E-I and I-R
#'
#' @details Accepts epidemiological + population parameters and model settings; runs version of SEIRV model which
#' uses a fixed time delay instead of a rate to move individuals from exposed (E) to infectious (I) and from
#' infectious to recovered (R). The model is run for one region over a specified time period for a number of
#' particles/threads and outputs time-dependent SEIRV values, infection numbers and total force of infection values.
#'
#' @param FOI_spillover Matrix of values of force of infection due to spillover from sylvatic reservoir
#'   (size depends on mode_time)
#' @param R0 Matrix of values of basic reproduction number for urban spread of infection (size depends on mode_time)
#' @param vacc_data Projected vaccination-based immunity (assuming vaccine_efficacy = 1) by region. age group and year
#' @param pop_data Population by region, age group and year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param year0 First year in population/vaccination data
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param time_inc Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination \cr
#'  If mode_start = 0, only vaccinated individuals \cr
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age) \cr
#'  If mode_start = 2, use SEIRV input in list from previous run(s) (TBD) \cr
#' @param start_SEIRV SEIRV data from end of a previous run to use as input (if mode_start = 2)
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
#' @param seed Random seed (set to NULL if not to be used)
#'
#' @export
#'
Model_Run_Delay_New <- function(FOI_spillover = c(), R0 = c(), vacc_data = list(), pop_data = list(), years_data = c(1940:1941),
                            year0 = 1940, vaccine_efficacy = 1.0, time_inc = 1.0, mode_start = 1, start_SEIRV = list(),
                            mode_time = 0, n_particles = 1, n_threads = 1, deterministic = FALSE, seed = NULL) {

  #TODO Add assert_that functions (NB - Some checks carried out in parameter_setup_old)
  assert_that(n_particles <=  20, msg = "Number of particles must be 20 or less")

  n_regions = dim(pop_data)[[1]]
  N_age = dim(pop_data)[[3]] #Number of age groups
  nd1 <- (t_incubation + t_latent)/time_inc
  nd2 <- t_infectious/time_inc
  step_begin = ((years_data[1] - year0) * (365/time_inc)) #Step at which data starts being saved for final output
  step_end = ((max(years_data) + 1 - year0) * (365/time_inc)) - 1 #Step at which to end
  t_pts_out = step_end - step_begin + 1 #Number of time points in final output data

  pars = YEP::parameter_setup(FOI_spillover, R0, vacc_data, pop_data, years_data, year0,
                         vaccine_efficacy, time_inc, mode_start, start_SEIRV, mode_time)

  #Carrying forward delay from previous run may cause errors
  #if(mode_start == 2){pars$E_delay0 = start_SEIRV$E_delay} else {pars$E_delay0 = rep(0, nd1 * N_age)}
  pars$np_E_delay = nd1 * N_age
  pars$np_I_delay = nd2 * N_age
  pars$E_delay0 = array(0, dim=c(n_regions, pars$np_E_delay))
  pars$I_delay0 = array(0, dim=c(n_regions,pars$np_I_delay))

  x <- dust_system_create(SEIRVModelDelay_new, pars = pars,
                          n_particles = n_particles, n_threads = n_particles, time = 0,
                          dt = 1, seed = seed, deterministic = deterministic,
                          preserve_particle_dimension = TRUE)
  dust_system_set_state_initial(x)
  x_res <- dust_system_simulate(x, times = c(step_begin:step_end))
  index = dust_unpack_index(x)

  dim1 = c(n_regions, N_age, n_particles, t_pts_out)
  output_data = list(day = x_res[1, 1, ], year = x_res[2, 1, ],
                     FOI_total = array(x_res[3, , ]/time_inc, dim = c(n_regions,n_particles,t_pts_out)),
                     S = array(x_res[index$S, , ], dim1), E = array(x_res[index$E, , ], dim1),
                     #E_delay = array(x_res[index$E_delay, , ], dim = c(N_age, nd1, n_particles, t_pts_out)),
                     I = array(x_res[index$I, , ], dim1), R = array(x_res[index$R, , ], dim1),
                     #I_delay = array(x_res[index$I_delay, , ], dim = c(N_age, nd2, n_particles, t_pts_out)),
                     V = array(x_res[index$V, , ], dim1), C = array(x_res[index$C, , ], dim1))

  return(output_data)
}
#-------------------------------------------------------------------------------
#' @title Model_Run_Delay_Reactive_New
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
#' @param seed Random seed (set to NULL if not to be used)
#'
#' @export
#'
Model_Run_Delay_Reactive_New <- function(FOI_spillover = c(), R0 = c(), vacc_data = list(), pop_data = list(), years_data = c(1940:1941),
                                year0 = 1940, vaccine_efficacy = 1.0, time_inc = 1.0, mode_start = 1, start_SEIRV = list(),
                                mode_time = 0, n_particles = 1, n_threads = 1, deterministic = FALSE,
                                response_delay = 56.0, p_rep = c(0.0, 0.0), case_threshold = Inf,
                                cluster_threshold = Inf, vacc_cov_cam = list(), t_cam = 0, seed = NULL) {

  #TODO Add assert_that functions (NB - Some checks carried out in parameter_setup_old)
  assert_that(n_particles <=  20, msg = "Number of particles must be 20 or less")

  n_regions = dim(pop_data)[[1]]
  N_age = dim(pop_data)[[3]] #Number of age groups
  nd1 <- (t_incubation + t_latent)/time_inc
  nd2 <- t_infectious/time_inc
  step_begin = ((years_data[1] - year0) * (365/time_inc)) #Step at which data starts being saved for final output
  step_end = ((max(years_data) + 1 - year0) * (365/time_inc)) - 1 #Step at which to end
  t_pts_out = step_end - step_begin + 1 #Number of time points in final output data

  pars = YEP::parameter_setup(FOI_spillover, R0, vacc_data, pop_data, years_data, year0,
                              vaccine_efficacy, time_inc, mode_start, start_SEIRV, mode_time)

  #Carrying forward delay from previous run may cause errors
  #if(mode_start == 2){pars$E_delay0 = start_SEIRV$E_delay} else {pars$E_delay0 = rep(0, nd1 * N_age)}
  pars$np_E_delay = nd1 * N_age
  pars$np_I_delay = nd2 * N_age
  pars$E_delay0 = array(0, dim=c(n_regions, pars$np_E_delay))
  pars$I_delay0 = array(0, dim=c(n_regions,pars$np_I_delay))
  pars$vacc_cov_cam = vacc_cov_cam
  pars$t_cam = t_cam
  pars$response_delay = response_delay
  pars$p_rep = p_rep
  pars$case_threshold = case_threshold
  pars$cluster_threshold = cluster_threshold

  #Check that there is no overlap between emergency campaign and other vaccination
  #(Not yet possible to adjust vaccine rates on the fly to deal with overlap)
  #TBC
  # for(i in 1:N_age){
  #   if(vacc_cov_cam[i]>0){assert_that(all(pars2$vacc_rate_daily[i, ] == 0))}
  # }

  x <- dust_system_create(SEIRVModelDelayReactive_new, pars = pars,
                          n_particles = n_particles, n_threads = n_particles, time = 0,
                          dt = 1, seed = seed, deterministic = deterministic,
                          preserve_particle_dimension = TRUE)
  dust_system_set_state_initial(x)
  x_res <- dust_system_simulate(x, times = c(step_begin:step_end))
  index = dust_unpack_index(x)

  dim1 = c(n_regions, N_age, n_particles, t_pts_out)
  dim2 = c(n_regions,n_particles,t_pts_out)
  output_data = list(day = x_res[1, 1, ], year = x_res[2, 1, ],
                     FOI_total = array(x_res[3, , ]/time_inc, dim2),
                     S = array(x_res[index$S, , ], dim1), E = array(x_res[index$E, , ], dim1),
                     #E_delay = array(x_res[index$E_delay, , ], dim = c(N_age, nd1, n_particles, t_pts_out)),
                     I = array(x_res[index$I, , ], dim1), R = array(x_res[index$R, , ], dim1),
                     #I_delay = array(x_res[index$I_delay, , ], dim = c(N_age, nd2, n_particles, t_pts_out)),
                     V = array(x_res[index$V, , ], dim1),
                     C_rep_total = array(x_res[index$C_rep_total, , ], dim2),
                     flag_cam = array(x_res[index$flag4,,],dim2))

  return(output_data)
}
