#Test versions of model run functions for updated Odin files
#TODO - test
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

  N_age = length(pop_data[1, 1, ]) #Number of age groups
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
  pars$E_delay0 = rep(0, pars$np_E_delay)
  pars$I_delay0 = rep(0, pars$np_I_delay)

  x <- dust_system_create(SEIRVModelDelay, pars = pars,
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
                     E_delay = array(x_res[index$E_delay, , ], dim = c(N_age, nd1, n_particles, t_pts_out)),
                     I = array(x_res[index$I, , ], dim1), R = array(x_res[index$R, , ], dim1),
                     I_delay = array(x_res[index$I_delay, , ], dim = c(N_age, nd2, n_particles, t_pts_out)),
                     V = array(x_res[index$V, , ], dim1), C = array(x_res[index$C, , ], dim1))

  return(output_data)
}
