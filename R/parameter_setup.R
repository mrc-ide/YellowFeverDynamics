#Parameter setup functions copied from YEP version 0.2
#-------------------------------------------------------------------------------
#' @title Parameter setup
#'
#' @description Set up parameters to input into model
#'
#' @details Takes in multiple inputs, outputs list for use by odin SEIRV model.
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
#' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' @param mode_time Type of time dependence of FOI_spillover and R0 to be used: \cr
#'  If mode_time = 0, no time variation (constant values)\cr
#'  If mode_time = 1, FOI/R0 vary annually without seasonality (number of values = number of years to consider) \cr
#'  If mode_time = 2, FOI/R0 vary with monthly seasonality without inter - annual variation (number of values = 12) \cr
#'  If mode_time = 3, FOI/R0 vary with daily seasonality without inter - annual variation (number of values = 365/time_inc) \cr
#'  If mode_time = 4, FOI/R0 vary annually with monthly seasonality (number of values = 12*number of years to consider) \cr
#'  If mode_time = 5, FOI/R0 vary annually with daily seasonality (number of values = (365/time_inc)*number of years to consider)
#' '
#' @export
#'
parameter_setup_old <- function(FOI_spillover = c(), R0 = c(), vacc_data = list(), pop_data = list(), years_data = c(), year0 = 1940,
                            vaccine_efficacy = 1.0, time_inc = 1.0, mode_start = 0, start_SEIRV = list(), mode_time = 0){

  assert_that(mode_start %in% c(0, 1, 2), msg = "mode_start must have value 0, 1, or 2")
  #if(mode_start == 3){mode_start = 1} #Temporary fix until mode_start harmonized across all functions/examples
  if(mode_start == 2){assert_that(is.null(start_SEIRV$S) == FALSE,
                                  msg = "When mode_start = 2, start_SEIRV data required")}
  assert_that(mode_time %in% c(0:5), msg = "mode_time must be an integer between 0 and 5")
  assert_that(all(FOI_spillover >= 0.0))
  assert_that(all(R0 >= 0.0))
  assert_that(length(pop_data[, 1]) > 1, msg = "Need population data for multiple years")
  assert_that(length(pop_data[1, ]) > 1, msg = "Need population data for multiple age groups")
  n_years = length(pop_data[, 1]) - 1
  N_age = length(pop_data[1, ])
  assert_that(length(vacc_data[, 1]) == n_years + 1,
              msg = "Population and vaccination data must be for same time periods")
  assert_that(length(vacc_data[1, ]) == N_age, msg = "No. age groups in population and vaccination data must match")
  assert_that(between(vaccine_efficacy,0.0,1.0),msg = "Vaccine efficacy must be between 0-1")
  assert_that(years_data[1] >= year0, msg = "First data year must be greater than or equal to year0")
  assert_that(max(years_data) + 1 - year0 <= n_years, msg = "Period of years_data must lie within population data")
  assert_that(time_inc %in% c(1, 2.5, 5), msg = "time_inc must have value 1, 2.5 or 5 days")
  pts_year = 365.0/time_inc
  n_t_pts = n_years*pts_year
  n_req = switch(mode_time + 1, 1, n_years, 12, pts_year, n_years*12, n_t_pts)
  assert_that(length(FOI_spillover) == n_req && length(R0) == n_req,
              msg = "Spillover FOI and R0 must be correct length for mode_time")
  inv_365 = 1.0/365.0

  date_values = switch(mode_time + 1,
                       rep(1, n_t_pts),
                       sort(rep(c(1:n_years),pts_year)),
                       1 + (floor(12*time_inc*inv_365*c(0:(n_t_pts - 1))) %% 12),
                       1 + (floor(time_inc*c(0:(n_t_pts - 1))) %% pts_year),
                       1 + (floor(12*time_inc*inv_365*c(0:(n_t_pts - 1))) %% 12) + (12*sort(rep(c(1:n_years) - 1,
                                                                                                pts_year))),
                       c(1:n_t_pts))

  FOI_spillover_t = FOI_spillover[date_values]
  R0_t = R0[date_values]

  P0 = S_0 = E_0 = I_0 = R_0 = V_0 = rep(0, N_age)
  dP1_all = dP2_all = vacc_rates = array(NA, dim = c(N_age, n_years))
  for(i in 1:N_age){ P0[i] = max(1.0, pop_data[1, i]) } #Set all population values to nonzero minimum to avoid NaN values
  for(n_year in 1:n_years){
    for(i in 1:N_age){
      dP1_all[i, n_year] = max(1.0, pop_data[n_year + 1, i])*inv_365
      dP2_all[i, n_year] = max(1.0, pop_data[n_year, i])*inv_365
      if(i == 1){
        vacc_rates[i, n_year] = vacc_data[n_year + 1, i]*inv_365
      } else {
        vacc_rates[i, n_year] = max(0.0, vacc_data[n_year + 1, i] - vacc_data[n_year, i - 1])*inv_365
      }
    }
  }

  vacc_initial = vacc_data[1, ]
  if(mode_start == 2){ #Use supplied SEIRV data
    S_0 = start_SEIRV$S
    E_0 = start_SEIRV$E
    I_0 = start_SEIRV$I
    R_0 = start_SEIRV$R
    V_0 = start_SEIRV$V
  } else {
    V_0 = P0*vacc_initial
    if(mode_start == 0){ #No initial immunity
      S_0 = P0*(1.0 - vacc_initial)
    } else { #Stratified herd immunity profile based on notional FOI (averaged over first year)
      R0_year0=mean(R0_t[c(1:pts_year)])
      FOI_spillover_year0=mean(FOI_spillover_t[c(1:pts_year)])
      ages = c(1:N_age) - 1
      if(R0_year0 <= 1.0){
        FOI_estimate = FOI_spillover_year0*365.0
      } else {
        estimation_results = nlm(imm_fraction_function_yfd, p =  - 4, R0_year0, ages, P0/sum(P0))
        FOI_estimate = min(0.1, (FOI_spillover_year0*365.0) + exp(estimation_results$estimate))
      }
      herd_immunity = 1.0 - (exp( - FOI_estimate*(ages + 0.5)))

      for(i in 1:N_age){
        if(vacc_initial[i]<herd_immunity[i]){
          R_0[i] = P0[i]*(herd_immunity[i] - vacc_initial[i])
          S_0[i] = P0[i]*(1.0 - herd_immunity[i])
        } else {
          S_0[i] = P0[i]*(1.0 - vacc_initial[i])
        }
      }
    }
  }

  return(list(FOI_spillover = FOI_spillover_t, R0 = R0_t, vacc_rate_daily = vacc_rates, N_age = N_age,
              S_0 = S_0, E_0 = E_0, I_0 = I_0, R_0 = R_0, V_0 = V_0, dP1_all = dP1_all, dP2_all = dP2_all,
              n_years = n_years, year0 = year0, vaccine_efficacy = vaccine_efficacy, time_inc = time_inc,
              t_incubation = t_incubation, t_latent = t_latent, t_infectious = t_infectious, n_t_pts = n_t_pts))
}
#-------------------------------------------------------------------------------
#' @title imm_fraction_function_yfd
#'
#' @description Function to estimate notional FOI for herd immunity based on R0 and population age distribution
#'
#' @details [TBA]
#'
#' @param log_lambda Natural logarithm of force of infection
#' @param R0 Basic reproduction number
#' @param ages List of age values
#' @param pop_fraction Population of each age group as proportion of total
#' '
#' @export
#'
imm_fraction_function_yfd <- function(log_lambda =  - 4, R0 = 1.0, ages = c(0:100), pop_fraction = rep(1/101, 101)){
  #TODO  -  Add assert_that functions

  lambda = exp(log_lambda)
  immunity = 1.0 - (exp( - lambda*(ages + 0.5)))
  imm_mean = sum(immunity*pop_fraction)
  imm_mean_target = 1.0 - (1.0/R0)
  dev = abs(imm_mean_target - imm_mean)

  return(dev)
}
