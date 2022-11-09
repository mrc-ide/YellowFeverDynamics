# R file for functions relating to annual case and death data in YellowFeverDynamics package
#-------------------------------------------------------------------------------
#' @title case_data_generate
#'
#' @description Take in single set of population data and model parameters, output infection data only
#'
#' @details Accepts epidemiological + population parameters and model settings; runs full version of SEIRV model
#' for one region over a specified time period for a number of particles/threads and outputs infection numbers at each
#' time increment only; optimized for running a large number of repetitions
#'
#' @param FOI_spillover = Force of infection due to spillover from sylvatic reservoir
#' @param R0 = Reproduction number for urban spread of infection
#' @param vacc_data Vaccination coverage in each age group by year
#' @param pop_data Population in each age group by year
#' @param year0 First year in population/vaccination data
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param n_reps Number of runs
#' @param year_end year to run up to
#' @param year_data_begin year to begin saving data
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' @param dt Time increment in days to use in model (should be either 1.0 or 5.0 days)
#' '
#' @export
#'
case_data_generate <- function(FOI_spillover=0.0,R0=1.0,vacc_data=list(),pop_data=list(),year0=1940,
                               mode_start=0,n_reps=1,year_end=2000,year_data_begin=1999,
                               vaccine_efficacy=vaccine_efficacy,start_SEIRV=list(),dt=1.0) {

  assert_that(n_reps>0)

  division=10
  n_particles0=min(division,n_reps)
  n_threads=min(10,n_particles0)
  n_divs=ceiling(n_reps/division)
  if(n_divs==1){
    n_particles_list=n_particles0
  } else {
    n_particles_list=c(rep(n_particles0,n_divs-1),n_reps-(division*(n_divs-1)))
  }

  pars=parameter_setup(FOI_spillover,R0,vacc_data,pop_data,year0,mode_start,year_end,
                       year_data_begin,vaccine_efficacy,start_SEIRV,dt)

  n=4 #Number of non-vector outputs
  N_age=length(pop_data[1,])
  t_pts=c(1:((year_end-year0)*(365/dt)))
  n_data_pts=(6*N_age)+n
  n_steps=length(t_pts)
  step0=(year_data_begin-year0)*(365/dt)
  steps=n_steps-step0
  results_data=list(year=sort(rep(c(year_data_begin:(year_end-1)),(365/dt))),C=array(data=rep(0,n_reps*steps),
                                                                                     dim=c(n_reps,steps)))
  pts_select=c(((5*N_age)+1+n):((6*N_age)+n))

  for(div in 1:n_divs){
    n_particles=n_particles_list[div]
    reps=c(1:n_particles)+((div-1)*division)
    x <- FullModelOD$new(pars=pars,time = 1,n_particles = n_particles,n_threads = n_threads)

    x_res <- array(NA, dim = c(n_data_pts, n_particles))
    for(t in step0:n_steps){
      x_res <- x$run(t)
      if(n_particles==1){
        results_data$C[reps,t-step0]=sum(x_res[pts_select])
      } else {
        results_data$C[reps,t-step0]=colSums(x_res[pts_select,],dims=1)
      }
    }
  }

  return(results_data)
}
#-------------------------------------------------------------------------------
# [TODO: CHANGE AND IF NECESSARY CREATE SEPARATE FUNCTION FOR OTHER DATA FORMATS]
#' @title deaths_compare
#'
#' @description Compare modelled and observed deaths using negative binomial
#'
#' @details Compares modelled data on fatal cases per year and compared with observed data, calculating logarithmic
#' likelihood of observing the latter given the former, using a negative binomial formula.
#'
#' @param model_data Modelled data in data frame format (list of years and number of reported deaths in each)
#' @param obs_data Data frame containing observed death data
#' '
#' @export
#'
deaths_compare <- function(model_data=list(),obs_data=list()){

  assert_that(is.null(model_data$rep_deaths)==FALSE)
  assert_that(is.null(obs_data$deaths)==FALSE)
  assert_that(length(model_data$rep_deaths)==length(obs_data$deaths))
  model_data$rep_deaths[model_data$rep_deaths==0]=0.1

  like_values=dnbinom(x=obs_data$deaths,mu=model_data$rep_deaths,size=rep(1,length(obs_data$deaths)),log=TRUE)
  LogLikelihood=sum(like_values,na.rm=TRUE)

  return(LogLikelihood)
}
#-------------------------------------------------------------------------------
# [TODO: CHANGE AND IF NECESSARY CREATE SEPARATE FUNCTION FOR OTHER DATA FORMATS]
#' @title cases_compare
#'
#' @description Compare modelled and observed severe cases using negative binomial
#'
#' @details Compares modelled data on severe cases per year and compared with observed data, calculating logarithmic
#' likelihood of observing the latter given the former, using a negative binomial formula.
#'
#' @param model_data Modelled data in data frame format (list of years and number of reported severe cases in each)
#' @param obs_data Data frame containing annual observed case data
#' '
#' @export
#'
cases_compare <- function(model_data=list(),obs_data=list()){

  assert_that(is.null(model_data$rep_cases)==FALSE)
  assert_that(is.null(obs_data$cases)==FALSE)
  assert_that(length(model_data$rep_cases)==length(obs_data$cases))
  model_data$rep_cases[model_data$rep_cases==0]=0.1

  like_values=dnbinom(x=obs_data$cases,mu=model_data$rep_cases,size=rep(1,length(obs_data$cases)),log=TRUE)
  LogLikelihood=sum(like_values,na.rm=TRUE)

  return(LogLikelihood)
}
