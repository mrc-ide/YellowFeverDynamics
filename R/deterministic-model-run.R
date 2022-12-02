#-------------------------------------------------------------------------------
#' @title Full_Model_Run_Deterministic
#'
#' @description Run deterministic version of SEIRV model
#'
#' @details Accepts epidemiological/population parameters and model settings; runs deterministic version of SEIRV model
#' for one region over a specified time period for a number of particles/threads and outputs time-dependent SEIRV
#' values, infection numbers and total force of infection values.
#'
#' @param FOI_spillover Force of infection due to spillover from sylvatic reservoir
#' @param R0 Reproduction number for urban spread of infection
#' @param vacc_data Vaccination coverage in each age group by year
#' @param pop_data Population in each age group by year
#' @param year0 First year in population/vaccination data
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRVC input in list from previous run(s)
#' @param n_particles number of particles to use
#' @param n_threads number of threads to use
#' @param year_end year to run up to
#' @param year_data_begin year to begin saving data
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' @param dt Time increment in days to use in model (should be either 1.0 or 5.0 days)
#' '
#' @export
#'
Full_Model_Run_Deterministic <- function(FOI_spillover=0.0,R0=1.0,vacc_data=list(),pop_data=list(),year0=1940,
                                         mode_start=0,n_particles=1,n_threads=1,year_end=2000,year_data_begin=1999,
                                         vaccine_efficacy=1.0,start_SEIRV=list(),dt=1.0) {

  assert_that(n_particles %in% c(1:20),msg="Number of particles must be an integer between 1 and 20")
  assert_that(n_threads<=n_particles,msg="Number of threads must be equal to or less than number of particles")
  assert_that(n_threads>0,msg="Number of threads must be between 1 and number of particles")

  x <- FullModelODDeterministic$new(pars=parameter_setup(FOI_spillover,R0,vacc_data,pop_data,year0,mode_start,year_end,
                                                         year_data_begin,vaccine_efficacy,start_SEIRV,dt),
                                    time = 1,n_particles = n_particles,n_threads = n_threads)

  n_nv=4 #Number of non-vector outputs
  N_age=length(pop_data[1,]) #Number of age groups
  t_pts_all=c(1:((year_end-year0)*(365/dt))) #All output time points
  n_data_pts=(6*N_age)+n_nv #Number of data values per time point in output
  n_steps=length(t_pts_all) #Total number of output time points
  step0=(year_data_begin-year0)*(365/dt) #Step at which data starts being saved for final output
  t_pts_out=n_steps-step0 #Number of time points in final output data
  x_res <- array(NA, dim = c(n_data_pts, n_particles, t_pts_out))
  for(t in step0:n_steps){
    x_res[,,t-step0] <- x$run(t)
  }

  return(list(day=array(x_res[2,,],dim=c(n_particles,t_pts_out)),
              year=array(x_res[3,,],dim=c(n_particles,t_pts_out)),
              FOI_total=array(x_res[4,,],dim=c(n_particles,t_pts_out)),
              S=array(x_res[c((1+n_nv):(N_age+n_nv)),,],dim=c(N_age,n_particles,t_pts_out)),
              E=array(x_res[c((N_age+1+n_nv):((2*N_age)+n_nv)),,],dim=c(N_age,n_particles,t_pts_out)),
              I=array(x_res[c(((2*N_age)+1+n_nv):((3*N_age)+n_nv)),,],dim=c(N_age,n_particles,t_pts_out)),
              R=array(x_res[c(((3*N_age)+1+n_nv):((4*N_age)+n_nv)),,],dim=c(N_age,n_particles,t_pts_out)),
              V=array(x_res[c(((4*N_age)+1+n_nv):((5*N_age)+n_nv)),,],dim=c(N_age,n_particles,t_pts_out)),
              C=array(x_res[c(((5*N_age)+1+n_nv):((6*N_age)+n_nv)),,],dim=c(N_age,n_particles,t_pts_out))))
}
#-------------------------------------------------------------------------------
#' @title case_data_generate_deterministic
#'
#' @description Take in single set of population data and model parameters, output infection data only (deterministic)
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
case_data_generate_deterministic <- function(FOI_spillover=0.0,R0=1.0,vacc_data=list(),pop_data=list(),year0=1940,
                               mode_start=0,n_reps=1,year_end=2000,year_data_begin=1999,
                               vaccine_efficacy=vaccine_efficacy,start_SEIRV=list(),dt=1.0) {

  assert_that(n_reps>0,msg="Number of runs must be positive")

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
    x <- FullModelODDeterministic$new(pars=pars,time = 1,n_particles = n_particles,n_threads = n_threads)

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
