# R file for general functions in YellowFeverDynamics package
#------------------------------------------------
#Global variables
p_severe_inf=0.12 #Probability that an infection is severe
p_death_severe_inf=0.39 #Probability that a severe infection becomes fatal
t_incubation <- 5 #Time for cases to incubate in mosquito
t_latent <- 5 #Latent period before cases become infectious
t_infectious <- 5 #Time cases remain infectious
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.
#' @useDynLib YellowFeverDynamics, .registration = TRUE
#' @importFrom assertthat assert_that
#' @import dde
#' @importFrom mvtnorm rmvnorm
#' @importFrom Rmisc CI
#' @importFrom stats cov dexp dnbinom dnorm prop.test rbinom runif
#' @importFrom tgp lhs
#' @importFrom truncdist dtrunc
#' @importFrom utils write.csv
#' @import YEP
#------------------------------------------------
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
#' @details Accepts epidemiological + population parameters and model settings; runs SEIRV delay model
#' for one region over a specified time period for a number of particles/threads and outputs time-dependent SEIRV
#' values, infection numbers and total force of infection values.
#'
#' @param FOI_spillover Force of infection due to spillover from sylvatic reservoir
#' @param R0 Basic reproduction number for urban spread of infection
#' @param vacc_data Vaccination coverage in each age group by year
#' @param pop_data Population in each age group by year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' @param year0 First year in population/vaccination data
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param dt Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param n_particles number of particles to use
#' @param n_threads number of threads to use
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' '
#' @export
#'
Model_Run_Delay <- function(FOI_spillover = 0.0,R0 = 1.0,vacc_data = list(),pop_data = list(),years_data = c(1940:1941),
                      start_SEIRV = list(), year0 = 1940, mode_start = 0,
                      vaccine_efficacy = 1.0, dt = 1.0, n_particles = 1, n_threads = 1, deterministic = FALSE) {

  #TODO Add assert_that functions

  n_nv=3 #Number of non-vector outputs
  N_age=length(pop_data[1,]) #Number of age groups
  nd1 <- (t_incubation+t_latent)/dt
  nd2 <- t_infectious/dt
  nd=nd1+nd2
  n_data_pts=((6+nd1+nd2)*N_age)+n_nv #Number of data values per time point in output
  step_begin=((years_data[1]-year0)*(365/dt)) #Step at which data starts being saved for final output
  step_end=((max(years_data)+1-year0)*(365/dt))-1 #Step at which to end
  t_pts_out=step_end-step_begin+1 #Number of time points in final output data

  x <- SEIRVModelDelay$new(pars=parameter_setup(FOI_spillover,R0,vacc_data,pop_data,year0,years_data,mode_start,
                                                 vaccine_efficacy,start_SEIRV,dt),
                            time = 0, n_particles = n_particles, n_threads = n_threads, deterministic = deterministic)

  x_res <- array(NA, dim = c(n_data_pts, n_particles, t_pts_out))
  for(step in step_begin:step_end){
    x_res[,,step-step_begin+1] <- x$run(step)
  }
  if(step_begin==0){x_res[2,,1]=rep(year0,n_particles)}

  dimensions=c(N_age,n_particles,t_pts_out)
  output_data=list(day=x_res[1,1,],year=x_res[2,1,])
  output_data$FOI_total=array(x_res[3,,]/dt,dim=c(n_particles,t_pts_out))
  output_data$S=array(x_res[c((1+n_nv):(N_age+n_nv)),,],dim=dimensions)
  output_data$E=array(x_res[c((N_age+1+n_nv):((2*N_age)+n_nv)),,],dim=dimensions)
  output_data$E_delay=array(x_res[c(((2*N_age)+1+n_nv):(((2+nd1)*N_age)+n_nv)),,],dim=c(N_age,nd1,n_particles,t_pts_out))
  output_data$I=array(x_res[c((((2+nd1)*N_age)+1+n_nv):(((3+nd1)*N_age)+n_nv)),,],dim=dimensions)
  output_data$I_delay=array(x_res[c((((2+nd1)*N_age)+1+n_nv):(((2+nd)*N_age)+n_nv)),,],dim=c(N_age,nd2,n_particles,t_pts_out))
  output_data$R=array(x_res[c((((3+nd)*N_age)+1+n_nv):(((4+nd)*N_age)+n_nv)),,],dim=dimensions)
  output_data$V=array(x_res[c((((4+nd)*N_age)+1+n_nv):(((5+nd)*N_age)+n_nv)),,],dim=dimensions)
  output_data$C=array(x_res[c((((5+nd)*N_age)+1+n_nv):(((6+nd)*N_age)+n_nv)),,],dim=dimensions)

  return(output_data)
}
#-------------------------------------------------------------------------------
#' @title Model_Run_Reactive
#'
#' @description Runs reactive version of SEIRV model
#'
#' @details Accepts epidemiological + population parameters and model settings; runs SEIRV model
#' for one region over a specified time period for a number of particles/threads and outputs time-dependent SEIRV
#' values, infection numbers and total force of infection values. Alternate version incorporating case reporting
#' and reactive surveillance/control measures based on case numbers
#'
#' @param FOI_spillover Force of infection due to spillover from sylvatic reservoir
#' @param R0 Basic reproduction number for urban spread of infection
#' @param vacc_data1 [TBA]
#' @param vacc_data2 [TBA]
#' @param pop_data Population in each age group by year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' @param year0 First year in population/vaccination data
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param dt Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param n_particles number of particles to use
#' @param n_threads number of threads to use
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' @param p_rep Probabilities of an infection being reported as a case under different conditions (TBA)
#' @param outbreak_threshold1 Threshold total no. reported cases to trigger outbreak flag 1
#' @param cluster_threshold1 Threshold current infectious fraction to trigger cluster flag 1
#' '
#' @export
#'
Model_Run_Reactive <- function(FOI_spillover = 0.0,R0 = 1.0,vacc_data1 = list(), vacc_data2 = list(),pop_data = list(),
                               years_data = c(1940:1941), start_SEIRV = list(), year0 = 1940, mode_start = 0,vaccine_efficacy = 1.0,
                               dt = 1.0, n_particles = 1, n_threads = 1, deterministic = FALSE, p_rep = 1.0, outbreak_threshold1 = 1,
                               cluster_threshold1 = 1.0) {

  #TODO Add assert_that functions

  n_nv=11 #Number of non-vector outputs
  N_age=length(pop_data[1,]) #Number of age groups
  n_data_pts=(7*N_age)+n_nv #Number of data values per time point in output
  step_begin=((years_data[1]-year0)*(365/dt)) #Step at which data starts being saved for final output
  step_end=((max(years_data)+1-year0)*(365/dt))-1 #Step at which to end
  t_pts_out=step_end-step_begin+1 #Number of time points in final output data

  pars1=parameter_setup(FOI_spillover,R0,vacc_data1,pop_data,year0,years_data,mode_start,
                       vaccine_efficacy,start_SEIRV,dt)
  n_years=length(pop_data[,1])-1
  inv_365=1.0/365.0
  vacc_rate_annual2=array(NA,dim=c(N_age,n_years,2))
  vacc_rate_annual2[,,1]=pars1$vacc_rate_annual
  for(n_year in 1:n_years){
    for(i in 1:N_age){
      if(i==1){
        vacc_rate_annual2[i,n_year,2]=vacc_data2[n_year+1,i]*inv_365
      } else {
        vacc_rate_annual2[i,n_year,2]=max(0.0,vacc_data2[n_year+1,i]-vacc_data2[n_year,i-1])*inv_365
      }
    }
  }
  pars2=list(FOI_spillover=pars1$FOI_spillover,R0=pars1$R0,vacc_rate_annual=vacc_rate_annual2,
             Cas0=pars1$Cas0,Exp0=pars1$Exp0,Inf0=pars1$Inf0,N_age=pars1$N_age,Rec0=pars1$Rec0,Sus0=pars1$Sus0,Vac0=pars1$Vac0,
             dP1_all=pars1$dP1_all,dP2_all=pars1$dP2_all,n_years=pars1$n_years,year0=pars1$year0,vaccine_efficacy=pars1$vaccine_efficacy,
             dt=pars1$dt,t_incubation=pars1$t_incubation,t_latent=pars1$t_latent,t_infectious=pars1$t_infectious,
             p_rep=p_rep,outbreak_threshold1=outbreak_threshold1,cluster_threshold1=cluster_threshold1)

  x <- SEIRVModelReactive$new(pars=pars2,time = 0, n_particles = n_particles, n_threads = n_threads, deterministic = deterministic)

  x_res <- array(NA, dim = c(n_data_pts, n_particles, t_pts_out))
  for(step in step_begin:step_end){
    x_res[,,step-step_begin+1] <- x$run(step)
  }
  if(step_begin==0){x_res[2,,1]=rep(year0,n_particles)}

  dims1=c(n_particles,t_pts_out)
  dims2=c(N_age,n_particles,t_pts_out)
  output_data=list(day=x_res[1,1,],year=x_res[2,1,],FOI_total=array(x_res[3,,]/dt,dim=dims1),C_rep_total=array(x_res[4,,],dim=dims1),
                   flag1a=array(x_res[5,,],dim=dims1),flag1b=array(x_res[6,,],dim=dims1),
                   flag2a=array(x_res[7,,],dim=dims1),flag2b=array(x_res[8,,],dim=dims1),
                   flag3=array(x_res[9,,],dim=dims1),report_rate=array(x_res[10,,],dim=dims1),
                   VR_check=array(x_res[11,,],dim=dims1))
  output_data$S=array(x_res[c((1+n_nv):(N_age+n_nv)),,],dim=dims2)
  output_data$E=array(x_res[c((N_age+1+n_nv):((2*N_age)+n_nv)),,],dim=dims2)
  output_data$I=array(x_res[c(((2*N_age)+1+n_nv):((3*N_age)+n_nv)),,],dim=dims2)
  output_data$R=array(x_res[c(((3*N_age)+1+n_nv):((4*N_age)+n_nv)),,],dim=dims2)
  output_data$V=array(x_res[c(((4*N_age)+1+n_nv):((5*N_age)+n_nv)),,],dim=dims2)
  output_data$C=array(x_res[c(((5*N_age)+1+n_nv):((6*N_age)+n_nv)),,],dim=dims2)
  output_data$C_rep=array(x_res[c(((6*N_age)+1+n_nv):((6*N_age)+n_nv)),,],dim=dims2)

  return(output_data)
}
#-------------------------------------------------------------------------------
#' @title Model_Run_Split
#'
#' @description Run full version of SEIRV model with daily infection output split into sylvatic and urban
#'
#' @details Accepts epidemiological + population parameters and model settings; runs full version of SEIRV model
#' for one region over a specified time period for a number of particles/threads and outputs time-dependent SEIRV
#' values, infection numbers and total force of infection values. Daily infection output split into sylvatic and urban.
#'
#' @param FOI_spillover Force of infection due to spillover from sylvatic reservoir
#' @param R0 Basic reproduction number for urban spread of infection
#' @param vacc_data Vaccination coverage in each age group by year
#' @param pop_data Population in each age group by year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' @param year0 First year in population/vaccination data
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param dt Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param n_particles number of particles to use
#' @param n_threads number of threads to use
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' '
#' @export
#'
Model_Run_Split <- function(FOI_spillover = 0.0,R0 = 1.0,vacc_data = list(),pop_data = list(),years_data = c(1940:1941),
                            start_SEIRV = list(), year0 = 1940, mode_start = 0,
                            vaccine_efficacy = 1.0, dt = 1.0, n_particles = 1, n_threads = 1, deterministic = FALSE) {

  #TODO Add assert_that functions

  n_nv=4 #Number of non-vector outputs
  N_age=length(pop_data[1,]) #Number of age groups
  n_data_pts=(9*N_age)+n_nv #Number of data values per time point in output
  step_begin=((years_data[1]-year0)*(365/dt)) #Step at which data starts being saved for final output
  step_end=((max(years_data)+1-year0)*(365/dt))-1 #Step at which to end
  t_pts_out=step_end-step_begin+1 #Number of time points in final output data

  x <- SEIRVModelSplitInfection$new(pars=parameter_setup(FOI_spillover,R0,vacc_data,pop_data,year0,years_data,mode_start,
                                                 vaccine_efficacy,start_SEIRV,dt),
                            time = 0, n_particles = n_particles, n_threads = n_threads, deterministic = deterministic)

  x_res <- array(NA, dim = c(n_data_pts, n_particles, t_pts_out))
  for(step in step_begin:step_end){
    x_res[,,step-step_begin+1] <- x$run(step)
  }
  if(step_begin==0){x_res[2,,1]=rep(year0,n_particles)}

  dimensions=c(N_age,n_particles,t_pts_out)
  output_data=list(day=x_res[1,1,],year=x_res[2,1,],FOI_sylvatic=array(x_res[3,,]/dt,dim=c(n_particles,t_pts_out)),
                   FOI_urban=array(x_res[4,,]/dt,dim=c(n_particles,t_pts_out)))
  output_data$S=array(x_res[c((1+n_nv):(N_age+n_nv)),,],dim=dimensions)
  output_data$E_sylvatic=array(x_res[c((N_age+1+n_nv):((2*N_age)+n_nv)),,],dim=dimensions)
  output_data$E_urban = array(x_res[c(((2*N_age)+1+n_nv):((3*N_age)+n_nv)),,],dim=dimensions)
  output_data$I_sylvatic=array(x_res[c(((3*N_age)+1+n_nv):((4*N_age)+n_nv)),,],dim=dimensions)
  output_data$I_urban = array(x_res[c(((4*N_age)+1+n_nv):((5*N_age)+n_nv)),,],dim=dimensions)
  output_data$R=array(x_res[c(((5*N_age)+1+n_nv):((6*N_age)+n_nv)),,],dim=dimensions)
  output_data$V=array(x_res[c(((6*N_age)+1+n_nv):((7*N_age)+n_nv)),,],dim=dimensions)
  output_data$C_sylvatic=array(x_res[c(((7*N_age)+1+n_nv):((8*N_age)+n_nv)),,],dim=dimensions)
  output_data$C_urban = array(x_res[c(((8*N_age)+1+n_nv):((9*N_age)+n_nv)),,],dim=dimensions)

  return(output_data)
}

#-------------------------------------------------------------------------------
#' @title Model_Run_VarFR
#'
#' @description Run alternate version of SEIRV model with annually varying FOI_spillover and R0
#'
#' @details Accepts epidemiological + population parameters and model settings; runs full version of SEIRV model
#' for one region over a specified time period for a number of particles/threads and outputs time-dependent SEIRV
#' values, infection numbers and total force of infection values.
#'
#' @param FOI_spillover Vector of annual values of force of infection due to spillover from sylvatic reservoir
#' @param R0 Vector of annual values of basic reproduction number for urban spread of infection
#' @param vacc_data Vaccination coverage in each age group by year
#' @param pop_data Population in each age group by year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' @param year0 First year in population/vaccination data
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param dt Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param n_particles number of particles to use
#' @param n_threads number of threads to use
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' '
#' @export
#'
Model_Run_VarFR <- function(FOI_spillover = c(),R0 = c(),vacc_data = list(),pop_data = list(),years_data = c(1940:1941),
                            start_SEIRV = list(), year0 = 1940, mode_start = 0,
                            vaccine_efficacy = 1.0, dt = 1.0, n_particles = 1, n_threads = 1, deterministic = FALSE) {

  #TODO Add assert_that functions

  n_nv=3 #Number of non-vector outputs
  N_age=length(pop_data[1,]) #Number of age groups
  n_data_pts=(6*N_age)+n_nv #Number of data values per time point in output
  step_begin=((years_data[1]-year0)*(365/dt)) #Step at which data starts being saved for final output
  step_end=((max(years_data)+1-year0)*(365/dt))-1 #Step at which to end
  t_pts_out=step_end-step_begin+1 #Number of time points in final output data

  pars1=parameter_setup(FOI_spillover[1],R0[1],vacc_data,pop_data,year0,years_data,mode_start,
                        vaccine_efficacy,start_SEIRV,dt)
  pars2=list(FOI_spillover=FOI_spillover,R0=R0,vacc_rate_annual=pars1$vacc_rate_annual,
             Cas0=pars1$Cas0,Exp0=pars1$Exp0,Inf0=pars1$Inf0,N_age=pars1$N_age,Rec0=pars1$Rec0,Sus0=pars1$Sus0,Vac0=pars1$Vac0,
             dP1_all=pars1$dP1_all,dP2_all=pars1$dP2_all,n_years=pars1$n_years,year0=pars1$year0,vaccine_efficacy=pars1$vaccine_efficacy,
             dt=pars1$dt,t_incubation=pars1$t_incubation,t_latent=pars1$t_latent,t_infectious=pars1$t_infectious)

  x <- SEIRVModelVarFR$new(pars=pars2,time = 0, n_particles = n_particles, n_threads = n_threads, deterministic = deterministic)

  x_res <- array(NA, dim = c(n_data_pts, n_particles, t_pts_out))
  for(step in step_begin:step_end){
    x_res[,,step-step_begin+1] <- x$run(step)
  }
  if(step_begin==0){x_res[2,,1]=rep(year0,n_particles)}

  dimensions=c(N_age,n_particles,t_pts_out)
  output_data=list(day=x_res[1,1,],year=x_res[2,1,])
  output_data$FOI_total=array(x_res[3,,]/dt,dim=c(n_particles,t_pts_out))
  output_data$S=array(x_res[c((1+n_nv):(N_age+n_nv)),,],dim=dimensions)
  output_data$E=array(x_res[c((N_age+1+n_nv):((2*N_age)+n_nv)),,],dim=dimensions)
  output_data$I=array(x_res[c(((2*N_age)+1+n_nv):((3*N_age)+n_nv)),,],dim=dimensions)
  output_data$R=array(x_res[c(((3*N_age)+1+n_nv):((4*N_age)+n_nv)),,],dim=dimensions)
  output_data$V=array(x_res[c(((4*N_age)+1+n_nv):((5*N_age)+n_nv)),,],dim=dimensions)
  output_data$C=array(x_res[c(((5*N_age)+1+n_nv):((6*N_age)+n_nv)),,],dim=dimensions)

  return(output_data)
}
#-------------------------------------------------------------------------------
#' @title total_burden_estimate
#'
#' @description Function to calculate annual yellow fever burden across multiple regions based on derived parameters
#'
#' @details Function to take in parameter sets derived from MCMC fitting and use to calculate annual total and reported
#' case and death numbers for multiple regions to compare with external data
#'
#' @param type Type of parameter set (FOI only, FOI+R0, FOI and/or R0 coefficients associated with environmental
#'   covariates); choose from "FOI","FOI+R0","FOI enviro","FOI+R0 enviro"
#' @param param_dist Data frame of values of input parameters, one set per row
#' @param input_data List of population and vaccination data for multiple regions
#' @param start_SEIRV SEIRV data to use as input
#' @param years_data Vector of years for which to output data
#' @param n_reps Number of repeats over which to average results
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param flag_reporting Flag indicating whether to output number of reported severe and fatal cases
#' @param dt Time increment in days to use in model (should be either 1.0 or 5.0 days)
#' @param enviro_data enviro_data Data frame containing values of environmental covariates; set to NULL if not in use
#' @param R0_fixed_values Values of R0 to use if not being taken from parameter distribution
#' @param vaccine_efficacy Vaccine efficacy (set to NULL if being varied as a parameter)
#' @param p_rep_severe Probability of observation of severe infection (set to NULL if being varied as a parameter)
#' @param p_rep_death Probability of observation of death (set to NULL if being varied as a parameter)
#' @param m_FOI_Brazil Multiplier of spillover FOI for Brazil regions (set to NULL if being varied as a parameter)
#'
#' @export
#'
total_burden_estimate <- function(type="FOI+R0 enviro",param_dist=list(),input_data=list(),start_SEIRV=NULL,
                                  years_data=c(),n_reps=1,mode_start=1,flag_reporting=FALSE,dt=5.0,
                                  enviro_data=NULL,R0_fixed_values=NULL,vaccine_efficacy=NULL,
                                  p_rep_severe=NULL,p_rep_death=NULL,m_FOI_Brazil=1.0){

  assert_that(input_data_check(input_data),
              msg="Input data must be in standard format (see https://mrc-ide.github.io/YellowFeverDynamics/articles/CGuideAInputs.html )")
  assert_that(all(input_data$region_labels==enviro_data$region)==TRUE) #TODO - msg
  assert_that(min(years_data)>=input_data$years_labels[1]) #TODO - msg
  assert_that(type %in% c("FOI+R0","FOI","FOI+R0 enviro","FOI enviro"))
  assert_that(is.logical(flag_reporting))
  assert_that(all(param_dist>0.0),msg="All parameter values in distribution must be positive")

  n_param_sets=nrow(param_dist)
  n_years=length(years_data)
  year_data_begin=years_data[1]
  year_end=max(years_data)+1
  regions=input_data$region_labels
  n_regions=length(regions)
  case_ar1=death_ar1=array(NA,dim=c(n_years,n_regions,n_param_sets,n_reps))
  case_ar2=death_ar2=array(NA,dim=c(n_years,n_regions,n_param_sets))
  case_ar3=death_ar3=array(NA,dim=c(n_years,n_param_sets))
  FOI_values=R0_values=rep(NA,n_regions)
  if(type %in% c("FOI+R0 enviro","FOI enviro")){n_env_vars=ncol(enviro_data)-1}

  cat("\nSets:\n")
  for(n_param_set in 1:n_param_sets){
    cat(" ",n_param_set,sep="")
    if(n_param_set %% 10 == 0){cat("\n")}
    params=param_dist[n_param_set,]
    names(params)=colnames(param_dist)
    if(is.null(vaccine_efficacy)){vaccine_efficacy_set=params$vaccine_efficacy} else {vaccine_efficacy_set=vaccine_efficacy}
    if(is.null(p_rep_severe)){p_rep_severe_set=params$p_rep_severe} else {p_rep_severe_set=p_rep_severe}
    if(is.null(p_rep_death)){p_rep_death_set=params$p_rep_death} else {p_rep_death_set=p_rep_death}
    if(is.null(m_FOI_Brazil)){m_FOI_Brazil_set=params$m_FOI_Brazil} else {m_FOI_Brazil_set=m_FOI_Brazil}

    if(type %in% c("FOI+R0 enviro","FOI enviro")){
      if(type=="FOI+R0 enviro"){enviro_coeffs=params[c(1:(2*n_env_vars))]} else {enviro_coeffs=params[c(1:n_env_vars)]}
      for(n_region in 1:n_regions){
        model_params=param_calc_enviro(enviro_coeffs,
                                       as.numeric(enviro_data[enviro_data$region==regions[n_region],1+c(1:n_env_vars)]))
        FOI_values[n_region]=model_params$FOI
        if(substr(regions[n_region],1,3)=="BRA"){FOI_values[n_region]=FOI_values[n_region]*m_FOI_Brazil_set}
        if(type=="FOI+R0 enviro"){R0_values[n_region]=model_params$R0} else {
          R0_values[n_region]=R0_fixed_values[n_region]}
      }
    }
    if(type %in% c("FOI+R0","FOI")){
      FOI_values=as.numeric(params)[c(1:n_regions)]
      if(type=="FOI+R0"){R0_values=as.numeric(params[c((n_regions+1):(2*n_regions))])
      } else {R0_values=R0_fixed_values}
    }

    for(n_region in 1:n_regions){

      if(mode_start==2){ #TODO - Edit
        start_SEIRV_set=list(S=start_SEIRV$S[,n_region,n_param_set],
                         E=start_SEIRV$E[,n_region,n_param_set],I=start_SEIRV$I[,n_region,n_param_set],
                         R=start_SEIRV$R[,n_region,n_param_set],V=start_SEIRV$V[,n_region,n_param_set])
      } else {start_SEIRV_set=NULL}
      case_data <- YEP::Model_Run(FOI_values[n_region],R0_values[n_region],
                                  vacc_data=input_data$vacc_data[n_region,,],pop_data=input_data$pop_data[n_region,,],
                                  years_data=c(year_data_begin:year_end),start_SEIRV=start_SEIRV_set,output_type="case",
                                  year0=input_data$years_labels[1],mode_start,vaccine_efficacy,dt,n_reps,n_reps,FALSE)
      for(n_year in 1:n_years){
        for(rep in 1:n_reps){
          infs=floor(sum(case_data$C[rep,case_data$year==years_data[n_year]]))
          severe_infs=rbinom(1,infs,p_severe_inf)
          deaths=rbinom(1,severe_infs,p_death_severe_inf)
          case_ar1[n_year,n_region,n_param_set,rep]=severe_infs
          death_ar1[n_year,n_region,n_param_set,rep]=deaths
        }
        case_ar2[n_year,n_region,n_param_set]=sum(case_ar1[n_year,n_region,n_param_set,])/n_reps
        death_ar2[n_year,n_region,n_param_set]=sum(death_ar1[n_year,n_region,n_param_set,])/n_reps
      }
    }

    for(n_year in 1:n_years){
      case_ar3[n_year,n_param_set]=sum(case_ar2[n_year,,n_param_set])
      death_ar3[n_year,n_param_set]=sum(death_ar2[n_year,,n_param_set])
    }
  }

  if(flag_reporting){
    obs_case_ar2=obs_death_ar2=array(NA,dim=c(n_years,n_regions,n_param_sets))
    obs_case_ar3=obs_death_ar3=array(NA,dim=c(n_years,n_param_sets))

    for(n_param_set in 1:n_param_sets){
      for(n_region in 1:n_regions){
        for(n_year in 1:n_years){
          cases=case_ar2[n_year,n_region,n_param_set]
          deaths=death_ar2[n_year,n_region,n_param_set]
          obs_deaths=rbinom(1,floor(deaths),p_rep_death_set)
          obs_cases=obs_deaths+rbinom(1,floor(cases-deaths),p_rep_severe_set)
          obs_case_ar2[n_year,n_region,n_param_set]=obs_cases
          obs_death_ar2[n_year,n_region,n_param_set]=obs_deaths
        }
      }
      for(n_year in 1:n_years){
        obs_case_ar3[n_year,n_param_set]=sum(obs_case_ar2[n_year,,n_param_set])
        obs_death_ar3[n_year,n_param_set]=sum(obs_death_ar2[n_year,,n_param_set])
      }
    }
    plot_frame1=data.frame(year=rep(years_data,n_regions*n_param_sets),
                           region=rep(sort(rep(regions,n_years)),n_param_sets),
                           set=sort(rep(c(1:n_param_sets),n_years*n_regions)),
                           cases=as.vector(case_ar2),deaths=as.vector(death_ar2),
                           obs_cases=as.vector(obs_case_ar2),obs_deaths=as.vector(obs_death_ar2))
    plot_frame2=data.frame(year=rep(years_data,n_param_sets),set=sort(rep(c(1:n_param_sets),n_years)),
                           cases=as.vector(case_ar3),deaths=as.vector(death_ar3),
                           obs_cases=as.vector(obs_case_ar3),obs_deaths=as.vector(obs_death_ar3))

  } else {
    plot_frame1=data.frame(year=rep(years_data,n_regions*n_param_sets),
                           region=rep(sort(rep(regions,n_years)),n_param_sets),
                           set=sort(rep(c(1:n_param_sets),n_years*n_regions)),
                           cases=as.vector(case_ar2),deaths=as.vector(death_ar2))
    plot_frame2=data.frame(year=rep(years_data,n_param_sets),set=sort(rep(c(1:n_param_sets),n_years)),
                           cases=as.vector(case_ar3),deaths=as.vector(death_ar3))
  }


  return(list(by_region=plot_frame1,all=plot_frame2))
}
