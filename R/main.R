# R file for general functions in YellowFeverDynamics package
#------------------------------------------------
#Global variables
t_incubation <- 5 #Time for cases to incubate in mosquito
t_latent <- 5 #Latent period before cases become infectious
t_infectious <- 5 #Time cases remain infectious
#------------------------------------------------
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
#' @details Accepts epidemiological + population parameters and model settings; runs version of SEIRV model which
#' uses a fixed time delay instead of a rate to move individuals from exposed (E) to infectious (I) and from
#' infectious to recovered (R). The model is run for one region over a specified time period for a number of
#' particles/threads and outputs time-dependent SEIRV values, infection numbers and total force of infection values.
#'
#' @param FOI_spillover Force of infection due to spillover from sylvatic reservoir
#' @param R0 Basic reproduction number for urban spread of infection
#' @param vacc_data Projected vaccination-based immunity (assuming vaccine_efficacy=1) by age group and year
#' @param pop_data Population by age group and year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param start_SEIRV SEIRV data (including E_delay and I_delay) from end of a previous run to use as input
#' @param year0 First year in population/vaccination data
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param time_inc Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param n_particles number of particles to use
#' @param n_threads number of threads to use
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' '
#' @export
#'
Model_Run_Delay <- function(FOI_spillover = 0.0, R0 = 1.0, vacc_data = list(), pop_data = list(), years_data = c(1940:1941),
                      start_SEIRV = list(), year0 = 1940, mode_start = 0, vaccine_efficacy = 1.0,
                      time_inc = 1.0, n_particles = 1, n_threads = 1, deterministic = FALSE) {

  #TODO Add assert_that functions (NB - Some checks carried out in parameter_setup)
  assert_that(n_particles <= 20, msg = "Number of particles must be 20 or less")

  N_age=length(pop_data[1,]) #Number of age groups
  nd1 <- (t_incubation+t_latent)/time_inc
  nd2 <- t_infectious/time_inc
  nd=nd1+nd2
  step_begin=((years_data[1]-year0)*(365/time_inc)) #Step at which data starts being saved for final output
  step_end=((max(years_data)+1-year0)*(365/time_inc))-1 #Step at which to end
  t_pts_out=step_end-step_begin+1 #Number of time points in final output data

  pars=parameter_setup(FOI_spillover,R0,vacc_data,pop_data,year0,years_data,mode_start,vaccine_efficacy,start_SEIRV,time_inc)

  #Carrying forward delay from previous run may cause errors
  #if(mode_start==2){pars$E_delay0=start_SEIRV$E_delay} else {pars$E_delay0=rep(0,nd1*N_age)}
  pars$np_E_delay=nd1*N_age
  pars$np_I_delay=nd2*N_age
  pars$E_delay0=rep(0,pars$np_E_delay)
  pars$I_delay0=rep(0,pars$np_I_delay)

  x <- dust_system_create(SEIRVModelDelay, pars = pars,
                          n_particles = n_particles, n_threads = n_threads, time = 0, dt = 1,
                          deterministic = deterministic, preserve_particle_dimension = TRUE)
  index=dust_unpack_index(x)
  dust_system_set_state_initial(x)
  x_res <- dust_system_simulate(x, c(step_begin:step_end))

  dimensions=c(N_age,n_particles,t_pts_out)
  output_data=list(day=x_res[1,1,],year=x_res[2,1,])
  output_data$FOI_total=array(x_res[3,,]/time_inc,dim=c(n_particles,t_pts_out))
  output_data$S=array(x_res[index$S,,],dim=dimensions)
  output_data$E=array(x_res[index$E,,],dim=dimensions)
  output_data$E_delay=array(x_res[index$E_delay,,],dim=c(N_age,nd1,n_particles,t_pts_out))
  output_data$I=array(x_res[index$I,,],dim=dimensions)
  output_data$I_delay=array(x_res[index$I_delay,,],dim=c(N_age,nd2,n_particles,t_pts_out))
  output_data$R=array(x_res[index$R,,],dim=dimensions)
  output_data$V=array(x_res[index$V,,],dim=dimensions)
  output_data$C=array(x_res[index$C,,],dim=dimensions)

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
#' @param FOI_spillover Force of infection due to spillover from sylvatic reservoir
#' @param R0 Basic reproduction number for urban spread of infection
#' @param vacc_data Projected vaccination-based immunity (assuming vaccine_efficacy=1 and in absence of emergency campaign)
#'   by age group and year
#' @param pop_data Population by age group and year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param start_SEIRV SEIRV data (including E_delay and I_delay) from end of a previous run to use as input
#' @param year0 First year in population/vaccination data
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param time_inc Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
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
Model_Run_Delay_Reactive <- function(FOI_spillover = 0.0,R0 = 1.0,vacc_data = list(), pop_data = list(),
                                     years_data = c(1940:1941), start_SEIRV = list(), year0 = 1940, mode_start = 0,
                                     vaccine_efficacy = 1.0, time_inc = 1.0, n_particles = 1, n_threads = 1, deterministic = FALSE,
                                     response_delay = 56.0, p_rep = c(0.0,0.0), case_threshold = Inf,
                                     cluster_threshold = Inf, vacc_cov_cam = c(), t_cam = 0) {

  #TODO Add assert_that functions (NB - Some checks carried out in parameter_setup)
  assert_that(n_particles <= 20, msg = "Number of particles must be 20 or less")

  N_age=length(pop_data[1,]) #Number of age groups
  assert_that(length(vacc_cov_cam)==N_age)
  nd1 <- (t_incubation+t_latent)/time_inc
  nd2 <- t_infectious/time_inc
  nd <- nd1 + nd2
  step_begin=((years_data[1]-year0)*(365/time_inc)) #Step at which data starts being saved for final output
  step_end=((max(years_data)+1-year0)*(365/time_inc))-1 #Step at which to end
  t_pts_out=step_end-step_begin+1 #Number of time points in final output data

  pars1=parameter_setup(FOI_spillover,R0,vacc_data,pop_data,year0,years_data,mode_start,
                        vaccine_efficacy,start_SEIRV,time_inc)
  n_years=length(pop_data[,1])-1
  inv_365=1.0/365.0
  pars2=list(FOI_spillover=pars1$FOI_spillover,R0=pars1$R0,vacc_rate_daily=pars1$vacc_rate_daily,
             vacc_cov_cam=vacc_cov_cam, t_cam=t_cam,N_age=pars1$N_age,
             S_0=pars1$S_0,E_0=pars1$E_0,I_0=pars1$I_0,R_0=pars1$R_0,V_0=pars1$V_0,
             dP1_all=pars1$dP1_all,dP2_all=pars1$dP2_all,n_years=pars1$n_years,year0=pars1$year0,vaccine_efficacy=pars1$vaccine_efficacy,
             time_inc=pars1$time_inc,t_incubation=pars1$t_incubation,t_latent=pars1$t_latent,t_infectious=pars1$t_infectious,response_delay=response_delay,
             p_rep=p_rep,case_threshold=case_threshold,cluster_threshold=cluster_threshold)

  #Carrying forward delay from previous run may cause errors
  #if(mode_start==2){pars2$E_delay0=start_SEIRV$E_delay} else {pars2$E_delay0=rep(0,nd1*N_age)}
  pars2$np_E_delay=nd1*N_age
  pars2$np_I_delay=nd2*N_age
  pars2$E_delay0=rep(0,pars2$np_E_delay)
  pars2$I_delay0=rep(0,pars2$np_I_delay)

  #Check that there is no overlap between emergency campaign and other vaccination
  #(Not yet possible to adjust vaccine rates on the fly to deal with overlap)
  for(i in 1:N_age){
    if(vacc_cov_cam[i]>0){assert_that(all(pars2$vacc_rate_daily[i,]==0))}
  }

  x <- dust_system_create(SEIRVModelDelayReactive, pars = pars2,
                          n_particles = n_particles, n_threads = n_threads, time = 0, dt = 1,
                          deterministic = deterministic, preserve_particle_dimension = TRUE)
  index=dust_unpack_index(x)
  dust_system_set_state_initial(x)
  x_res <- dust_system_simulate(x, c(step_begin:step_end))

  dims1=c(n_particles,t_pts_out)
  dims2=c(N_age,n_particles,t_pts_out)
  output_data=list(day=x_res[1,1,],year=x_res[2,1,],FOI_total=array(x_res[3,,]/time_inc,dim=dims1),C_rep_total=array(x_res[4,,],dim=dims1),
                   flag1=array(x_res[5,,],dim=dims1),flag2=array(x_res[6,,],dim=dims1),
                   flag3=array(x_res[7,,],dim=dims1),flag4=array(x_res[8,,],dim=dims1),
                   report_rate=array(x_res[9,,],dim=dims1))
  output_data$S=array(x_res[index$S,,],dim=dims2)
  output_data$E=array(x_res[index$E,,],dim=dims2)
  output_data$E_delay=array(x_res[index$E_delay,,],dim=c(N_age,nd1,n_particles,t_pts_out))
  output_data$I=array(x_res[index$I,,],dim=dims2)
  output_data$I_delay=array(x_res[index$I_delay,,],dim=c(N_age,nd2,n_particles,t_pts_out))
  output_data$R=array(x_res[index$R,,],dim=dims2)
  output_data$V=array(x_res[index$V,,],dim=dims2)
  output_data$C=array(x_res[index$C,,],dim=dims2)
  #output_data$C_rep=array(x_res[index$C_rep,,],dim=dims2)

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
#' @param FOI_spillover Force of infection due to spillover from sylvatic reservoir
#' @param R0 Basic reproduction number for urban spread of infection
#' @param vacc_data Projected vaccination-based immunity (assuming vaccine_efficacy=1) by age group and year
#' @param pop_data Population by age group and year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' @param year0 First year in population/vaccination data
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param time_inc Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param n_particles number of particles to use
#' @param n_threads number of threads to use
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' '
#' @export
#'
Model_Run_Split <- function(FOI_spillover = 0.0,R0 = 1.0,vacc_data = list(),pop_data = list(),years_data = c(1940:1941),
                            start_SEIRV = list(), year0 = 1940, mode_start = 0,
                            vaccine_efficacy = 1.0, time_inc = 1.0, n_particles = 1, n_threads = 1, deterministic = FALSE) {

  #TODO Add assert_that functions (NB - Some checks carried out in parameter_setup)
  assert_that(n_particles <= 20, msg = "Number of particles must be 20 or less")

  N_age=length(pop_data[1,]) #Number of age groups
  step_begin=((years_data[1]-year0)*(365/time_inc)) #Step at which data starts being saved for final output
  step_end=((max(years_data)+1-year0)*(365/time_inc))-1 #Step at which to end
  t_pts_out=step_end-step_begin+1 #Number of time points in final output data

  x <- dust_system_create(SEIRVModelSplitInfection, pars = parameter_setup(FOI_spillover,R0,vacc_data,pop_data,
                                                                           year0,years_data,mode_start,
                                                                           vaccine_efficacy,start_SEIRV,time_inc),
                          n_particles = n_particles, n_threads = n_threads, time = 0, dt = 1,
                          deterministic = deterministic, preserve_particle_dimension = TRUE)
  index=dust_unpack_index(x)
  dust_system_set_state_initial(x)
  x_res <- dust_system_simulate(x, c(step_begin:step_end))

  dimensions=c(N_age,n_particles,t_pts_out)
  output_data=list(day=x_res[1,1,],year=x_res[2,1,],FOI_sylvatic=array(x_res[3,,]/time_inc,dim=c(n_particles,t_pts_out)),
                   FOI_urban=array(x_res[4,,]/time_inc,dim=c(n_particles,t_pts_out)))
  output_data$S=array(x_res[index$S,,],dim=dimensions)
  output_data$E_sylvatic=array(x_res[index$E_sylvatic,,],dim=dimensions)
  output_data$E_urban = array(x_res[index$E_urban,,],dim=dimensions)
  output_data$I_sylvatic=array(x_res[index$I_sylvatic,,],dim=dimensions)
  output_data$I_urban = array(x_res[index$I_urban,,],dim=dimensions)
  output_data$R=array(x_res[index$R,,],dim=dimensions)
  output_data$V=array(x_res[index$V,,],dim=dimensions)
  output_data$C_sylvatic=array(x_res[index$C_sylvatic,,],dim=dimensions)
  output_data$C_urban = array(x_res[index$C_urban,,],dim=dimensions)

  return(output_data)
}
#-------------------------------------------------------------------------------
#TODO - Remove once split-p-rep MCMC brought into line with YEP MCMC
#' @title create_param_labels2
#'
#' @description Apply names to the parameters in a set used for data matching and parameter fitting
#'
#' @details Takes in input list and environmental data along with names of additional parameters (vaccine efficacy
#' and reporting probabilities) and generates list of names for parameter set to use as input for fitting functions
#'
#' @param type Type of parameter set (FOI only, FOI+R0, FOI and/or R0 coefficients associated with environmental
#'   covariates); choose from "FOI", "FOI+R0", "FOI enviro", "FOI+R0 enviro"
#' @param input_data List of population and vaccination data for multiple regions (created using data input creation
#' code and usually loaded from RDS file)
#' @param enviro_data Environmental data frame, containing only relevant environmental variables
#' @param extra_estimated_params Vector of strings listing variable parameters besides ones determining FOI/R0 (may include
#' vaccine efficacy and/or infection/death reporting probabilities and/or Brazil FOI adjustment factor)
#'
#' @export
#'
create_param_labels2 <- function(type = "FOI", input_data = list(), enviro_data = NULL, extra_estimated_params = c("vacc_eff")){

  assert_that(type %in% c("FOI", "FOI+R0", "FOI enviro", "FOI+R0 enviro"))
  assert_that(input_data_check(input_data), msg = paste("Input data must be in standard format",
                                                        " (see https://mrc-ide.github.io/YEP/articles/CGuideAInputs.html)"))

  n_extra = length(extra_estimated_params)

  if(type %in% c("FOI", "FOI+R0")){
    regions = input_data$region_labels
    n_regions = length(regions)
    if(type == "FOI"){n_params = n_regions+n_extra} else {n_params = (2*n_regions)+n_extra}
    param_names = rep("", n_params)
    for(i in 1:n_regions){
      param_names[i] = paste("FOI_", regions[i], sep = "")
      if(type == "FOI+R0"){param_names[i+n_regions] = paste("R0_", regions[i], sep = "")}
    }
  } else {
    assert_that(is.data.frame(enviro_data)) #TODO - msg
    env_vars = colnames(enviro_data)[c(2:ncol(enviro_data))]
    n_env_vars = length(env_vars)
    if(type == "FOI enviro"){n_params = n_env_vars+n_extra} else {n_params = (2*n_env_vars)+n_extra}
    param_names = rep("", n_params)
    for(i in 1:n_env_vars){
      param_names[i] = paste("FOI_", env_vars[i], sep = "")
      if(type == "FOI+R0 enviro"){param_names[i+n_env_vars] = paste("R0_", env_vars[i], sep = "")}
    }
  }
  if(n_extra>0){param_names[(n_params-n_extra+1):n_params] = extra_estimated_params}

  return(param_names)
}
#-------------------------------------------------------------------------------
#TODO - Remove once split-p-rep MCMC brought into line with YEP MCMC
#' @title param_calc_enviro2
#'
#' @description Parameter calculation from environmental covariates
#'
#' @details Takes in set of coefficients of environmental covariates and covariate values and calculates values of
#'   spillover force of infection and reproduction number.
#'
#' @param enviro_coeffs Values of environmental coefficients
#' @param enviro_covar_values Values of environmental covariates
#' '
#' @export
#'
param_calc_enviro2 <- function(enviro_coeffs = c(), enviro_covar_values = c()){

  assert_that(all(enviro_coeffs >= 0), msg = "All environmental coefficients must have positive values")
  n_env_vars = length(enviro_covar_values)
  assert_that(length(enviro_coeffs) %in% c(n_env_vars, 2*n_env_vars), msg = "Wrong number of environmental coefficients")

  output = list(FOI = NA, R0 = NA)
  output$FOI = sum(enviro_coeffs[c(1:n_env_vars)]*enviro_covar_values)
  if(length(enviro_coeffs) == 2*n_env_vars){
    output$R0 = sum(enviro_coeffs[c(1:n_env_vars)+n_env_vars]*enviro_covar_values)
  }

  return(output)
}
