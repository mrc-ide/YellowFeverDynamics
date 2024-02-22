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
#' @importFrom stats cov dexp dnbinom prop.test rbinom runif
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
#' @param dt Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param n_particles number of particles to use
#' @param n_threads number of threads to use
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' '
#' @export
#'
Model_Run_Delay <- function(FOI_spillover = 0.0, R0 = 1.0, vacc_data = list(), pop_data = list(), years_data = c(1940:1941),
                      start_SEIRV = list(), year0 = 1940, mode_start = 0, vaccine_efficacy = 1.0,
                      dt = 1.0, n_particles = 1, n_threads = 1, deterministic = FALSE) {

  #TODO Add assert_that functions
  assert_that(length(FOI_spillover)==1,msg="Spillover FOI must be singular value")
  assert_that(length(R0)==1,msg="R0 must be singular value")

  n_nv=3 #Number of non-vector outputs
  N_age=length(pop_data[1,]) #Number of age groups
  nd1 <- (t_incubation+t_latent)/dt
  nd2 <- t_infectious/dt
  nd=nd1+nd2
  n_data_pts=((6+nd1+nd2)*N_age)+n_nv #Number of data values per time point in output
  step_begin=((years_data[1]-year0)*(365/dt)) #Step at which data starts being saved for final output
  step_end=((max(years_data)+1-year0)*(365/dt))-1 #Step at which to end
  t_pts_out=step_end-step_begin+1 #Number of time points in final output data

  pars=parameter_setup(FOI_spillover,R0,vacc_data,pop_data,year0,years_data,mode_start,vaccine_efficacy,start_SEIRV,dt)
  #Carrying forward delay seems to result in errors - further investigation may be needed
  # if(mode_start==2){
  #   pars2$E_delay0=start_SEIRV$E_delay
  #   pars2$I_delay0=start_SEIRV$I_delay
  # } else {
  pars2$E_delay0=rep(0,nd1*N_age)
  pars2$I_delay0=rep(0,nd2*N_age)
  # }

  x <- SEIRVModelDelay$new(pars,time = 0, n_particles = n_particles, n_threads = n_threads, deterministic = deterministic)

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
#' @title Model_Run_Delay_Reactive
#'
#' @description Runs delay+reactive version of SEIRV model
#'
#' @details Accepts epidemiological + population parameters and model settings; runs delay/reactive SEIRV model
#' for one region over a specified time period for a number of particles/threads and outputs time-dependent SEIRV
#' values, infection numbers and total force of infection values. This version of the model differs from the standard
#' one in taking two sets of input vaccination data, one the default one and one applied after one or more cases have
#' been reported (as well as using delay instead of rate for incubation, infectious period etc.. Case reporting is governed
#' by an additional parameter p_rep which can also change after one or more cases have been reported in order to reflect changes
#' in surveillance.
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
#' @param dt Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param n_particles number of particles to use
#' @param n_threads number of threads to use
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' @param response_delay Delay time in days between a threshold being reached and emergency conditions coming into effect
#' @param p_rep Probabilities of an infection being reported as a case before emergency conditions triggered (1st value) or
#'   after emergency conditions triggered (2nd value)
#' @param case_threshold Threshold total no. reported cases to trigger emergency conditions
#' @param cluster_threshold Threshold current infectious fraction to trigger emergency conditions
#' @param vacc_rate_cam Vaccination rate by age group during emergency vaccination campaign
#' @param t_cam Duration in days of emergency vaccination campaign
#' '
#' @export
#'
Model_Run_Delay_Reactive <- function(FOI_spillover = 0.0,R0 = 1.0,vacc_data = list(), pop_data = list(),
                                     years_data = c(1940:1941), start_SEIRV = list(), year0 = 1940, mode_start = 0,
                                     vaccine_efficacy = 1.0, dt = 1.0, n_particles = 1, n_threads = 1, deterministic = FALSE,
                                     response_delay = 56.0, p_rep = c(0.0,0.0), case_threshold = Inf,
                                     cluster_threshold = Inf, vacc_rate_cam = c(), t_cam = 0) {

  #TODO Add assert_that functions
  assert_that(length(FOI_spillover)==1,msg="Spillover FOI must be singular value")
  assert_that(length(R0)==1,msg="R0 must be singular value")

  n_nv=9 #Number of non-vector outputs
  N_age=length(pop_data[1,]) #Number of age groups
  assert_that(length(vacc_rate_cam)==N_age)
  nd1 <- (t_incubation+t_latent)/dt
  nd2 <- t_infectious/dt
  nd <- nd1 + nd2
  n_data_pts=((6+nd1+nd2)*N_age)+n_nv #Number of data values per time point in output
  step_begin=((years_data[1]-year0)*(365/dt)) #Step at which data starts being saved for final output
  step_end=((max(years_data)+1-year0)*(365/dt))-1 #Step at which to end
  t_pts_out=step_end-step_begin+1 #Number of time points in final output data

  pars1=parameter_setup(FOI_spillover,R0,vacc_data,pop_data,year0,years_data,mode_start,
                        vaccine_efficacy,start_SEIRV,dt)
  n_years=length(pop_data[,1])-1
  inv_365=1.0/365.0
  pars2=list(FOI_spillover=pars1$FOI_spillover,R0=pars1$R0,vacc_rate_daily=pars1$vacc_rate_daily,
             vacc_rate_cam=vacc_rate_cam, t_cam=t_cam,N_age=pars1$N_age,
             S_0=pars1$S_0,E_0=pars1$E_0,I_0=pars1$I_0,R_0=pars1$R_0,V_0=pars1$V_0,
             dP1_all=pars1$dP1_all,dP2_all=pars1$dP2_all,n_years=pars1$n_years,year0=pars1$year0,vaccine_efficacy=pars1$vaccine_efficacy,
             dt=pars1$dt,t_incubation=pars1$t_incubation,t_latent=pars1$t_latent,t_infectious=pars1$t_infectious,response_delay=response_delay,
             p_rep=p_rep,case_threshold=case_threshold,cluster_threshold=cluster_threshold)
  #Carrying forward delay seems to result in errors - further investigation may be needed
  # if(mode_start==2){
  #   pars2$E_delay0=start_SEIRV$E_delay
  #   pars2$I_delay0=start_SEIRV$I_delay
  # } else {
    pars2$E_delay0=rep(0,nd1*N_age)
    pars2$I_delay0=rep(0,nd2*N_age)
  # }

  #Check that there is no overlap between emergency campaign and other vaccination
  #(Not yet possible to adjust vaccine rates on the fly to deal with overlap)
  for(i in 1:N_age){
    if(vacc_rate_cam[i]>0){assert_that(all(pars2$vacc_rate_daily[i,]==0))}
  }

  x <- SEIRVModelDelayReactive$new(pars=pars2,time = 0, n_particles = n_particles, n_threads = n_threads, deterministic = deterministic)

  x_res <- array(NA, dim = c(n_data_pts, n_particles, t_pts_out))
  for(step in step_begin:step_end){
    x_res[,,step-step_begin+1] <- x$run(step)
  }
  if(step_begin==0){x_res[2,,1]=rep(year0,n_particles)}

  dims1=c(n_particles,t_pts_out)
  dims2=c(N_age,n_particles,t_pts_out)
  output_data=list(day=x_res[1,1,],year=x_res[2,1,],FOI_total=array(x_res[3,,]/dt,dim=dims1),C_rep_total=array(x_res[4,,],dim=dims1),
                   flag1=array(x_res[5,,],dim=dims1),flag2=array(x_res[6,,],dim=dims1),
                   flag3=array(x_res[7,,],dim=dims1),flag4=array(x_res[8,,],dim=dims1),
                   report_rate=array(x_res[9,,],dim=dims1))
  output_data$S=array(x_res[c((1+n_nv):(N_age+n_nv)),,],dim=dims2)
  output_data$E=array(x_res[c((N_age+1+n_nv):((2*N_age)+n_nv)),,],dim=dims2)
  output_data$E_delay=array(x_res[c(((2*N_age)+1+n_nv):(((2+nd1)*N_age)+n_nv)),,],dim=c(N_age,nd1,n_particles,t_pts_out))
  output_data$I=array(x_res[c((((2+nd1)*N_age)+1+n_nv):(((3+nd1)*N_age)+n_nv)),,],dim=dims2)
  output_data$I_delay=array(x_res[c((((3+nd1)*N_age)+1+n_nv):(((3+nd)*N_age)+n_nv)),,],dim=c(N_age,nd2,n_particles,t_pts_out))
  output_data$R=array(x_res[c((((3+nd)*N_age)+1+n_nv):(((4+nd)*N_age)+n_nv)),,],dim=dims2)
  output_data$V=array(x_res[c((((4+nd)*N_age)+1+n_nv):(((5+nd)*N_age)+n_nv)),,],dim=dims2)
  output_data$C=array(x_res[c((((5+nd)*N_age)+1+n_nv):(((6+nd)*N_age)+n_nv)),,],dim=dims2)
  #output_data$C_rep=array(x_res[c((((6+nd)*N_age)+1+n_nv):(((7+nd)*N_age)+n_nv)),,],dim=dims2)

  return(output_data)
}
#-------------------------------------------------------------------------------
#' @title Model_Run_Reactive
#'
#' @description Runs reactive version of SEIRV model
#'
#' @details Accepts epidemiological + population parameters and model settings; runs reactive SEIRV model
#' for one region over a specified time period for a number of particles/threads and outputs time-dependent SEIRV
#' values, infection numbers and total force of infection values. This version of the model differs from the standard
#' one in taking two sets of input vaccination data, one the default one and one applied after one or more cases have
#' been reported. Case reporting is governed by an additional parameter p_rep which can also change after one or more
#' cases have been reported in order to reflect changes in surveillance.
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
#' @param dt Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param n_particles number of particles to use
#' @param n_threads number of threads to use
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' @param response_delay Delay time in days between a threshold being reached and emergency conditions coming into effect
#' @param p_rep Probabilities of an infection being reported as a case before emergency conditions triggered (1st value) or
#'   after emergency conditions triggered (2nd value)
#' @param case_threshold Threshold total no. reported cases to trigger emergency conditions
#' @param cluster_threshold Threshold current infectious fraction to trigger emergency conditions
#' @param vacc_rate_cam Vaccination rate by age group during emergency vaccination campaign
#' @param t_cam Duration in days of emergency vaccination campaign
#' '
#' @export
#'
Model_Run_Reactive <- function(FOI_spillover = 0.0,R0 = 1.0,vacc_data = list(), pop_data = list(),
                               years_data = c(1940:1941), start_SEIRV = list(), year0 = 1940, mode_start = 0,
                               vaccine_efficacy = 1.0, dt = 1.0, n_particles = 1, n_threads = 1, deterministic = FALSE,
                               response_delay = 56.0, p_rep = c(0.0,0.0), case_threshold = Inf,
                               cluster_threshold = Inf, vacc_rate_cam = c(), t_cam = 0) {

  #TODO Add assert_that functions
  assert_that(length(FOI_spillover)==1,msg="Spillover FOI must be singular value")
  assert_that(length(R0)==1,msg="R0 must be singular value")

  pars1=parameter_setup(FOI_spillover,R0,vacc_data,pop_data,year0,years_data,mode_start,
                        vaccine_efficacy,start_SEIRV,dt)
  n_years=length(pop_data[,1])-1
  inv_365=1.0/365.0
  pars2=list(FOI_spillover=pars1$FOI_spillover,R0=pars1$R0,vacc_rate_daily=pars1$vacc_rate_daily,
             vacc_rate_cam=vacc_rate_cam, t_cam=t_cam,N_age=pars1$N_age,
             S_0=pars1$S_0,E_0=pars1$E_0,I_0=pars1$I_0,R_0=pars1$R_0,V_0=pars1$V_0,
             dP1_all=pars1$dP1_all,dP2_all=pars1$dP2_all,n_years=pars1$n_years,year0=pars1$year0,vaccine_efficacy=pars1$vaccine_efficacy,
             dt=pars1$dt,t_incubation=pars1$t_incubation,t_latent=pars1$t_latent,t_infectious=pars1$t_infectious,response_delay=response_delay,
             p_rep=p_rep,case_threshold=case_threshold,cluster_threshold=cluster_threshold)

  n_nv=9 #Number of non-vector outputs
  N_age=length(pop_data[1,]) #Number of age groups
  assert_that(length(vacc_rate_cam)==N_age)
  n_data_pts=(6*N_age)+n_nv #Number of data values per time point in output
  step_begin=((years_data[1]-year0)*(365/dt)) #Step at which data starts being saved for final output
  step_end=((max(years_data)+1-year0)*(365/dt))-1 #Step at which to end
  t_pts_out=step_end-step_begin+1 #Number of time points in final output data

  #Check that there is no overlap between emergency campaign and other vaccination
  #(Not yet possible to adjust vaccine rates on the fly to deal with overlap)
  for(i in 1:N_age){
    if(vacc_rate_cam[i]>0){assert_that(all(pars2$vacc_rate_daily[i,]==0))}
  }

  x <- SEIRVModelReactive$new(pars=pars2,time = 0, n_particles = n_particles, n_threads = n_threads, deterministic = deterministic)

  x_res <- array(NA, dim = c(n_data_pts, n_particles, t_pts_out))
  for(step in step_begin:step_end){
    x_res[,,step-step_begin+1] <- x$run(step)
  }
  if(step_begin==0){x_res[2,,1]=rep(year0,n_particles)}

  dims1=c(n_particles,t_pts_out)
  dims2=c(N_age,n_particles,t_pts_out)
  output_data=list(day=x_res[1,1,],year=x_res[2,1,],FOI_total=array(x_res[3,,]/dt,dim=dims1),C_rep_total=array(x_res[4,,],dim=dims1),
                   flag1=array(x_res[5,,],dim=dims1),flag2=array(x_res[6,,],dim=dims1),
                   flag3=array(x_res[7,,],dim=dims1),flag4=array(x_res[8,,],dim=dims1),
                   report_rate=array(x_res[9,,],dim=dims1))
  output_data$S=array(x_res[c((1+n_nv):(N_age+n_nv)),,],dim=dims2)
  output_data$E=array(x_res[c((N_age+1+n_nv):((2*N_age)+n_nv)),,],dim=dims2)
  output_data$I=array(x_res[c(((2*N_age)+1+n_nv):((3*N_age)+n_nv)),,],dim=dims2)
  output_data$R=array(x_res[c(((3*N_age)+1+n_nv):((4*N_age)+n_nv)),,],dim=dims2)
  output_data$V=array(x_res[c(((4*N_age)+1+n_nv):((5*N_age)+n_nv)),,],dim=dims2)
  output_data$C=array(x_res[c(((5*N_age)+1+n_nv):((6*N_age)+n_nv)),,],dim=dims2)
  #output_data$C_rep=array(x_res[c(((6*N_age)+1+n_nv):((7*N_age)+n_nv)),,],dim=dims2)

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
  assert_that(length(FOI_spillover)==1,msg="Spillover FOI must be singular value")
  assert_that(length(R0)==1,msg="R0 must be singular value")

  x <- SEIRVModelSplitInfection$new(pars=parameter_setup(FOI_spillover,R0,vacc_data,pop_data,year0,years_data,mode_start,
                                                 vaccine_efficacy,start_SEIRV,dt),
                            time = 0, n_particles = n_particles, n_threads = n_threads, deterministic = deterministic)

  n_nv=4 #Number of non-vector outputs
  N_age=length(pop_data[1,]) #Number of age groups
  n_data_pts=(9*N_age)+n_nv #Number of data values per time point in output
  step_begin=((years_data[1]-year0)*(365/dt)) #Step at which data starts being saved for final output
  step_end=((max(years_data)+1-year0)*(365/dt))-1 #Step at which to end
  t_pts_out=step_end-step_begin+1 #Number of time points in final output data

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
#' @details Accepts epidemiological + population parameters and model settings; runs version of SEIRV model with
#' time-varying FOI_spillover and R0 values for one region over a specified time period for a number of
#' particles/threads and outputs time-dependent SEIRV values, infection numbers and total force of infection values.
#'
#' @param FOI_spillover Vector of values of force of infection due to spillover from sylvatic reservoir (size depends on mode_time)
#' @param R0 Vector of values of basic reproduction number for urban spread of infection (size depends on mode_time)
#' @param vacc_data Projected vaccination-based immunity (assuming vaccine_efficacy=1) by age group and year
#' @param pop_data Population by age group and year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' @param output_type Type of data to output:
#'   "full" = SEIRVC + FOI for all steps and ages
#'   "case" = annual total new infections (C) summed across all ages
#'   "sero" = annual SEIRV
#'   "case+sero" = annual SEIRVC, cases summed across all ages
#'   "case_alt" = annual total new infections not combined by age
#'   "case_alt2" = total new infections combined by age for all steps
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
#' @param mode_time Type of time dependence of FOI_spillover and R0 to be used:
#'  If mode_time=0, no time variation (identical to Model_Run())
#'  If mode_time=1, FOI/R0 vary annually without seasonality (number of values = number of years to consider)
#'  If mode_time=2, FOI/R0 vary with monthly seasonality without inter-annual variation (number of values = 12)
#'  If mode_time=3, FOI/R0 vary with daily seasonality without inter-annual variation (number of values = 365/dt)
#'  If mode_time=4, FOI/R0 vary annually with monthly seasonality (number of values = 12*number of years to consider)
#'  If mode_time=5, FOI/R0 vary annually with daily seasonality (number of values = (365/dt)*number of years to consider)
#' '
#' @export
#'
Model_Run_VarFR <- function(FOI_spillover = c(), R0 = c(), vacc_data = list(),pop_data = list(),years_data = c(1940:1941),
                            start_SEIRV = list(), output_type = "full", year0 = 1940, mode_start = 0, vaccine_efficacy = 1.0,
                            dt = 1.0, n_particles = 1, n_threads = 1, deterministic = FALSE, mode_time = 1) {

  #TODO Add additional assert_that functions (NB - Some checks carried out in parameter_setup)
  assert_that(n_particles<=20,msg="Number of particles must be 20 or less")
  pts_year=365/dt
  assert_that(mode_time %in% c(0:5))

  pars1=parameter_setup(FOI_spillover[1],R0[1],vacc_data,pop_data,year0,years_data,mode_start,
                        vaccine_efficacy,start_SEIRV,dt)
  pts_total=pars1$n_years*pts_year

  #Check length of FOI/R0 input and generate values for every time point
  if(mode_time==0){
    assert_that(length(FOI_spillover)==1 && length(R0)==1,msg="Spillover FOI and R0 must be single values if mode_time=0")
    FOI_spillover_t=rep(FOI_spillover,pts_total)
    R0_t=rep(R0,pts_total)
  } else {
    FOI_spillover_t=R0_t=rep(NA,pts_total)
  }
  if(mode_time==1){
    assert_that(length(FOI_spillover)==pars1$n_years && length(R0)==pars1$n_years,
                msg="Spillover FOI and R0 must be vectors of length equal to no.years considered if mode_time=1")
    for(i in 1:pars1$n_years){
      FOI_spillover_t[c(1:pts_year)+(i-1)*pts_year]=rep(FOI_spillover[i],pts_year)
      R0_t[c(1:pts_year)+(i-1)*pts_year]=rep(R0[i],pts_year)
    }
  }
  if(mode_time==2){
    assert_that(length(FOI_spillover)==12 && length(R0)==12, msg="Spillover FOI and R0 must be vectors of length 12 if mode_time=2")
    month_values=(1+floor(((12*dt)/365)*c(0:(pts_total-1))) %% 12)
    FOI_spillover_t=FOI_spillover[month_values]
    R0_t=R0[month_values]
  }
  if(mode_time==3){
    assert_that(length(FOI_spillover)==pts_year && length(R0)==pts_year,
                msg="Spillover FOI and R0 must be vectors of length (365/dt) if mode_time=3")
    date_values=1+(floor(dt*c(0:(pts_total-1))) %% (365/dt))
    FOI_spillover_t=FOI_spillover[date_values]
    R0_t=R0[date_values]
  }
  if(mode_time==4){
    assert_that(length(FOI_spillover)==pars1$n_years*12 && length(R0)==pars1$n_years*12,
                msg="Spillover FOI and R0 must be vectors of length equal to 12*no.years considered if mode_time=4")
    #month_values=1+floor(((12*dt)/365)*c(0:(pts_total-1)))
    month_values=(1+floor(((12*dt)/365)*c(0:(pts_total-1))) %% 12)
    date_values=month_values+sort(rep(c(1:pars1$n_years)-1,pts_year))
    FOI_spillover_t=FOI_spillover[date_values]
    R0_t=R0[date_values]
  }
  if(mode_time==5){
    assert_that(length(FOI_spillover)==pars1$n_years*pts_year && length(R0)==pars1$n_years*pts_year,
                msg="Spillover FOI and R0 must be vectors of length equal to (365/dt)*no.years considered if mode_time=5")
    FOI_spillover_t=FOI_spillover
    R0_t=R0
  }

  pars2=list(FOI_spillover=FOI_spillover_t,R0=R0_t,vacc_rate_daily=pars1$vacc_rate_daily,N_age=pars1$N_age,
             S_0=pars1$S_0,E_0=pars1$E_0,I_0=pars1$I_0,R_0=pars1$R_0,V_0=pars1$V_0,dP1_all=pars1$dP1_all,dP2_all=pars1$dP2_all,
             n_years=pars1$n_years, n_t_pts=pts_total,year0=pars1$year0,vaccine_efficacy=pars1$vaccine_efficacy,
             dt=pars1$dt,t_incubation=pars1$t_incubation,t_latent=pars1$t_latent,t_infectious=pars1$t_infectious)

  n_nv=3 #Number of non-vector outputs
  N_age=length(pop_data[1,]) #Number of age groups
  n_data_pts=(6*N_age)+n_nv #Number of data values per time point in output
  step_begin=(years_data[1]-year0)*pts_year #Step at which data starts being saved for final output
  step_end=((max(years_data)+1-year0)*pts_year)-1 #Step at which to end
  t_pts_out=step_end-step_begin+1 #Number of time points in final output data

  x <- SEIRVModelVarFR$new(pars=pars2,time = 0, n_particles = n_particles, n_threads = n_threads, deterministic = deterministic)

  x_res <- array(NA, dim = c(n_data_pts, n_particles, t_pts_out))
  for(step in step_begin:step_end){
    x_res[,,step-step_begin+1] <- x$run(step)
  }
  if(step_begin==0){x_res[2,,1]=rep(year0,n_particles)}

  if(output_type=="full"){
    dimensions=c(N_age,n_particles,t_pts_out)
    output_data=list(day=x_res[1,1,],year=x_res[2,1,])
    output_data$FOI_total=array(x_res[3,,]/dt,dim=c(n_particles,t_pts_out))
    output_data$S=array(x_res[c((1+n_nv):(N_age+n_nv)),,],dim=dimensions)
    output_data$E=array(x_res[c((N_age+1+n_nv):((2*N_age)+n_nv)),,],dim=dimensions)
    output_data$I=array(x_res[c(((2*N_age)+1+n_nv):((3*N_age)+n_nv)),,],dim=dimensions)
    output_data$R=array(x_res[c(((3*N_age)+1+n_nv):((4*N_age)+n_nv)),,],dim=dimensions)
    output_data$V=array(x_res[c(((4*N_age)+1+n_nv):((5*N_age)+n_nv)),,],dim=dimensions)
    output_data$C=array(x_res[c(((5*N_age)+1+n_nv):((6*N_age)+n_nv)),,],dim=dimensions)
  } else {
    if(output_type=="case_alt2"){
      output_data=list(day=x_res[1,1,],year=x_res[2,1,])
      output_data$C=array(0,dim=c(n_particles,t_pts_out))
      for(pt in 1:t_pts_out){
        for(n_p in 1:n_particles){
          output_data$C[n_p,pt]=sum(x_res[c(((5*N_age)+1+n_nv):((6*N_age)+n_nv)),n_p,pt])
        }
      }
    }  else {
      n_years=length(years_data)
      output_data=list(year=years_data)
      if(output_type=="case+sero" || output_type=="sero"){
        output_data$V=output_data$R=output_data$I=output_data$E=output_data$S=array(0,dim=c(N_age,n_particles,n_years))
        for(n_year in 1:n_years){
          pts=c(1:t_pts_out)[x_res[2,1,]==years_data[n_year]]
          for(n_p in 1:n_particles){
            output_data$S[,n_p,n_year]=rowMeans(x_res[c((1+n_nv):(N_age+n_nv)),n_p,pts])
            output_data$E[,n_p,n_year]=rowMeans(x_res[c((N_age+1+n_nv):((2*N_age)+n_nv)),n_p,pts])
            output_data$I[,n_p,n_year]=rowMeans(x_res[c(((2*N_age)+1+n_nv):((3*N_age)+n_nv)),n_p,pts])
            output_data$R[,n_p,n_year]=rowMeans(x_res[c(((3*N_age)+1+n_nv):((4*N_age)+n_nv)),n_p,pts])
            output_data$V[,n_p,n_year]=rowMeans(x_res[c(((4*N_age)+1+n_nv):((5*N_age)+n_nv)),n_p,pts])
          }
        }
      }
      if(output_type=="case+sero" || output_type=="case"){
        output_data$C=array(0,dim=c(n_particles,n_years))
        for(n_year in 1:n_years){
          pts=c(1:t_pts_out)[x_res[2,1,]==years_data[n_year]]
          for(n_p in 1:n_particles){
            output_data$C[n_p,n_year]=sum(x_res[c(((5*N_age)+1+n_nv):((6*N_age)+n_nv)),n_p,pts])
          }
        }
      }
      if(output_type=="case_alt"){
        output_data$C=array(0,dim=c(N_age,n_particles,n_years))
        for(n_year in 1:n_years){
          pts=c(1:t_pts_out)[x_res[2,1,]==years_data[n_year]]
          for(n_p in 1:n_particles){
            output_data$C[,n_p,n_year]=rowSums(x_res[c(((5*N_age)+1+n_nv):((6*N_age)+n_nv)),n_p,pts])
          }
        }
      }
    }
  }

  return(output_data)
}
