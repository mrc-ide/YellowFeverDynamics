# R file for general functions in YellowFeverDynamics package
#------------------------------------------------
#Global variables
p_severe_inf=0.12 #Probability that an infection is severe
p_death_severe_inf=0.47 #Probability that a severe infection becomes fatal
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.
#' @useDynLib YellowFeverDynamics, .registration = TRUE
#' @import assertthat
#' @import ggplot2
#' @import graphics
#' @import maptree
#' @import mvtnorm
#' @import qpdf
#' @import Rmisc
#' @import stats
#' @import testthat
#' @import tgp
#' @import truncdist
#' @import utils
#------------------------------------------------
# unload DLL when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("YellowFeverDynamics", libpath)
}
#-------------------------------------------------------------------------------
#' @title Full_Model_Run
#'
#' @description Run full version of SEIRV model
#'
#' @details Accepts epidemiological + population parameters and model settings; runs full version of SEIRV model
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
Full_Model_Run <- function(FOI_spillover=0.0,R0=1.0,vacc_data=list(),pop_data=list(),year0=1940,mode_start=0,
                           n_particles=1,n_threads=1,year_end=2000,year_data_begin=1999,vaccine_efficacy=1.0,
                           start_SEIRV=list(),dt=1.0) {

  assert_that(n_particles>0)
  assert_that(n_particles<=20)
  assert_that(n_threads<=n_particles)
  assert_that(n_threads>0)

  x <- FullModelOD$new(pars=parameter_setup(FOI_spillover,R0,vacc_data,pop_data,year0,mode_start,year_end,
                                          year_data_begin,vaccine_efficacy,start_SEIRV,dt),
                     step = 1,n_particles = n_particles,n_threads = n_threads)

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
#' @title Basic_Model_Run
#'
#' @description Run basic version of SEIRV model
#'
#' @details Accepts epidemiological + population parameters and model settings; runs basic version of SEIRV model
#' for one region over a specified time period for a number of particles/threads and outputs time-dependent SEIRV
#' values. (Differs from Full_Model_Run in not recording new infection numbers or total force of infection.)
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
Basic_Model_Run <- function(FOI_spillover=0.0,R0=1.0,vacc_data=list(),pop_data=list(),year0=1940,mode_start=0,
                            n_particles=1,n_threads=1,year_end=2000,year_data_begin=1999,vaccine_efficacy=1.0,
                            start_SEIRV=list(),dt=1.0) {

  assert_that(n_particles>0)
  assert_that(n_particles<=20)
  assert_that(n_threads<=n_particles)
  assert_that(n_threads>0)

  x <- BasicModelOD$new(pars=parameter_setup(FOI_spillover,R0,vacc_data,pop_data,year0,mode_start,year_end,
                                          year_data_begin,vaccine_efficacy,start_SEIRV,dt),
                     step = 1,n_particles = n_particles,n_threads = n_threads)

  n_nv=3 #Number of non-vector outputs
  N_age=length(pop_data[1,]) #Number of age groups
  t_pts_all=c(1:((year_end-year0)*(365/dt))) #All output time points
  n_data_pts=(5*N_age)+n_nv #Number of data values per time point in output
  n_steps=length(t_pts_all) #Total number of output time points
  step0=(year_data_begin-year0)*(365/dt) #Step at which data starts being saved for final output
  t_pts_out=n_steps-step0 #Number of time points in final output data
  x_res <- array(NA, dim = c(n_data_pts, n_particles, t_pts_out))
  for(t in step0:n_steps){
    x_res[,,t-step0] <- x$run(t)
  }

  return(list(day=array(x_res[2,,],dim=c(n_particles,t_pts_out)),
              year=array(x_res[3,,],dim=c(n_particles,t_pts_out)),
              S=array(x_res[c((1+n_nv):(N_age+n_nv)),,],dim=c(N_age,n_particles,t_pts_out)),
              E=array(x_res[c((N_age+1+n_nv):((2*N_age)+n_nv)),,],dim=c(N_age,n_particles,t_pts_out)),
              I=array(x_res[c(((2*N_age)+1+n_nv):((3*N_age)+n_nv)),,],dim=c(N_age,n_particles,t_pts_out)),
              R=array(x_res[c(((3*N_age)+1+n_nv):((4*N_age)+n_nv)),,],dim=c(N_age,n_particles,t_pts_out)),
              V=array(x_res[c(((4*N_age)+1+n_nv):((5*N_age)+n_nv)),,],dim=c(N_age,n_particles,t_pts_out))))
}
#-------------------------------------------------------------------------------
#' @title Parameter setup
#'
#' @description Set up parameters to input into model
#'
#' @details Takes in multiple inputs, outputs list for use by odin.dust SEIRV model versions.
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
#' @param year_end year to run up to
#' @param year_data_begin year to begin saving data
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' @param dt Time increment in days to use in model (should be either 1.0 or 5.0 days)
#' '
#' @export
#'
parameter_setup <- function(FOI_spillover=0.0,R0=1.0,vacc_data=list(),pop_data=list(),year0=1940,mode_start=0,
                            year_end=2000,year_data_begin=1999,vaccine_efficacy=1.0,start_SEIRV=list(),dt=1.0){

  assert_that(length(pop_data[,1])>1)
  assert_that(length(pop_data[1,])>1)
  n_years=length(pop_data[,1])-1
  N_age=length(pop_data[1,])
  assert_that(length(vacc_data[,1])==n_years+1)
  assert_that(length(vacc_data[1,])==N_age)
  assert_that(mode_start %in% c(0,1,2))
  if(mode_start==2){assert_that(is.null(start_SEIRV$S)==FALSE)}
  assert_that(year_data_begin>=year0)
  assert_that(year_data_begin<year_end)
  assert_that(year_end-year0<=n_years)
  vacc_initial=vacc_data[1,]
  assert_that(dt %in% c(1,5))
  inv_365=1.0/365.0

  P0=Cas0=Sus0=Exp0=Inf0=Rec0=Vac0=rep(0,N_age)
  dP1_all=dP2_all=vacc_rates=array(rep(0,N_age*n_years),dim=c(N_age,n_years))
  for(i in 1:N_age){
    P0[i]=max(1.0,pop_data[1,i]) #Set all population values to a minimum of 1 to avoid NaN values appearing
  }
  for(n_year in 1:n_years){
    for(i in 1:N_age){
      dP1_all[i,n_year]=max(1.0,pop_data[n_year+1,i])*inv_365
      dP2_all[i,n_year]=max(1.0,pop_data[n_year,i])*inv_365
      if(i==1){
        vacc_rates[i,n_year]=vacc_data[n_year+1,i]*inv_365
      } else {
        vacc_rates[i,n_year]=max(0.0,vacc_data[n_year+1,i]-vacc_data[n_year,i-1])*inv_365
      }
    }
  }

  if(mode_start==0){
    Sus0=P0*(1.0-vacc_initial)
  }
  if(mode_start==1)
  {
    if(R0>1.0){
      herd_immunity=1.0-(1.0/R0)
    } else {
      herd_immunity=0.0
    }
    for(i in 1:N_age){
      if(vacc_initial[i]<herd_immunity){
        Rec0[i]=P0[i]*(herd_immunity-vacc_initial[i])
        Sus0[i]=P0[i]*(1.0-herd_immunity)
      } else {
        Sus0[i]=P0[i]*(1.0-vacc_initial[i])
      }
    }
  }
  if(mode_start==2){
    Sus0=start_SEIRV$S
    Exp0=start_SEIRV$E
    Inf0=start_SEIRV$I
    Rec0=start_SEIRV$R
    Vac0=start_SEIRV$V
    Cas0=rep(0,N_age)
  } else {
    assert_that(length(vacc_initial)==N_age)
    Vac0=P0*vacc_initial
  }

  return(list(FOI_spillover=FOI_spillover,R0=R0,vacc_rate_annual=vacc_rates,
              Cas0=Cas0,Exp0=Exp0,Inf0=Inf0,N_age=N_age,Rec0=Rec0,Sus0=Sus0,Vac0=Vac0,
              dP1_all=dP1_all,dP2_all=dP2_all,n_years=n_years,year0=year0,vaccine_efficacy=vaccine_efficacy,dt=dt))
}
#-------------------------------------------------------------------------------
#' @title mcmc_FOI_R0_setup
#'
#' @description Set up FOI and R0 values and calculate some prior probability values for MCMC calculation
#'
#' @details Takes in parameter values used for Markov Chain Monte Carlo fitting, calculates spillover force of
#' infection and (optionally) reproduction number values either directly or from environmental covariates. Also
#' calculates related components of prior probability.
#'
#' @param type Type of parameter set (FOI only, FOI+R0, FOI and/or R0 coefficients associated with environmental
#'   covariates); choose from "FOI","FOI+R0","FOI enviro","FOI+R0 enviro"
#' @param prior_type Text indicating which type of calculation to use for prior probability
#'  If prior_type = "zero", prior probability is always zero
#'  If prior_type = "flat", prior probability is zero if FOI/R0 in designated ranges, -Inf otherwise
#'  If prior_type = "exp", prior probability is given by dexp calculation on FOI/R0 values
#'  If prior_type = "norm", prior probability is given by dnorm calculation on parameter values
#' @param regions Vector of region names
#' @param param_prop Proposed parameter values
#' @param enviro_data Environmental data frame, containing only relevant environmental variables
#' @param R0_fixed_values Values of R0 to use if not being fitted
#' @param pars_min Lower limits of parameter values if specified
#' @param pars_max Upper limits of parameter values if specified
#' '
#' @export
#'
mcmc_FOI_R0_setup <- function(type="",prior_type="",regions="",param_prop=c(),enviro_data=list(),R0_fixed_values=c(),
                              pars_min=c(),pars_max=c()){

  n_params=length(param_prop)
  n_regions=length(regions)
  if(type %in% c("FOI+R0 enviro","FOI enviro")){n_env_vars=ncol(enviro_data)-1}
  FOI_values=R0_values=rep(0,n_regions)

  if(type %in% c("FOI+R0 enviro","FOI enviro")){
    for(i in 1:n_regions){
      model_params=param_calc_enviro(param=param_prop,enviro_data=enviro_data[enviro_data$adm1==regions[i],])
      FOI_values[i]=model_params$FOI
      if(type=="FOI+R0 enviro"){R0_values[i]=model_params$R0} else {R0_values[i]=R0_fixed_values[i]}
    }
  }
  if(type %in% c("FOI+R0","FOI")){
    FOI_values=exp(param_prop[c(1:n_regions)])
    if(type=="FOI+R0"){R0_values=exp(param_prop[c((n_regions+1):(2*n_regions))])
    } else {R0_values=R0_fixed_values}
  }

  prior=0
  if(prior_type=="exp"){
    prior_FOI=dexp(FOI_values,rate=1,log=TRUE)
    if(type %in% c("FOI+R0","FOI+R0 enviro")){prior_R0=dexp(R0_values,rate=1,log=TRUE)} else {prior_R0=0}
    prior = prior+sum(prior_FOI)+sum(prior_R0)
  }
  if(prior_type=="flat"){
    if(is.null(pars_min)==FALSE){
      for(i in 1:n_params){
        if(param_prop[i]<pars_min[i]){prior=-Inf}
      }
    }
    if(is.null(pars_max)==FALSE){
      for(i in 1:n_params){
        if(param_prop[i]>pars_max[i]){prior=-Inf}
      }
    }
  }
  if(prior_type=="norm"){
    if(type=="FOI"){n_params_check=n_regions}
    if(type=="FOI+R0"){n_params_check=2*n_regions}
    if(type=="FOI enviro"){n_params_check=n_env_vars}
    if(type=="FOI+R0 enviro"){n_params_check=2*n_env_vars}
    prior=sum(dnorm(param_prop[c(1:n_params_check)],mean = 0,sd = 30,log = TRUE))
  }

  output=list(FOI_values=FOI_values,R0_values=R0_values,prior=prior)
  return(output)
}
#-------------------------------------------------------------------------------
#' @title param_calc_enviro
#'
#' @description Parameter calculation from environmental variables
#'
#' @details Takes in parameter set used for Markov Chain Monte Carlo fitting and calculates values of spillover
#' force of infection and reproduction number.
#'
#' @param param Parameter values
#' @param enviro_data Environmental data frame line, containing only relevant environmental variables
#' '
#' @export
#'
param_calc_enviro <- function(param=c(),enviro_data=c()){

  n_vars=dim(enviro_data)[2]-1
  variable_names=names(enviro_data)[c(2:(n_vars+1))]
  output=list(FOI=0.0,R0=0.0)

  for(i in 1:n_vars){
    variable=variable_names[i]
    variable_value=as.numeric(enviro_data[[variable]])
    output$FOI=output$FOI+(variable_value*exp(param[i]))
    if(i+n_vars<=length(param)){output$R0=output$R0+(variable_value*exp(param[i+n_vars]))}
  }

  return(output)
}
#-------------------------------------------------------------------------------
#' @title param_prop_setup
#'
#' @description Set up proposed new parameter values for next step in chain
#'
#' @details Takes in current values of parameter set used for Markov Chain Monte Carlo fitting and proposes new values
#' from multivariate normal distribution where the existing values form the mean and the standard deviation is
#' based on the chain covariance or (if the flag "adapt" is set to 1) a flat value based on the number of parameters.
#'
#' @param param Previous parameter values used as input
#' @param chain_cov Covariance calculated from previous steps in chain
#' @param adapt 0/1 flag indicating which type of calculation to use for proposition value
#' '
#' @export
#'
param_prop_setup <- function(param=c(),chain_cov=1,adapt=0){

  n_params = length(param)
  if (adapt==1) {
    sigma = (2.38 ^ 2) * chain_cov / n_params #'optimal' scaling of chain covariance
    param_prop_a = rmvnorm(n = 1, mean = param, sigma = sigma)
  } else {
    sigma = ((1e-2) ^ 2) * diag(n_params) / n_params #this is an inital proposal covariance, see [Mckinley et al 2014]
    param_prop_a = rmvnorm(n = 1, mean = param, sigma = sigma)
  }
  param_prop = param_prop_a[1,]
  names(param_prop)=names(param)

  return(param_prop)
}
#-------------------------------------------------------------------------------
#' @title create_param_labels
#'
#' @description Apply names to the parameters in the set used for Markov Chain Monte Carlo fitting
#'
#' @details Takes in input list and environmental data along with names of additional parameters (vaccine efficacy
#' and reporting probabilities) and generates list of names for parameter set to use as input for fitting functions
#'
#' @param type Type of parameter set (FOI only, FOI+R0, FOI and/or R0 coefficients associated with environmental
#'   covariates); choose from "FOI","FOI+R0","FOI enviro","FOI+R0 enviro"
#' @param input_data = List of population and vaccination data for multiple regions (created using data input creation
#' code and usually loaded from RDS file)
#' @param enviro_data = Environmental data frame, containing only relevant environmental variables
#' @param extra_params = Vector of strings listing parameters besides ones determining FOI/R0 (may include vaccine
#' efficacy and/or infection/death reporting probabilities)
#' '
#' @export
#'
create_param_labels <- function(type="FOI",input_data=list(),enviro_data=NULL,extra_params=c("vacc_eff")){
  #TODO - Add assert_that functions

  n_extra=length(extra_params)

  if(type %in% c("FOI","FOI+R0")){
    regions=input_data$region_labels
    n_regions=length(regions)
    if(type=="FOI"){n_params=n_regions+n_extra} else {n_params=(2*n_regions)+n_extra}
    param_names=rep("",n_params)
    for(i in 1:n_regions){
      param_names[i]=paste("FOI_",regions[i],sep="")
      if(type=="FOI+R0"){param_names[i+n_regions]=paste("R0_",regions[i],sep="")}
    }
  } else {
    env_vars=colnames(enviro_data)[c(2:ncol(enviro_data))]
    n_env_vars=length(env_vars)
    if(type=="FOI enviro"){n_params=n_env_vars+n_extra} else {n_params=(2*n_env_vars)+n_extra}
    param_names=rep("",n_params)
    for(i in 1:n_env_vars){
      param_names[i]=paste("FOI_",env_vars[i],sep="")
      if(type=="FOI+R0 enviro"){param_names[i+n_env_vars]=paste("R0_",env_vars[i],sep="")}
    }
  }
  if(n_extra>0){param_names[(n_params-n_extra+1):n_params]=extra_params}

  return(param_names)
}
