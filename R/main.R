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

      if(mode_start==2){
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
