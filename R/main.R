# R file for general functions in YellowFeverDynamics package
#------------------------------------------------
#Global variables
p_severe_inf=0.12 #Probability that an infection is severe
p_death_severe_inf=0.39 #Probability that a severe infection becomes fatal
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the NAMESPACE file.
#' @useDynLib YellowFeverDynamics, .registration = TRUE
#' @import assertthat
#' @import coda
#' @import ggplot2
#' @import graphics
#' @import grDevices
#' @import maptree
#' @import mvtnorm
#' @import qpdf
#' @import rgdal
#' @import Rmisc
#' @import sp
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
#'  If mode_start=2, use SEIRV input in list from previous run(s)
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
#' @title Generate_Dataset
#'
#' @description Generate serological, annual case/death and/or annual outbreak risk data
#'
#' @details This function is used to generate serological, annual case/death and/or annual outbreak risk data based on
#'   observed or dummy data sets; it is normally used by single_like_calc() and data_match_single() functions
#'
#' @param input_data List of population and vaccination data for multiple regions
#' @param FOI_values Values for each region of the force of infection due to spillover from sylvatic reservoir
#' @param R0_values Values for each region of the basic reproduction number for human-human transmission
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param obs_outbreak_data Outbreak Y/N data for comparison, by region and year, in format 0 = no outbreaks,
#'   1 = 1 or more outbreak(s)
#' @param vaccine_efficacy Fractional vaccine efficacy
#' @param p_rep_severe Probability of reporting of a severe but non-fatal infection
#' @param p_rep_death Probability of reporting of a fatal infection
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param n_reps Number of repeats over which to average results
#' @param dt Time increment in days to use in model (should be either 1.0 or 5.0 days)
#' '
#' @export
#'
Generate_Dataset <- function(input_data=list(),FOI_values=c(),R0_values=c(),
                             obs_sero_data=NULL,obs_case_data=NULL,obs_outbreak_data=NULL,
                             vaccine_efficacy=1.0,p_rep_severe=1.0,p_rep_death=1.0,mode_start=1,n_reps=1,dt=1.0){

  assert_that(input_data_check(input_data))
  assert_that(any(is.null(obs_sero_data)==FALSE,is.null(obs_case_data)==FALSE,is.null(obs_outbreak_data)==FALSE),
              msg="Need at least one of obs_sero_data, obs_case_data or obs_outbreak_data")
  assert_that(vaccine_efficacy >=0.0 && vaccine_efficacy <=1.0)
  if(is.null(obs_case_data)==FALSE ||is.null(obs_outbreak_data)==FALSE){
    assert_that(p_rep_severe >=0.0 && p_rep_severe <=1.0)
    assert_that(p_rep_death >=0.0 && p_rep_death <=1.0)}

  n_regions=length(input_data$region_labels)
  assert_that(length(FOI_values)==n_regions)
  assert_that(length(R0_values)==n_regions)

  if(is.null(input_data$flag_sero)){
    input_data=input_data_process(input_data,obs_sero_data,obs_case_data,obs_outbreak_data)
  }
  frac=1.0/n_reps

  #Set up data structures to take modelled data corresponding to observed data
  if(is.null(obs_sero_data)==FALSE){
    blank=rep(0,nrow(obs_sero_data))
    model_sero_data=data.frame(samples=blank,positives=blank,sero=blank)
  } else {model_sero_data=NULL}

  if(is.null(obs_case_data)==FALSE){
    model_case_values=model_death_values=rep(0,nrow(obs_case_data))
  } else {model_case_values=model_death_values=NA}

  if(is.null(obs_outbreak_data)==FALSE){
    model_outbreak_risk_values=rep(0,nrow(obs_outbreak_data))
  } else {model_outbreak_risk_values=NA}

  #Model all regions and save relevant output data
  for(n_region in 1:n_regions){

    #Get information on which observed data types are available for considered region
    flag_sero=input_data$flag_sero[n_region]
    flag_case=input_data$flag_case[n_region]
    flag_outbreak=input_data$flag_outbreak[n_region]

    #Get input data on region
    vacc_data=input_data$vacc_data[n_region,,]
    pop_data=input_data$pop_data[n_region,,]
    year_end=input_data$year_end[n_region]
    year_data_begin=input_data$year_data_begin[n_region]

    #Run model
    if(flag_sero==1){
      model_output=Full_Model_Run(FOI_values[n_region],R0_values[n_region],vacc_data,pop_data,
                                  year0=min(input_data$years_labels),mode_start,
                                  n_particles=n_reps,n_threads=n_reps,year_end,
                                  year_data_begin,vaccine_efficacy,dt=dt)
    } else {
      model_output=case_data_generate(FOI_values[n_region],R0_values[n_region],vacc_data,pop_data,
                                      year0=min(input_data$years_labels),mode_start,n_reps,
                                      year_end,year_data_begin,vaccine_efficacy,dt=dt)
    }

    #Compile outbreak/case data if needed
    if(max(flag_case,flag_outbreak)==1){
      if(flag_case==1){
        case_line_list=input_data$case_line_list[[n_region]]
        years_outbreak=obs_case_data$year[case_line_list]
        n_years_outbreak=length(case_line_list)
      } else {
        years_outbreak=obs_outbreak_data$year[input_data$outbreak_line_list[[n_region]]]
        n_years_outbreak=length(input_data$outbreak_line_list[[n_region]])
      }
      blank=array(data=rep(0,n_reps*n_years_outbreak),dim=c(n_reps,n_years_outbreak))
      annual_data=list(rep_cases=blank,rep_deaths=blank)
      for(i in 1:n_reps){
        for(n_year in 1:n_years_outbreak){
          year=years_outbreak[n_year]
          if(flag_sero==1){
            infs=sum(model_output$C[,i,model_output$year[1,]==year])
          } else {
            infs=sum(model_output$C[i,model_output$year==year])
          }
          severe_infs=rbinom(1,floor(infs),p_severe_inf)
          deaths=rbinom(1,severe_infs,p_death_severe_inf)
          annual_data$rep_deaths[i,n_year]=rbinom(1,deaths,p_rep_death)
          annual_data$rep_cases[i,n_year]=annual_data$rep_deaths[i,n_year]+rbinom(1,severe_infs-deaths,
                                                                                  p_rep_severe)
        }
      }

      if(flag_case==1){
        for(rep in 1:n_reps){
          model_case_values[case_line_list]=model_case_values[case_line_list]+annual_data$rep_cases[rep,]
          model_death_values[case_line_list]=model_death_values[case_line_list]+annual_data$rep_deaths[rep,]
        }
      }

      if(flag_outbreak==1){
        outbreak_risk=rep(0,n_years_outbreak)
        for(rep in 1:n_reps){
          for(n_year in 1:n_years_outbreak){
            if(annual_data$rep_cases[rep,n_year]>0){outbreak_risk[n_year]=outbreak_risk[n_year]+frac}
          }
        }
        for(n_year in 1:n_years_outbreak){
          if(outbreak_risk[n_year]<1.0e-4){outbreak_risk[n_year]=1.0e-4}
          if(outbreak_risk[n_year]>0.9999){outbreak_risk[n_year]=0.9999}
        }
        model_outbreak_risk_values[input_data$outbreak_line_list[[n_region]]]=outbreak_risk
      }
    }

    #Compile seroprevalence data if necessary
    if(flag_sero==1){
      sero_line_list=input_data$sero_line_list[[n_region]]
      for(i in 1:n_reps){
        sero_results=sero_calculate2(obs_sero_data[sero_line_list,],model_output,i)
        model_sero_data$samples[sero_line_list]=model_sero_data$samples[sero_line_list]+sero_results$samples
        model_sero_data$positives[sero_line_list]=model_sero_data$positives[sero_line_list]+sero_results$positives
      }
    }
    model_output<-NULL
  }

  if(is.null(obs_sero_data)==FALSE){model_sero_data$sero=model_sero_data$positives/model_sero_data$samples}
  if(is.null(obs_case_data)==FALSE){
    model_case_values=model_case_values*frac
    model_death_values=model_death_values*frac
  }

  return(list(model_sero_values=model_sero_data$sero,model_case_values=model_case_values,
              model_death_values=model_death_values,model_outbreak_risk_values=model_outbreak_risk_values))
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
  assert_that(vaccine_efficacy<=1.0 && vaccine_efficacy>=0.0)
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
#' @title param_calc_enviro
#'
#' @description Parameter calculation from environmental covariates
#'
#' @details Takes in set of coefficients of environmental covariates and calculates values of spillover
#' force of infection and reproduction number.
#'
#' @param enviro_coeffs Values of environmental coefficients
#' @param enviro_data Environmental data frame line, containing only relevant environmental covariates
#' '
#' @export
#'
param_calc_enviro <- function(enviro_coeffs=c(),enviro_data=c()){

  assert_that(all(enviro_coeffs>=0))
  n_env_vars=dim(enviro_data)[2]-1
  assert_that(length(enviro_coeffs) %in% c(n_env_vars,2*n_env_vars))
  env_vars_names=names(enviro_data)[c(2:(n_env_vars+1))]
  if(length(enviro_coeffs)==n_env_vars){
    n_type=1
    for(i in 1:n_env_vars){}
  } else {
    n_type=2
  }

  output=list(FOI=0.0,R0=0.0)
  for(i in 1:n_env_vars){
    variable=env_vars_names[i]
    output$FOI=sum(enviro_coeffs[c(1:n_env_vars)]*enviro_data[c(1:n_env_vars+1)])
    if(n_type==2){output$R0=sum(enviro_coeffs[c(1:n_env_vars)+n_env_vars]*enviro_data[c(1:n_env_vars+1)])}
  }

  return(output)
}
#-------------------------------------------------------------------------------
#' @title plot_model_output
#'
#' @description Plot output of Full_Model_Run() or Basic_Model_Run() functions
#'
#' @details Takes in the output of the Full_Model_Run() or Basic_Model_Run() functions and outputs a ggplot graph of
#'   S, R and V summed over all age groups, with error bars showing the spread of values from multiple repetitions
#'
#' @param model_output List of output data produced via Full_Model_Run() or Basic_Model_Run() functions
#' '
#' @export
#'
plot_model_output <- function(model_output=list()){

  assert_that(is.list(model_output))
  assert_that(is.null(model_output$day)==FALSE)

  values=label=values_mean=values_low=values_high=NULL
  N_age=dim(model_output$S)[1]
  n_particles=dim(model_output$S)[2]
  t_pts=dim(model_output$S)[3]
  date_values=model_output$year[1,1]+((model_output$day[1,]-model_output$day[1,1]+1)/365.0)
  blank=array(0,dim=c(n_particles,t_pts))
  S_sum=R_sum=V_sum=blank
  for(i in 1:N_age){
    S_sum=S_sum+model_output$S[i,,]
    R_sum=R_sum+model_output$R[i,,]
    V_sum=V_sum+model_output$V[i,,]
  }
  xlabels=c(model_output$year[1,1]:(model_output$year[1,t_pts]+1))
  ylabels=10^c(0:10)

  if(n_particles==1){
    data_combined=data.frame(date_values=rep(date_values,3),label=c(rep("S",t_pts),rep("R",t_pts),rep("V",t_pts)),
                             values=c(S_sum[1,],R_sum[1,],V_sum[1,]))
    plot1 <- ggplot(data=data_combined,aes(x=date_values,y=log(values),group=label))+theme_bw()
    plot1 <- plot1 + geom_line(aes(colour=label),size=1) + theme(legend.title=element_blank())
    plot1 <- plot1 + scale_x_continuous(name="",breaks=xlabels,labels=xlabels)
    plot1 <- plot1 + scale_y_continuous(name="",breaks=log(ylabels),labels=ylabels)
  } else {
    data_combined=data.frame(date_values=rep(date_values,3),label=c(rep("S",t_pts),rep("R",t_pts),rep("V",t_pts)),
                             values_mean=rep(NA,t_pts*3),values_low=rep(NA,t_pts*3),values_high=rep(NA,t_pts*3))
    for(i in 1:t_pts){
      S_CI=CI(S_sum[,i])
      R_CI=CI(R_sum[,i])
      V_CI=CI(V_sum[,i])
      data_combined$values_mean[i]=as.numeric(S_CI[2])
      data_combined$values_mean[i+t_pts]=as.numeric(R_CI[2])
      data_combined$values_mean[i+(2*t_pts)]=as.numeric(V_CI[2])
      data_combined$values_low[i]=as.numeric(S_CI[3])
      data_combined$values_low[i+t_pts]=as.numeric(R_CI[3])
      data_combined$values_low[i+(2*t_pts)]=as.numeric(V_CI[3])
      data_combined$values_high[i]=as.numeric(S_CI[1])
      data_combined$values_high[i+t_pts]=as.numeric(R_CI[1])
      data_combined$values_high[i+(2*t_pts)]=as.numeric(V_CI[1])
    }
    plot1 <- ggplot(data=data_combined,aes(x=date_values,y=log(values_mean),group=label))+theme_bw()
    plot1 <- plot1 + geom_line(aes(colour=label),size=1) + theme(legend.title=element_blank())
    plot1 <- plot1 + geom_errorbar(data=data_combined,aes(ymin=log(values_low),ymax=log(values_high),colour=label),
                                   width=(model_output$day[1,2]-model_output$day[1,1])/365.0)
    plot1 <- plot1 + scale_x_continuous(name="",breaks=xlabels,labels=xlabels)
    plot1 <- plot1 + scale_y_continuous(name="",breaks=log(ylabels),labels=ylabels)
  }

  return(plot1)
}
#-------------------------------------------------------------------------------
#' @title create_param_labels
#'
#' @description Apply names to the parameters in a set used for data matching and parameter fitting
#'
#' @details Takes in input list and environmental data along with names of additional parameters (vaccine efficacy
#' and reporting probabilities) and generates list of names for parameter set to use as input for fitting functions
#'
#' @param type Type of parameter set (FOI only, FOI+R0, FOI and/or R0 coefficients associated with environmental
#'   covariates); choose from "FOI","FOI+R0","FOI enviro","FOI+R0 enviro"
#' @param input_data List of population and vaccination data for multiple regions (created using data input creation
#' code and usually loaded from RDS file)
#' @param enviro_data Environmental data frame, containing only relevant environmental variables
#' @param extra_params Vector of strings listing parameters besides ones determining FOI/R0 (may include vaccine
#'   efficacy and/or infection/death reporting probabilities)
#'
#' @export
#'
create_param_labels <- function(type="FOI",input_data=list(),enviro_data=NULL,extra_params=c("vacc_eff")){

  assert_that(type %in% c("FOI","FOI+R0","FOI enviro","FOI+R0 enviro"))
  assert_that(input_data_check(input_data))

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
    assert_that(is.data.frame(enviro_data))
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
#' @param start_SEIRV0 SEIRV data to use as input
#' @param years_data Vector of years for which to output data
#' @param n_reps Number of repeats over which to average results
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param dt Time increment in days to use in model (should be either 1.0 or 5.0 days)
#' @param enviro_data enviro_data Data frame containing values of environmental covariates; set to NULL if not in use
#' @param R0_fixed_values Values of R0 to use if not being taken from parameter distribution
#' @param vaccine_efficacy0 Vaccine efficacy (set to NULL if being varied as a parameter)
#' @param p_rep_severe0 Probability of observation of severe infection (set to NULL if being varied as a parameter)
#' @param p_rep_death0 Probability of observation of death (set to NULL if being varied as a parameter)
#' @param m_FOI_Brazil0 Multiplier of spillover FOI for Brazil regions (set to NULL if being varied as a parameter)
#' @param flag_reporting Flag indicating whether to output number of reported severe and fatal cases
#'
#' @export
#'
total_burden_estimate <- function(type="FOI+R0 enviro",param_dist=list(),input_data=list(),start_SEIRV0=NULL,
                                  years_data=c(),n_reps=1,mode_start=1,dt=5.0,
                                  enviro_data=NULL,R0_fixed_values=NULL,vaccine_efficacy0=NULL,
                                  p_rep_severe0=NULL,p_rep_death0=NULL,m_FOI_Brazil0=NULL,flag_reporting=TRUE){

  assert_that(input_data_check(input_data))
  assert_that(all(input_data$region_labels==enviro_data$region)==TRUE)
  assert_that(min(years_data)>=input_data$years_labels[1])
  assert_that(type %in% c("FOI+R0","FOI","FOI+R0 enviro","FOI enviro"))
  assert_that(is.logical(flag_reporting))
  assert_that(all(param_dist>0.0))

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
    if(is.null(vaccine_efficacy0)){vaccine_efficacy=params$vaccine_efficacy} else {vaccine_efficacy=vaccine_efficacy0}
    if(is.null(p_rep_severe0)){p_rep_severe=params$p_rep_severe} else {p_rep_severe=p_rep_severe0}
    if(is.null(p_rep_death0)){p_rep_death=params$p_rep_death} else {p_rep_death=p_rep_death0}
    if(is.null(m_FOI_Brazil)){m_FOI_Brazil=params$m_FOI_Brazil} else {m_FOI_Brazil=m_FOI_Brazil0}

    if(type %in% c("FOI+R0 enviro","FOI enviro")){
      if(type=="FOI+R0 enviro"){enviro_coeffs=params[c(1:(2*n_env_vars))]} else {enviro_coeffs=params[c(1:n_env_vars)]}
      for(n_region in 1:n_regions){
        model_params=param_calc_enviro(enviro_coeffs,enviro_data[enviro_data$region==regions[n_region],])
        FOI_values[n_region]=model_params$FOI
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
        start_SEIRV=list(S=start_SEIRV0$S[,n_region,n_param_set],
                         E=start_SEIRV0$E[,n_region,n_param_set],I=start_SEIRV0$I[,n_region,n_param_set],
                         R=start_SEIRV0$R[,n_region,n_param_set],V=start_SEIRV0$V[,n_region,n_param_set])
      } else {start_SEIRV=NULL}
      case_data <- case_data_generate(FOI_values[n_region]*m_FOI_Brazil,R0_values[n_region],
                                  vacc_data=input_data$vacc_data[n_region,,],pop_data=input_data$pop_data[n_region,,],
                                  year0=input_data$years_labels[1],mode_start,n_reps,year_end,year_data_begin,
                                  vaccine_efficacy,start_SEIRV,dt)
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
          obs_deaths=rbinom(1,floor(deaths),p_rep_death)
          obs_cases=obs_deaths+rbinom(1,floor(cases-deaths),p_rep_severe)
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
