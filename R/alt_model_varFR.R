#-------------------------------------------------------------------------------
#' @title Full_Model_Run_VarFR
#'
#' @description Run alternate full version of SEIRV model with annually varying FOI_spillover and R0
#'
#' @details Accepts epidemiological + population parameters and model settings; runs full version of SEIRV model
#' for one region over a specified time period for a number of particles/threads and outputs time-dependent SEIRV
#' values, infection numbers and total force of infection values.
#'
#' @param FOI_spillover Vector of annual values of force of infection due to spillover from sylvatic reservoir
#' @param R0 Vector of annual values of basic reproduction number for urban spread of infection
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
Full_Model_Run_VarFR <- function(FOI_spillover=0.0,R0=1.0,vacc_data=list(),pop_data=list(),year0=1940,mode_start=0,
                           n_particles=1,n_threads=1,year_end=2000,year_data_begin=1999,vaccine_efficacy=1.0,
                           start_SEIRV=list(),dt=1.0) {

  assert_that(n_particles %in% c(1:20),msg="Number of particles must be an integer between 1 and 20")
  assert_that(n_threads<=n_particles,msg="Number of threads must be equal to or less than number of particles")
  assert_that(n_threads>0,msg="Number of threads must be between 1 and number of particles")

  x <- FullModelODVarFR$new(pars=parameter_setup_varFR(FOI_spillover,R0,vacc_data,pop_data,year0,mode_start,year_end,
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
# R file for functions relating to annual case and death data in YellowFeverDynamics package
#-------------------------------------------------------------------------------
#' @title case_data_generate_FR
#'
#' @description Take in single set of population data and model parameters, output infection data only
#'
#' @details Accepts epidemiological + population parameters and model settings; runs full version of SEIRV model
#' for one region over a specified time period for a number of particles/threads and outputs infection numbers at each
#' time increment only; optimized for running a large number of repetitions
#'
#' @param FOI_spillover = Force of infection due to spillover from sylvatic reservoir (variable by year)
#' @param R0 = Basic reproduction number for urban spread of infection (variable by year)
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
case_data_generate_varFR <- function(FOI_spillover=c(),R0=c(),vacc_data=list(),pop_data=list(),year0=1940,
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

  pars=parameter_setup_varFR(FOI_spillover,R0,vacc_data,pop_data,year0,mode_start,year_end,
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
    x <- FullModelODVarFR$new(pars=pars,time = 1,n_particles = n_particles,n_threads = n_threads)

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
#' @title parameter_setup_varFR
#'
#' @description Set up parameters to input into model (alternate version for use with model version with annually
#'   varying FOI_spillover and R0)
#'
#' @details Takes in multiple inputs, outputs list for use by odin.dust SEIRV model versions.
#'
#' @param FOI_spillover Vector of annual values of force of infection due to spillover from sylvatic reservoir
#' @param R0 Vector of annual values of basic reproduction number for urban spread of infection
#' @param vacc_data Vaccination coverage in each age group by year
#' @param pop_data Population in each age group by year
#' @param year0 First year in population/vaccination data
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param year_end year to run up to
#' @param year_data_begin year to begin saving data
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' @param dt Time increment in days to use in model (should be either 1.0 or 5.0 days)
#' '
#' @export
#'
parameter_setup_varFR <- function(FOI_spillover=0.0,R0=1.0,vacc_data=list(),pop_data=list(),year0=1940,mode_start=0,
                            year_end=2000,year_data_begin=1999,vaccine_efficacy=1.0,start_SEIRV=list(),dt=1.0){

  assert_that(length(pop_data[,1])>1) #TODO - msg
  assert_that(length(pop_data[1,])>1) #TODO - msg
  n_years=length(pop_data[,1])-1
  assert_that(length(FOI_spillover)==n_years) #TODO - msg
  assert_that(length(R0)==n_years) #TODO - msg
  N_age=length(pop_data[1,])
  assert_that(length(vacc_data[,1])==n_years+1,msg="Population and vaccination data must be for same time periods")
  assert_that(length(vacc_data[1,])==N_age,msg="Number of age groups in population and vaccination data must match")
  assert_that(mode_start %in% c(0,1,2),msg="mode_start must have value 0, 1 or 2")
  assert_that(vaccine_efficacy<=1.0 && vaccine_efficacy>=0.0,msg="Vaccine efficacy must be between 0 and 1")
  if(mode_start==2){assert_that(is.null(start_SEIRV$S)==FALSE,msg="When mode_start=2, start_SEIRV data is required")}
  assert_that(year_data_begin>=year0,msg="year_data_begin must be greater than or equal to year0")
  assert_that(year_data_begin<year_end,msg="year_data_begin must be less than year_end")
  assert_that(year_end-year0<=n_years,msg="Period year0->year_end must lie within population data")
  vacc_initial=vacc_data[1,]
  assert_that(dt %in% c(1,5),msg="dt must have value 1 or 5 days")
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
    if(R0[1]>1.0){
      herd_immunity=1.0-(1.0/R0[1])
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
    Vac0=P0*vacc_initial
  }

  return(list(FOI_spillover=FOI_spillover,R0=R0,vacc_rate_annual=vacc_rates,
              Cas0=Cas0,Exp0=Exp0,Inf0=Inf0,N_age=N_age,Rec0=Rec0,Sus0=Sus0,Vac0=Vac0,
              dP1_all=dP1_all,dP2_all=dP2_all,n_years=n_years,year0=year0,vaccine_efficacy=vaccine_efficacy,dt=dt))
}

#-------------------------------------------------------------------------------
#' @title total_burden_estimate_varFR
#'
#' @description Function to calculate annual yellow fever burden across multiple regions based on derived parameters
#'   (variation for annually varying FOI_spillover and R0)
#'
#' @details Function to take in parameter sets derived from MCMC fitting and use to calculate annual total and reported
#' case and death numbers for multiple regions to compare with external data
#'
#' @param FOI_spillover_values #Array of values of spillover force of infection by region (dim1) and year (dim2)
#' @param R0_values # #Array of values of R0 by region (dim1) and year (dim2)
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
#' @param vaccine_efficacy Vaccine efficacy (set to NULL if being varied as a parameter)
#' @param p_rep_severe Probability of observation of severe infection (set to NULL if being varied as a parameter)
#' @param p_rep_death Probability of observation of death (set to NULL if being varied as a parameter)
#'
#' @export
#'
total_burden_estimate_varFR <- function(FOI_spillover_values=list(),R0_values=list(),input_data=list(),
                                        start_SEIRV=NULL,years_data=c(),n_reps=1,mode_start=1,flag_reporting=FALSE,
                                        dt=5.0,vaccine_efficacy=NULL,p_rep_severe=NULL,p_rep_death=NULL){

  #TODO - Add additional assert_that checks
  assert_that(input_data_check(input_data),
              msg="Input data must be in standard format (see https://mrc-ide.github.io/YellowFeverDynamics/articles/CGuideAInputs.html )")
  assert_that(min(years_data)>=input_data$years_labels[1]) #TODO - msg
  assert_that(is.logical(flag_reporting))

  n_years=length(years_data)
  year_data_begin=years_data[1]
  year_end=max(years_data)+1
  regions=input_data$region_labels
  n_regions=length(regions)
  case_ar1=death_ar1=array(NA,dim=c(n_years,n_regions,n_reps))
  case_ar2=death_ar2=array(NA,dim=c(n_years,n_regions))
  case_ar3=death_ar3=array(NA,dim=c(n_years))

  for(n_region in 1:n_regions){

    start_SEIRV_region=list(S=start_SEIRV$S[,n_region,1],
                            E=start_SEIRV$E[,n_region,1],I=start_SEIRV$I[,n_region,1],
                            R=start_SEIRV$R[,n_region,1],V=start_SEIRV$V[,n_region,1])
    case_data <- case_data_generate_varFR(FOI_spillover_values[n_region,],R0_values[n_region,],
                                          vacc_data=input_data$vacc_data[n_region,,],
                                          pop_data=input_data$pop_data[n_region,,],year0=input_data$years_labels[1],
                                          mode_start,n_reps,year_end,year_data_begin,vaccine_efficacy,
                                          start_SEIRV_region,dt)
    for(n_year in 1:n_years){
      for(rep in 1:n_reps){
        infs=floor(sum(case_data$C[rep,case_data$year==years_data[n_year]]))
        severe_infs=rbinom(1,infs,p_severe_inf)
        deaths=rbinom(1,severe_infs,p_death_severe_inf)
        case_ar1[n_year,n_region,rep]=severe_infs
        death_ar1[n_year,n_region,rep]=deaths
      }
      case_ar2[n_year,n_region]=sum(case_ar1[n_year,n_region,])/n_reps
      death_ar2[n_year,n_region]=sum(death_ar1[n_year,n_region,])/n_reps
    }
  }

  for(n_year in 1:n_years){
    case_ar3[n_year]=sum(case_ar2[n_year,])
    death_ar3[n_year]=sum(death_ar2[n_year,])
  }

  if(flag_reporting){
    obs_case_ar2=obs_death_ar2=array(NA,dim=c(n_years,n_regions))
    obs_case_ar3=obs_death_ar3=array(NA,dim=c(n_years))

    for(n_region in 1:n_regions){
      for(n_year in 1:n_years){
        cases=case_ar2[n_year,n_region]
        deaths=death_ar2[n_year,n_region]
        obs_deaths=rbinom(1,floor(deaths),p_rep_death)
        obs_cases=obs_deaths+rbinom(1,floor(cases-deaths),p_rep_severe)
        obs_case_ar2[n_year,n_region]=obs_cases
        obs_death_ar2[n_year,n_region]=obs_deaths
      }
    }
    for(n_year in 1:n_years){
      obs_case_ar3[n_year]=sum(obs_case_ar2[n_year,])
      obs_death_ar3[n_year]=sum(obs_death_ar2[n_year,])
    }
    plot_frame1=data.frame(year=rep(years_data,n_regions),
                           region=rep(sort(rep(regions,n_years))),
                           cases=as.vector(case_ar2),deaths=as.vector(death_ar2),
                           obs_cases=as.vector(obs_case_ar2),obs_deaths=as.vector(obs_death_ar2))
    plot_frame2=data.frame(year=rep(years_data),
                           cases=as.vector(case_ar3),deaths=as.vector(death_ar3),
                           obs_cases=as.vector(obs_case_ar3),obs_deaths=as.vector(obs_death_ar3))

  } else {
    plot_frame1=data.frame(year=rep(years_data,n_regions),
                           region=rep(sort(rep(regions,n_years))),
                           cases=as.vector(case_ar2),deaths=as.vector(death_ar2))
    plot_frame2=data.frame(year=rep(years_data),
                           cases=as.vector(case_ar3),deaths=as.vector(death_ar3))
  }


  return(list(by_region=plot_frame1,all=plot_frame2))
}
