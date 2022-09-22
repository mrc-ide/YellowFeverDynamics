# R file for functions relating to annual case and outbreak data in YellowFeverDynamics package
#-------------------------------------------------------------------------------
#' @title get_outbreak_data
#'
#' @description Generate outbreak data from model output
#'
#' @details Takes in case (infection) data output by Full_Model_OD() or case_data_generate(), uses severe and fatal
#' infection rates and reporting probabilities to calculate the number of reported cases, and outputs a list of
#' reported outbreaks. An outbreak is assumed to be reported when one or more cases is reported, and continues until
#' no new cases have been reported for 10 days.
#'
#' @param case_data Vector of cases by time point summed over age
#' @param year_data Vector of year values corresponding to case data
#' @param p_rep_mild probability of a mild infection being reported
#' @param p_rep_severe probability of a severe infection being reported
#' @param p_rep_death probability of a death being reported
#' '
#' @export
#'
get_outbreak_data <- function(case_data=c(),year_data=c(),p_rep_mild=0.0,p_rep_severe=1.0,p_rep_death=1.0){
  assert_that(is.numeric(case_data))
  assert_that(is.numeric(year_data))
  assert_that(length(case_data)==length(year_data))
  assert_that(p_rep_mild>=0 && p_rep_mild<=1.0)
  assert_that(p_rep_severe>=0 && p_rep_severe<=1.0)
  assert_that(p_rep_death>=0 && p_rep_death<=1.0)

  year0=min(year_data)
  t_pts=length(year_data)
  n_years=length(table(year_data))
  dt=(n_years*365.0)/t_pts
  pt_severe_infs=pt_mild_infs=pt_rep_cases=pt_deaths=pt_rep_deaths=rep(0,t_pts)
  annual_rep_cases=annual_rep_deaths=rep(0,n_years)
  for(i in 1:t_pts){
    n_year=year_data[i]-year0+1
    pt_severe_infs[i]=rbinom(1,case_data[i],p_severe_inf)
    pt_mild_infs[i]=case_data[i]-pt_severe_infs[i]
    pt_deaths[i]=rbinom(1,pt_severe_infs[i],p_death_severe_inf)
    pt_rep_deaths[i]=rbinom(1,pt_deaths[i],p_rep_death)
    pt_rep_cases[i]=rbinom(1,pt_mild_infs,p_rep_mild)+rbinom(1,pt_severe_infs[i]-pt_rep_deaths[i],
                                                                    p_rep_severe)+pt_rep_deaths[i]
    annual_rep_cases[n_year]=annual_rep_cases[n_year]+pt_rep_cases[i]
    annual_rep_deaths[n_year]=annual_rep_deaths[n_year]+pt_rep_deaths[i]
  }

  n_outbreaks=0
  sizes=start_days=end_days=c()
  flag_outbreak=0
  caseless_days=0
  for(i in 1:t_pts){
    if(flag_outbreak==0){
      if(pt_rep_cases[i]>0){
        flag_outbreak=1
        n_outbreaks=n_outbreaks+1
        sizes=append(sizes,pt_rep_cases[i],after=length(sizes))
        start_days=append(start_days,i*dt,after=length(start_days))
        caseless_days=0
      }
    } else {
      if(pt_rep_cases[i]==0){
        caseless_days=caseless_days+dt
        if(caseless_days==10){
          flag_outbreak=0
          end_days=append(end_days,i*dt,after=length(end_days))
        }
      } else {
        caseless_days=0
        sizes[n_outbreaks]=sizes[n_outbreaks]+pt_rep_cases[i]
      }
    }
  }

  if(length(end_days)<length(start_days)){end_days=append(end_days,i,after=length(end_days))}
  start_years=((start_days-(start_days %% 365))/365)+year0
  end_years=((end_days-(end_days %% 365))/365)+year0

  outbreak_data=list(n_outbreaks=n_outbreaks,sizes=sizes,
                     start_days=start_days,end_days=end_days,start_years=start_years,end_years=end_years,
                     pt_rep_cases=pt_rep_cases,pt_rep_deaths=pt_rep_deaths,
                     annual_rep_cases=annual_rep_cases,annual_rep_deaths=annual_rep_deaths)

  return(outbreak_data)
}
#-------------------------------------------------------------------------------
#' @title outbreak_risk_compare
#'
#' @description Calculate likelihood of observed outbreak data based on modelled outbreak risk
#'
#' @details Takes in modelled outbreak risk data (probabilities of one or more outbreaks being reported in given
#' years) and compared with observed data on whether 1 or more outbreaks were reported in given years, calculating
#' logarithmic likelihood of observing the latter given the former.
#'
#' @param model_outbreak_risk Modelled annual outbreak risk (0-1)
#' @param obs_data Observed annual outbreak Y/N data
#' '
#' @export
#'
outbreak_risk_compare <- function(model_outbreak_risk=list(),obs_data=list()){
  assert_that(is.list(model_outbreak_risk))
  assert_that(is.list(obs_data))
  assert_that(length(model_outbreak_risk)==length(obs_data))

  for(i in 1:length(model_outbreak_risk)){
    if(model_outbreak_risk[i]>(1.0-1.0e-4)){model_outbreak_risk[i]=1.0-1.0e-4}
    if(model_outbreak_risk[i]<1.0e-4){model_outbreak_risk[i]=1.0e-4}
  }

  LogLikelihood=0
  for(i in 1:length(model_outbreak_risk)){
    if(obs_data[i]>0){
      LogLikelihood=LogLikelihood+log(model_outbreak_risk[i])
    } else {
      LogLikelihood=LogLikelihood+log(1.0-model_outbreak_risk[i])
    }
  }

  return(LogLikelihood)
}

#-------------------------------------------------------------------------------
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
  assert_that(length(model_data$rep_deaths[1,])==length(obs_data$deaths))
  model_data$rep_deaths[model_data$rep_deaths==0]=0.1

  like_values=dnbinom(x=obs_data$deaths,mu=model_data$rep_deaths,size=rep(1,length(obs_data$deaths)),log=TRUE)
  #like_values[like_values==-Inf]=NA
  LogLikelihood=sum(like_values,na.rm=TRUE)

  return(LogLikelihood)
}
#-------------------------------------------------------------------------------
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
  assert_that(length(model_data$rep_cases[1,])==length(obs_data$cases))
  model_data$rep_cases[model_data$rep_cases==0]=0.1

  like_values=dnbinom(x=obs_data$cases,mu=model_data$rep_cases,size=rep(1,length(obs_data$cases)),log=TRUE)
  LogLikelihood=sum(like_values,na.rm=TRUE)

  return(LogLikelihood)
}
#-------------------------------------------------------------------------------
#' @title case_data_generate
#'
#' @description Take in single set of population data and model parameters, output infection data only
#'
#' @details Accepts epidemiological + population parameters and model settings; runs full version of SEIRV model
#' for one region over a specified time period for a number of particles/threads and infection numbers at each time
#' increment only; optimized for running a large number of repetitions
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
    x <- FullModelOD$new(pars=pars,step = 1,n_particles = n_particles,n_threads = n_threads)

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
#' @title outbreak_risk_generate
#'
#' @description Generate annual outbreak risk data
#'
#' @details This function is used to calculate the annual risk of one or more outbreaks in a region over one or more
#'   years, based on epidemiological parameters
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
#' @param n_sets Number of datasets to run pre-data period for and output
#' @param n_reps Number of repetitions to run per dataset
#' @param year_end year to run up to
#' @param year_data_begin year to begin saving data
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' @param dt Time increment in days to use in model (should be either 1.0 or 5.0 days)
#' @param p_rep_severe Probability of a severe case being reported
#' @param p_rep_death Probability of a fatal case being reported
#' '
#' @export
#'
outbreak_risk_generate <- function(FOI_spillover=0.0,R0=1.0,vacc_data=list(),pop_data=list(),year0=1940,
                              mode_start=0,n_sets=1,n_reps=1,year_end=2000,year_data_begin=1999,
                              vaccine_efficacy=1.0,start_SEIRV=list(),dt=1.0,p_rep_severe=1.0,p_rep_death=1.0){
  assert_that(p_rep_severe>=0.0 && p_rep_severe<=1.0)
  assert_that(p_rep_death>=0.0 && p_rep_death<=1.0)

  years=c(year_data_begin:(year_end-1))
  n_years=length(years)
  frac=1.0/n_reps
  outbreak_risk=array(0,dim=c(n_sets,n_years))

  #Generate n_sets sets of starting data
  start_data <- Basic_Model_Run(FOI_spillover,R0,vacc_data,pop_data,year0,mode_start,n_particles=n_sets,
                                n_threads=n_sets,year_end=year_data_begin+1,year_data_begin,vaccine_efficacy,
                                start_SEIRV,dt)

  n_end=dim(start_data$S)[3]
  for(set in 1:n_sets){
    start_SEIRV_set=list(S=start_data$S[,set,n_end],E=start_data$E[,set,n_end],I=start_data$I[,set,n_end],
                     R=start_data$R[,set,n_end],V=start_data$V[,set,n_end])

    #Generate n_reps sets of data for period of interest for each set
    case_data_all <- case_data_generate(FOI_spillover,R0,vacc_data,pop_data,year0=year_data_begin,mode_start=2,
                                        n_reps=n_reps,year_end,year_data_begin,vaccine_efficacy,
                                        start_SEIRV=start_SEIRV_set,dt)

    for(rep in 1:n_reps){
      for(n_year in 1:n_years){
        annual_cases=floor(sum(case_data_all$C[rep,case_data_all$year==years[n_year]]))
        severe_infs=rbinom(1,annual_cases,p_severe_inf)
        deaths=rbinom(1,severe_infs,p_death_severe_inf)
        rep_deaths=rbinom(1,deaths,p_rep_death)
        rep_cases=rbinom(1,severe_infs-rep_deaths,p_rep_severe)+rep_deaths
        if(rep_cases>0){outbreak_risk[set,n_year]=outbreak_risk[set,n_year]+frac}
      }
    }
  }

  return(outbreak_risk)
}
