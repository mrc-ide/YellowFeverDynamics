# R file for functions relating to outbreak data in YellowFeverDynamics package, including annual outbreak risk and
# binary occurrence data used for MCMC parameter estimation
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
#' @param p_rep_severe probability of a severe infection being reported
#' @param p_rep_death probability of a death being reported
#' '
#' @export
#'
get_outbreak_data <- function(case_data=c(),year_data=c(),p_rep_severe=1.0,p_rep_death=1.0){
  assert_that(is.numeric(case_data)) #TODO - Improve case_data checking
  assert_that(is.numeric(year_data)) #TODO - Improve year_data checking
  assert_that(length(case_data)==length(year_data),msg="Number of entries in case data must match number of years")
  assert_that(p_rep_severe>=0 && p_rep_severe<=1.0,msg="Severe infection reporting probability must be between 0 and 1")
  assert_that(p_rep_death>=0 && p_rep_death<=1.0,msg="Fatal infection reporting probability must be between 0 and 1")

  year0=min(year_data)
  t_pts=length(year_data)
  n_years=length(table(year_data))
  dt=(n_years*365.0)/t_pts
  pt_severe_infs=pt_rep_cases=pt_deaths=pt_rep_deaths=rep(0,t_pts)
  annual_rep_cases=annual_rep_deaths=rep(0,n_years)
  for(i in 1:t_pts){
    n_year=year_data[i]-year0+1
    pt_severe_infs[i]=rbinom(1,floor(case_data[i]),p_severe_inf)
    pt_deaths[i]=rbinom(1,pt_severe_infs[i],p_death_severe_inf)
    pt_rep_deaths[i]=rbinom(1,pt_deaths[i],p_rep_death)
    pt_rep_cases[i]=rbinom(1,pt_severe_infs[i]-pt_rep_deaths[i],p_rep_severe)+pt_rep_deaths[i]
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

  outbreak_data=list(n_outbreaks=n_outbreaks,
                     outbreak_list=data.frame(size=sizes,start_day=start_days,end_day=end_days,start_year=start_years,
                                              end_year=end_years),
                     rep_pts=data.frame(pt=c(1:length(pt_rep_cases)),rep_cases=pt_rep_cases,rep_deaths=pt_rep_deaths),
                     rep_annual=data.frame(year=as.numeric(names(table(year_data))),rep_cases=annual_rep_cases,
                                           rep_deaths=annual_rep_deaths))

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
#' @param obs_data Observed annual outbreak Y/N data (0 or 1)
#' '
#' @export
#'
outbreak_risk_compare <- function(model_outbreak_risk=c(),obs_data=c()){
  assert_that(is.numeric(model_outbreak_risk))
  assert_that(is.numeric(obs_data))
  n_values=length(model_outbreak_risk)
  assert_that(n_values==length(obs_data),
              msg="Number of outbreak risk values must match number of lines in observed outbreak data")

  for(i in 1:n_values){
    if(model_outbreak_risk[i]>(1.0-1.0e-4)){model_outbreak_risk[i]=1.0-1.0e-4}
    if(model_outbreak_risk[i]<1.0e-4){model_outbreak_risk[i]=1.0e-4}
  }

  like_values=rep(NA,n_values)
  for(i in 1:n_values){
    if(obs_data[i]>0){
      like_values[i]=log(model_outbreak_risk[i])
    } else {
      like_values[i]=log(1.0-model_outbreak_risk[i])
    }
  }

  return(like_values)
}
#' #-------------------------------------------------------------------------------
#' #' @title outbreak_risk_generate
#' #'
#' #' @description Generate annual outbreak risk data
#' #'
#' #' @details Accepts epidemiological, population and case reporting parameters and model settings; generates annual
#' #' outbreak risk over a specified time period for one region. Runs basic SEIRV model (Basic_Model_Run function) for the
#' #' specified region and time period for the specified number of output datasets (n_sets), then runs case data
#' #' generation  function (case_data_generate) for a larger number of repetitions (n_reps) for each dataset, with annual
#' #' outbreak risk given by (number of repetitions where one or more cases are reported)/n_reps.
#' #'
#' #' @param FOI_spillover = Force of infection due to spillover from sylvatic reservoir
#' #' @param R0 = Reproduction number for urban spread of infection
#' #' @param vacc_data Vaccination coverage in each age group by year
#' #' @param pop_data Population in each age group by year
#' #' @param year0 First year in population/vaccination data
#' #' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#' #'  If mode_start=0, only vaccinated individuals
#' #'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#' #'  If mode_start=2, use SEIRV input in list from previous run(s)
#' #' @param n_sets Number of datasets to run pre-data period for and output
#' #' @param n_reps Number of repetitions to run per dataset
#' #' @param year_end year to run up to
#' #' @param year_data_begin year to begin saving data
#' #' @param vaccine_efficacy Proportional vaccine efficacy
#' #' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' #' @param dt Time increment in days to use in model (should be either 1.0 or 5.0 days)
#' #' @param p_rep_severe Probability of a severe infection being reported
#' #' @param p_rep_death Probability of a fatal infection being reported
#' #' '
#' #' @export
#' #'
#' outbreak_risk_generate <- function(FOI_spillover=0.0,R0=1.0,vacc_data=list(),pop_data=list(),year0=1940,
#'                               mode_start=0,n_sets=1,n_reps=1,year_end=2000,year_data_begin=1999,
#'                               vaccine_efficacy=1.0,start_SEIRV=list(),dt=1.0,p_rep_severe=1.0,p_rep_death=1.0){
#'   assert_that(p_rep_severe>=0.0 && p_rep_severe<=1.0,msg="Severe infection reporting probability must be between 0 and 1")
#'   assert_that(p_rep_death>=0.0 && p_rep_death<=1.0,msg="Fatal infection reporting probability must be between 0 and 1")
#'
#'   years=c(year_data_begin:(year_end-1))
#'   n_years=length(years)
#'   frac=1.0/n_reps
#'   outbreak_risk=array(0,dim=c(n_sets,n_years))
#'
#'   #Generate n_sets sets of starting data
#'   start_data <- Basic_Model_Run(FOI_spillover,R0,vacc_data,pop_data,year0,mode_start,n_particles=n_sets,
#'                                 n_threads=n_sets,year_end=year_data_begin+1,year_data_begin,vaccine_efficacy,
#'                                 start_SEIRV,dt)
#'
#'   n_end=dim(start_data$S)[3]
#'   for(set in 1:n_sets){
#'     start_SEIRV_set=list(S=start_data$S[,set,n_end],E=start_data$E[,set,n_end],I=start_data$I[,set,n_end],
#'                      R=start_data$R[,set,n_end],V=start_data$V[,set,n_end])
#'
#'     #Generate n_reps sets of data for period of interest for each set
#'     case_data_all <- case_data_generate(FOI_spillover,R0,vacc_data,pop_data,year0=year_data_begin,mode_start=2,
#'                                         n_reps=n_reps,year_end,year_data_begin,vaccine_efficacy,
#'                                         start_SEIRV=start_SEIRV_set,dt)
#'
#'     for(rep in 1:n_reps){
#'       for(n_year in 1:n_years){
#'         annual_cases=floor(sum(case_data_all$C[rep,case_data_all$year==years[n_year]]))
#'         severe_infs=rbinom(1,annual_cases,p_severe_inf)
#'         deaths=rbinom(1,severe_infs,p_death_severe_inf)
#'         rep_deaths=rbinom(1,deaths,p_rep_death)
#'         rep_cases=rbinom(1,severe_infs-rep_deaths,p_rep_severe)+rep_deaths
#'         if(rep_cases>0){outbreak_risk[set,n_year]=outbreak_risk[set,n_year]+frac}
#'       }
#'     }
#'   }
#'
#'   return(outbreak_risk)
#' }
