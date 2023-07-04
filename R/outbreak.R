# R file for functions relating to outbreak data in YellowFeverDynamics package
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
#' @param years_data Vector of year values corresponding to case data
#' @param p_severe_inf Probability of an infection being severe
#' @param p_death_severe_inf Probability of a severe infection resulting in death
#' @param p_rep_severe Probability of a severe infection being reported
#' @param p_rep_death Probability of a death being reported
#' @param max_case_interval Maximum number of days between reported cases before outbreak is declared over
#' @param flag_additional_output TRUE/FALSE flag indicating whether to output detailed case data
#'
#' @export
#'
get_outbreak_data <- function(case_data=c(),years_data=c(),p_severe_inf = 0.12, p_death_severe_inf = 0.39, p_rep_severe=1.0,p_rep_death=1.0,
                              max_case_interval=10,flag_additional_output=FALSE){
  assert_that(is.numeric(case_data)) #TODO - Improve case_data checking
  assert_that(is.numeric(years_data)) #TODO - Improve years_data checking
  assert_that(length(case_data)==length(years_data),msg="Number of entries in case data must match number of years")
  assert_that(p_severe_inf>=0.0 && p_severe_inf<=1.0,msg="Severe infection rate must be between 0-1")
  assert_that(p_death_severe_inf>=0.0 && p_death_severe_inf<=1.0,msg="Fatality rate of severe infections must be between 0-1")
  assert_that(p_rep_severe>=0 && p_rep_severe<=1.0,msg="Severe infection reporting probability must be between 0 and 1")
  assert_that(p_rep_death>=0 && p_rep_death<=1.0,msg="Fatal infection reporting probability must be between 0 and 1")
  assert_that(is.numeric(max_case_interval))

  year0=min(years_data)
  t_pts=length(years_data)
  n_years=length(table(years_data))
  dt=(n_years*365.0)/t_pts

  pt_severe_infs=pt_rep_cases=pt_deaths=pt_rep_deaths=rep(0,t_pts)
  annual_rep_cases=annual_rep_deaths=annual_outbreaks=annual_occurrence=rep(0,n_years)
  obs_cases=obs_deaths=severe_infs=deaths=start_days=end_days=c()
  flag_outbreak=caseless_days=n_outbreaks=0

  for(i in 1:t_pts){
    n_year=years_data[i]-year0+1
    pt_severe_infs[i]=rbinom(1,floor(case_data[i]),p_severe_inf)
    pt_deaths[i]=rbinom(1,pt_severe_infs[i],p_death_severe_inf)
    pt_rep_deaths[i]=rbinom(1,pt_deaths[i],p_rep_death)
    pt_rep_cases[i]=rbinom(1,pt_severe_infs[i]-pt_rep_deaths[i],p_rep_severe)+pt_rep_deaths[i]
    annual_rep_cases[n_year]=annual_rep_cases[n_year]+pt_rep_cases[i]
    annual_rep_deaths[n_year]=annual_rep_deaths[n_year]+pt_rep_deaths[i]
    if(flag_outbreak==0){
      if(pt_rep_cases[i]>0){
        flag_outbreak=1
        n_outbreaks=n_outbreaks+1
        obs_cases=append(obs_cases,pt_rep_cases[i],after=length(obs_cases))
        obs_deaths=append(obs_deaths,pt_rep_deaths[i],after=length(obs_deaths))
        severe_infs=append(severe_infs,pt_severe_infs[i],after=length(severe_infs))
        deaths=append(deaths,pt_deaths[i],after=length(deaths))
        start_days=append(start_days,i*dt,after=length(start_days))
        caseless_days=0
      }
    } else {
      severe_infs[n_outbreaks]=severe_infs[n_outbreaks]+pt_severe_infs[i]
      deaths[n_outbreaks]=deaths[n_outbreaks]+pt_deaths[i]
      if(pt_rep_cases[i]==0){
        caseless_days=caseless_days+dt
        if(caseless_days>max_case_interval){
          flag_outbreak=0
          end_days=append(end_days,i*dt,after=length(end_days))
        }
      } else {
        caseless_days=0
        obs_cases[n_outbreaks]=obs_cases[n_outbreaks]+pt_rep_cases[i]
        obs_deaths[n_outbreaks]=obs_deaths[n_outbreaks]+pt_rep_deaths[i]
      }
    }
  }

  if(length(end_days)<length(start_days)){end_days=append(end_days,NA,after=length(end_days))}
  start_years=((start_days-(start_days %% 365))/365)+year0
  end_years=((end_days-(end_days %% 365))/365)+year0
  #outbreak_durations=end_days-start_days

  years_table=table(start_years)
  annual_outbreaks=rep(0,n_years)
  for(n_year in 1:n_years){
    if(is.na(years_table[n_year][[1]])==FALSE){annual_outbreaks[n_year]=years_table[n_year][[1]]}
    if(annual_outbreaks[n_year]>0){annual_occurrence[n_year]=1}
  }

  outbreak_data=list(annual_outbreaks=annual_outbreaks,annual_occurrence=annual_occurrence,
                     outbreak_list=data.frame(obs_cases=obs_cases,obs_deaths=obs_deaths,severe_infs=severe_infs,deaths=deaths,start_day=start_days,
                                              end_day=end_days,start_year=start_years,end_year=end_years))

  if(flag_additional_output){
    outbreak_data$rep_pts=data.frame(day=c(1:t_pts)*dt,severe_infs=pt_severe_infs,deaths=pt_deaths,rep_cases=pt_rep_cases,rep_deaths=pt_rep_deaths)
    outbreak_data$rep_annual=data.frame(year=as.numeric(names(table(years_data))),rep_cases=annual_rep_cases,rep_deaths=annual_rep_deaths)
  }

  return(outbreak_data)
}
#-------------------------------------------------------------------------------
#' @title get_outbreak_data_multi
#'
#' @description Generate outbreak data from model output over multiple repetitions
#'
#' @details [TBA]
#'
#' @param case_data Vector of cases by time point and repetition summed over age
#' @param years_data Vector of year values corresponding to case data
#' @param p_severe_inf Probability of an infection being severe
#' @param p_death_severe_inf Probability of a severe infection resulting in death
#' @param p_rep_severe Probability of a severe infection being reported
#' @param p_rep_death Probability of a death being reported
#' @param max_case_interval Maximum number of days between reported cases before outbreak is declared over
#'
#' @export
#'
get_outbreak_data_multi <- function(case_data=c(),years_data=c(),p_severe_inf = 0.12, p_death_severe_inf = 0.39, p_rep_severe=1.0,p_rep_death=1.0,
                              max_case_interval=10){
  assert_that(length(dim(case_data))==2)
  #TODO - Additional assert_that functions

  n_reps=dim(case_data)[1]
  n_years=length(unique(years_data))
  outbreak_list_all=data.frame(rep=0,obs_cases=0,obs_deaths=0,severe_infs=0,deaths=0,start_day=0,end_day=0,start_year=0,end_year=0)
  annual_outbreaks_all=array(NA,dim=c(n_reps,n_years))
  annual_outbreak_probs=rep(0,n_years)

  for(rep in 1:n_reps){
    outbreak_data=get_outbreak_data(case_data=case_data[rep,],years_data,p_severe_inf,p_death_severe_inf,p_rep_severe,p_rep_death,max_case_interval,
                                    flag_additional_output=FALSE)
    outbreak_list_all=rbind(outbreak_list_all,cbind(rep=rep(rep,nrow(outbreak_data$outbreak_list)),outbreak_data$outbreak_list))
    annual_outbreaks_all[rep,]=outbreak_data$annual_outbreaks
    annual_outbreak_probs=annual_outbreak_probs+outbreak_data$annual_occurrence
  }
  annual_outbreak_probs=annual_outbreak_probs/n_reps
  outbreak_list_all=outbreak_list_all[c(2:nrow(outbreak_list_all)),]

  return(list(outbreak_list_all=outbreak_list_all,annual_outbreaks_all=annual_outbreaks_all,annual_outbreak_probs=annual_outbreak_probs))
}
