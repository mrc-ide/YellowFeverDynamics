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
#' @param p_severe_inf Probability of an infection being severe
#' @param p_death_severe_inf Probability of a severe infection resulting in death
#' @param p_rep_severe probability of a severe infection being reported
#' @param p_rep_death probability of a death being reported
#' '
#' @export
#'
get_outbreak_data <- function(case_data=c(),year_data=c(),p_severe_inf = 0.12, p_death_severe_inf = 0.39, p_rep_severe=1.0,p_rep_death=1.0){
  assert_that(is.numeric(case_data)) #TODO - Improve case_data checking
  assert_that(is.numeric(year_data)) #TODO - Improve year_data checking
  assert_that(length(case_data)==length(year_data),msg="Number of entries in case data must match number of years")
  assert_that(p_severe_inf>=0.0 && p_severe_inf<=1.0,msg="Severe infection rate must be between 0-1")
  assert_that(p_death_severe_inf>=0.0 && p_death_severe_inf<=1.0,msg="Fatality rate of severe infections must be between 0-1")
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
