# R file for functions relating to annual case and death data in YellowFeverDynamics package
#-------------------------------------------------------------------------------
# [TODO: CHANGE AND IF NECESSARY CREATE SEPARATE FUNCTION FOR OTHER DATA FORMATS]
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

  assert_that(is.null(model_data$rep_deaths)==FALSE,msg="Modelled data must include reported deaths")
  assert_that(is.null(obs_data$deaths)==FALSE,msg="Observed data must include numbers of deaths")
  assert_that(length(model_data$rep_deaths)==length(obs_data$deaths),
              msg="Numbers of entries in observed and modelled data must match")
  model_data$rep_deaths[model_data$rep_deaths==0]=0.1

  like_values=dnbinom(x=obs_data$deaths,mu=model_data$rep_deaths,size=rep(1,length(obs_data$deaths)),log=TRUE)
  LogLikelihood=sum(like_values,na.rm=TRUE)

  return(LogLikelihood)
}
#-------------------------------------------------------------------------------
# [TODO: CHANGE AND IF NECESSARY CREATE SEPARATE FUNCTION FOR OTHER DATA FORMATS]
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

  assert_that(is.null(model_data$rep_cases)==FALSE,msg="Modelled data must include reported cases")
  assert_that(is.null(obs_data$cases)==FALSE,msg="Observed data must include reported cases")
  assert_that(length(model_data$rep_cases)==length(obs_data$cases),
              msg="Numbers of entries in modelled and observed data must match")
  model_data$rep_cases[model_data$rep_cases==0]=0.1

  like_values=dnbinom(x=obs_data$cases,mu=model_data$rep_cases,size=rep(1,length(obs_data$cases)),log=TRUE)
  LogLikelihood=sum(like_values,na.rm=TRUE)

  return(LogLikelihood)
}
