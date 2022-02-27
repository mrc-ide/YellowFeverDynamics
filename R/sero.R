# R file for functions relating to serological data in YellowFeverDynamics package
#-------------------------------------------------------------------------------
#' @title sero_calculate
#'
#' @description Calculate seroprevalence in unvaccinated people from modelled data for one or more years and one age
#' range
#'
#' @details Takes in information on minimum and maximum ages of desired range, year(s) for which to calculate
#' seroprevalence, factor representing proportion of patients with unknown vaccine status, and SEIRV model output
#' data, and calculates seroprevalence in unvaccinated people in specified age range for specified year(s).
#'
#' @param age_min = Minimum age of age group
#' @param age_max = Maximum age of age group
#' @param years = Years for which to calculate average annual seroprevalence
#' @param vc_factor = Proportion of patients tested for whom vaccine status unknown
#' @param data = Output of Basic_Model_Run or Full_Model_Run
#' '
#' @export
#'
sero_calculate <- function(age_min=0,age_max=101,years=NULL,vc_factor=0,data=list()){

  assert_that(age_min>=0)
  assert_that(age_max>age_min)
  assert_that(is.null(years)==FALSE)
  assert_that(vc_factor>=0 && vc_factor<=1)
  assert_that(is.null(data$S)==FALSE)
  assert_that(is.null(data$E)==FALSE)
  assert_that(is.null(data$I)==FALSE)
  assert_that(is.null(data$R)==FALSE)
  assert_that(is.null(data$V)==FALSE)
  ages=c((age_min+1):age_max)
  sero_values=rep(0,length(years))

  for(i in 1:length(years)){
    days=which(data$year %in% years[i])
    if(vc_factor==0){
      samples=data$S[days,ages]+data$E[days,ages]+data$I[days,ages]+data$R[days,ages]
      positives=data$R[days,ages]
      sero_values[i]=sum(positives)/sum(samples)
    } else {
     if(vc_factor==1){
       samples=data$S[days,ages]+data$E[days,ages]+data$I[days,ages]+data$R[days,ages]+data$V[days,ages]
       positives=data$R[days,ages]+data$V[days,ages]
       sero_values[i]=sum(positives)/sum(samples)
     } else {
       samples=data$S[days,ages]+data$E[days,ages]+data$I[days,ages]+data$R[days,ages]
       positives=data$R[days,ages]
       sero_values[i]=((1.0-vc_factor)*sum(positives))/sum(samples)
       samples=samples+data$V[days,ages]
       positives=positives+data$V[days,ages]
       sero_values[i]=sero_values[i]+((vc_factor*sum(positives))/sum(samples))
     }
    }
  }

  return(sero_values)
}
#-------------------------------------------------------------------------------
#' @title sero_compare
#'
#' @description Take model results, calculate seroprevalence for comparison with observed seroprevalence and
#' calculate likelihood (single region, multiple years/age ranges)
#'
#' @details Takes in SEIRV model output data and observed seroprevalence data, calculates seroprevalence from modelled
#' data, and compares modelled and observed data, calculating logarithmic likelihood of observing the latter given the
#' former, using a binomial formula.
#'
#' @param model_data = Output of Basic_Model_Run_OD or Full_Model_Run_OD
#' @param obs_sero_data = Seroprevalence data for comparison, by year and age group, in format
#' no. samples/no. positives
#' '
#' @export
#'
sero_compare <- function(model_data=list(),obs_sero_data=list()){

  assert_that(is.null(model_data$S)==FALSE)
  assert_that(is.null(obs_sero_data$year)==FALSE)
  assert_that(is.null(obs_sero_data$age_min)==FALSE)
  assert_that(is.null(obs_sero_data$age_max)==FALSE)
  assert_that(is.null(obs_sero_data$samples)==FALSE)
  assert_that(is.null(obs_sero_data$positives)==FALSE)
  assert_that(is.null(obs_sero_data$vc_factor)==FALSE)

  n_lines=length(obs_sero_data$year)
  model_sero_data=rep(0,n_lines)

  for(i in 1:n_lines){
    model_sero_data[i]=sero_calculate(obs_sero_data$age_min[i],obs_sero_data$age_max[i],
                                        obs_sero_data$year[i],obs_sero_data$vc_factor[i],model_data)
  }
  model_sero_data[is.na(model_sero_data)]=0
  model_sero_data[is.infinite(model_sero_data)]=0

  LogLikelihood = sum(lgamma(obs_sero_data$samples+1)-lgamma(obs_sero_data$positives+1)
                      -lgamma(obs_sero_data$samples-obs_sero_data$positives+1) +
                        obs_sero_data$positives*log(model_sero_data) +
                        (obs_sero_data$samples-obs_sero_data$positives)*log(1.0-model_sero_data))

  return(LogLikelihood)
}
#-------------------------------------------------------------------------------
#' @title sero_compare_multiparticle
#'
#' @description Take model results from multi-particle run using OD version of model, calculate
#' seroprevalence for comparison with observed seroprevalence  and  calculate likelihood value (single region,
#' multiple years/age ranges) for each particle; return as vector of values
#'
#' @details Calculates logarithmic likelihood of observing observed seroprevalence given modelled data, similar to
#' sero_compare(), but using modelled data for multiple particles.
#'
#' @param model_data = Output of Basic_Model_Run or Full_Model_Run
#' @param obs_sero_data = Seroprevalence data for comparison, by year and age group, in format no. samples/no.
#'   positives
#' '
#' @export
#'
sero_compare_multiparticle <- function(model_data=list(),obs_sero_data=list()){

  dimensions=dim(model_data$S)
  if(length(dimensions)==3){
    N_age=dimensions[1]
    n_particles=dimensions[2]
    n_pts=dimensions[3]
  } else {
    N_age=dimensions[1]
    n_particles=1
    n_pts=dimensions[2]
  }

  like_values=rep(0,n_particles)

  for(i in 1:n_particles){
    data_1set=list(day=model_data$day[i,],year=model_data$year[i,],S=t(model_data$S[,i,]),E=t(model_data$E[,i,]),
                   I=t(model_data$I[,i,]),R=t(model_data$R[,i,]),V=t(model_data$V[,i,]))
    like_values[i]=sero_compare(data_1set,obs_sero_data)
  }

  return(like_values)
}
