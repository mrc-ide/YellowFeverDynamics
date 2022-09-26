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
#' @param n_p = Particle to select from data
#' '
#' @export
#'
sero_calculate <- function(age_min=0,age_max=101,years=NULL,vc_factor=0,data=list(),n_p=1){

  assert_that(age_min>=0)
  assert_that(age_max>age_min)
  assert_that(is.null(years)==FALSE)
  assert_that(vc_factor>=0 && vc_factor<=1)
  assert_that(is.null(data$S)==FALSE)
  assert_that(is.null(data$E)==FALSE)
  assert_that(is.null(data$I)==FALSE)
  assert_that(is.null(data$R)==FALSE)
  assert_that(is.null(data$V)==FALSE)
  assert_that(n_p>0)
  assert_that(n_p<=dim(data$S)[2])

  ages=c((age_min+1):age_max)
  sero_values=rep(0,length(years))

  for(i in 1:length(years)){
    n_t=which(data$year %in% years[i])
    if(vc_factor==0){
      samples=data$S[ages,n_p,n_t]+data$E[ages,n_p,n_t]+data$I[ages,n_p,n_t]+data$R[ages,n_p,n_t]
      positives=data$R[ages,n_p,n_t]
      sero_values[i]=sum(positives)/sum(samples)
    } else {
     if(vc_factor==1){
       samples=data$S[ages,n_p,n_t]+data$E[ages,n_p,n_t]+data$I[ages,n_p,n_t]+data$R[ages,n_p,n_t]+data$V[ages,n_p,n_t]
       positives=data$R[ages,n_p,n_t]+data$V[ages,n_p,n_t]
       sero_values[i]=sum(positives)/sum(samples)
     } else {
       samples=data$S[ages,n_p,n_t]+data$E[ages,n_p,n_t]+data$I[ages,n_p,n_t]+data$R[ages,n_p,n_t]
       positives=data$R[ages,n_p,n_t]
       sero_values[i]=((1.0-vc_factor)*sum(positives))/sum(samples)
       samples=samples+data$V[ages,n_p,n_t]
       positives=positives+data$V[ages,n_p,n_t]
       sero_values[i]=sero_values[i]+((vc_factor*sum(positives))/sum(samples))
     }
    }
  }

  return(sero_values)
}
#-------------------------------------------------------------------------------
#' @title sero_calculate2
#'
#' @description Calculate number of "samples" and number of "positives" from modelled data for specified age range(s)
#' and year(s)
#'
#' @details Takes in information on minimum and maximum ages of desired range(s), year(s) for which to calculate
#' number of "samples" (people eligible for testing) and "positives" (people who would test positive), plus vc_factor
#' (proportion of people for whom vaccination status unknown)
#'
#' @param sero_data = Data frame containing years, minimum and maximum ages, and values of vc_factor (proportion of
#' people for whom vaccination status unknown)
#' @param model_data = Output of Basic_Model_Run or Full_Model_Run
#' @param n_p = Particle to select from model_data
#' '
#' @export
#'
sero_calculate2 <- function(sero_data=list(),model_data=list(),n_p=1){
  assert_that(is.data.frame(sero_data))
  assert_that(is.list(model_data))
  assert_that(is.null(model_data$S)==FALSE)
  assert_that(n_p>0)
  assert_that(n_p<=dim(model_data$S)[2])

  nrows=nrow(sero_data)
  output_frame=data.frame(samples=rep(NA,nrows),positives=rep(NA,nrows))

  for(i in 1:nrows){
    ages=c((sero_data$age_min[i]+1):sero_data$age_max[i])
    year=sero_data$year[i]
    vc_factor=sero_data$vc_factor[i]
    n_t=which(model_data$year==year)
    S_sum=sum(model_data$S[ages,n_p,n_t])
    E_sum=sum(model_data$E[ages,n_p,n_t])
    I_sum=sum(model_data$I[ages,n_p,n_t])
    R_sum=sum(model_data$R[ages,n_p,n_t])
    samples=S_sum+E_sum+I_sum+R_sum
    if(vc_factor>0){
      V_sum=sum(model_data$V[ages,n_p,n_t])
      if(vc_factor==1){
        samples=samples+V_sum
        positives=R_sum+V_sum
      } else {
        T_sum=samples+V_sum
        positives=((1.0-vc_factor)*R_sum)+(vc_factor*(samples/T_sum)*(R_sum+V_sum))
      }
    } else {
      positives=R_sum
    }
    output_frame$samples[i]=samples
    output_frame$positives[i]=positives
  }

  return(output_frame)
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
#'   no. samples/no. positives
#' @param n_particle = Particle to select from model_data
#' '
#' @export
#'
sero_compare <- function(model_data=list(),obs_sero_data=list(),n_particle=1){

  assert_that(is.null(model_data$S)==FALSE)
  assert_that(is.null(model_data$E)==FALSE)
  assert_that(is.null(model_data$I)==FALSE)
  assert_that(is.null(model_data$R)==FALSE)
  assert_that(is.null(model_data$V)==FALSE)
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
                                        obs_sero_data$year[i],obs_sero_data$vc_factor[i],model_data,n_particle)
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
  N_age=dimensions[1]
  n_particles=dimensions[2]
  n_pts=dimensions[3]

  like_values=rep(0,n_particles)
  for(i in 1:n_particles){like_values[i]=sero_compare(model_data,obs_sero_data,i)}

  return(like_values)
}
