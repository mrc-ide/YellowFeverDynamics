# R file for functions used for Markov Chain Monte Carlo fitting (and preliminary maximum-likelihood fitting) in
# YellowFeverDynamics package
#-------------------------------------------------------------------------------
#' @title MCMC
#'
#' @description Combined MCMC Multi-Region - series of MCMC steps for one or more regions
#'
#' @details TBA
#'
#' @param type Type of parameter set (FOI only, FOI+R0, FOI and/or R0 coefficients associated with environmental
#'   covariates); choose from "FOI","FOI+R0","FOI enviro","FOI+R0 enviro"
#' @param pars_ini Initial log values of parameters
#' @param pars_min Lower limits of parameter values if specified
#' @param pars_max Upper limits of parameter values if specified
#' @param input_data List of population and vaccination data for multiple regions (created using data input creation
#'   code and usually loaded from RDS file)
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no.
#'   cases/no. deaths
#' @param obs_outbreak_data Outbreak Y/N data for comparison, by region and year, in format 0 = no outbreaks,
#'   1 = 1 or more outbreak(s)
#' @param n_reps Number of times to repeat calculations to get average likelihood at each step
#' @param Niter Total number of steps to run
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start = 0, only vaccinated individuals
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity
#' @param prior_type Text indicating which type of calculation to use for prior probability
#'  If prior_type = "zero", prior probability is always zero
#'  If prior_type = "flat", prior probability is zero if FOI/R0 in designated ranges, -Inf otherwise
#'  If prior_type = "exp", prior probability is given by dexp calculation on FOI/R0 values
#'  If prior_type = "norm", prior probability is given by dnorm calculation on parameter values
#' @param dt time increment in days (must be 1 or 5)
#' @param enviro_data Values of environmental variables (if in use)
#' @param R0_fixed_values Values of R0 to use if not being fitted
#' @param vaccine_efficacy Vaccine efficacy (set to NULL if being varied as a parameter)
#' @param p_obs_severe Probability of observation of severe infection (set to NULL if being varied as a parameter)
#' @param p_obs_death Probability of observation of death (set to NULL if being varied as a parameter)
#' @param filename_prefix Prefix of names for output files
MCMC <- function(type=NULL,pars_ini=c(),pars_min=NULL,pars_max=NULL,input_data=list(),
                          obs_sero_data=NULL,obs_case_data=NULL,obs_outbreak_data=NULL,n_reps=1,Niter=1,
                          mode_start=0,prior_type="zero",dt=1.0,enviro_data=NULL,R0_fixed_values=NULL,
                          vaccine_efficacy=NULL,p_obs_severe=NULL,p_obs_death=NULL,filename_prefix="Chain"){

  regions=names(table(input_data$region_labels))
  n_regions=length(regions)
  n_params=length(pars_ini)
  ###########################################
  checks<-mcmc_checks(pars_ini,n_params,type,prior_type,n_regions,enviro_data,R0_fixed_values,
                      vaccine_efficacy,p_obs_severe,p_obs_death)
  ###########################################

  ### find posterior probability at start ###
  out = MCMC_step(type,param=pars_ini,pars_min,pars_max,input_data,obs_sero_data,obs_case_data,
                           obs_outbreak_data,chain_cov=1,adapt=0,like_current=-Inf,n_reps,mode_start,prior_type,dt,
                           enviro_data,R0_fixed_values,vaccine_efficacy,p_obs_severe,p_obs_death)

  #mcmc setup
  chain=chain_prop=posterior_current=posterior_prop=flag_accept=chain_cov_all=NULL
  burnin = min(2*n_params, Niter)
  fileIndex = 0
  chain_cov=1

  for (iter in 1:Niter){
    #save current step
    param = out$param
    param_prop=out$param_prop
    like_current = out$like_current
    like_prop=out$like_prop
    accept = out$accept
    chain = rbind(chain, param)
    chain_prop=rbind(chain_prop,param_prop)
    posterior_current=rbind(posterior_current,like_current)
    posterior_prop=rbind(posterior_prop,like_prop)
    flag_accept = rbind(flag_accept, accept)
    chain_cov_all = rbind(chain_cov_all,max(chain_cov))

    if(iter==1){ #Set output headings
      colnames(chain) = names(pars_ini)
      colnames(chain_prop) = names(pars_ini)
      for(i in 1:n_params){
        colnames(chain_prop)[i]=paste("Test_",colnames(chain_prop)[i],sep="")
      }
      colnames(posterior_current) = "posterior_current"
      colnames(posterior_prop) = "posterior_prop"
      colnames(flag_accept) = "flag_accept"
      colnames(chain_cov_all) = "chain_cov_all"
    }

    if (iter %% 10 == 0){
      if (iter %% 10000 == 0){fileIndex  = iter/10000}

      filename=paste(filename_prefix,fileIndex,".csv",sep="")
      lines=min((fileIndex * 10000+1),iter):iter
      cat("\nIteration ",iter,sep="")
      data_out<-cbind(posterior_current,posterior_prop,exp(chain),flag_accept,exp(chain_prop),chain_cov_all)[lines,]
      write.csv(data_out,filename,row.names=FALSE)
    }

    #adapt?
    if (iter>burnin & runif(1)<0.9){ #adapt
      adapt = 1
      chain_cov  = cov(chain[max(nrow(chain)-10000, 1):nrow(chain),])
    } else {
      adapt = 0
      chain_cov = 1
    }

    #new step
    out = MCMC_step(type,param,pars_min,pars_max,input_data,obs_sero_data,obs_case_data,
                             obs_outbreak_data,chain_cov,adapt,like_current,n_reps,mode_start,prior_type,dt,
                             enviro_data,R0_fixed_values,vaccine_efficacy,p_obs_severe,p_obs_death)
  }

  param_out=exp(out$param)
  names(param_out)=names(pars_ini)

  return(param_out)
}
#-------------------------------------------------------------------------------
#' @title MCMC_step
#'
#' @description Single MCMC step - one or more regions
#'
#' @details TBA
#'
#' @param type Type of parameter set (FOI only, FOI+R0, FOI and/or R0 coefficients associated with environmental
#'   covariates); choose from "FOI","FOI+R0","FOI enviro","FOI+R0 enviro"
#' @param param Log values of parameters
#' @param pars_min Lower limits of parameter values if specified
#' @param pars_max Upper limits of parameter values if specified
#' @param input_data List of population and vaccination data for multiple regions (created using data input creation
#'   code and usually loaded from RDS file)
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param obs_outbreak_data Outbreak Y/N data for comparison, by region and year, in format 0 = no outbreaks,
#'   1 = 1 or more outbreak(s)
#' @param chain_cov = Chain covariance
#' @param adapt = 0/1 flag indicating which type of calculation to use for proposition value
#' @param like_current = Current accepted likelihood value
#' @param n_reps Number of times to repeat calculations to get average likelihood at each step
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start = 0, only vaccinated individuals
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity
#' @param prior_type Text indicating which type of calculation to use for prior probability
#'  If prior_type = "zero", prior probability is always zero
#'  If prior_type = "flat", prior probability is zero if FOI/R0 in designated ranges, -Inf otherwise
#'  If prior_type = "exp", prior probability is given by dexp calculation on FOI/R0 values
#'  If prior_type = "norm", prior probability is given by dnorm calculation on parameter values
#' @param dt time increment in days (must be 1 or 5)
#' @param enviro_data Values of environmental variables (if in use)
#' @param R0_fixed_values Values of R0 to use if not being fitted
#' @param vaccine_efficacy Vaccine efficacy (set to NULL if being varied as a parameter)
#' @param p_obs_severe Probability of observation of severe infection (set to NULL if being varied as a parameter)
#' @param p_obs_death Probability of observation of death (set to NULL if being varied as a parameter)
MCMC_step <- function(type=NULL,param=c(),pars_min=NULL,pars_max=NULL,input_data=list(),
                               obs_sero_data=NULL,obs_case_data=NULL,obs_outbreak_data=NULL,chain_cov=1,adapt=0,
                               like_current=-Inf,n_reps=1,mode_start=0,prior_type="zero",dt=1.0,enviro_data=NULL,
                               R0_fixed_values=c(),vaccine_efficacy=NULL,p_obs_severe=NULL,p_obs_death=NULL) {

  #Propose new parameter values
  param_prop=param_prop_setup(param,chain_cov,adapt)

  #Calculate likelihood using single_like_calc function
  like_prop=single_like_calc(type,param_prop,pars_min,pars_max,input_data,obs_sero_data,obs_case_data,
                             obs_outbreak_data,n_reps,mode_start,prior_type,dt,enviro_data,R0_fixed_values,
                             vaccine_efficacy,p_obs_severe,p_obs_death)

  if(is.finite(like_prop)==FALSE) {
    p_accept = -Inf
  } else {
    p_accept= like_prop - like_current
    if(is.na(p_accept) ){ p_accept = -Inf}
  }

  ## accept/reject step:
  tmp = runif(1)
  if(tmp<min(exp(p_accept),1)) { # accept:
    param = param_prop
    like_current = like_prop
    accept = 1
  } else { # reject:
    accept = 0
  }

  output=list(param=param, param_prop=param_prop, like_current=like_current,like_prop=like_prop,accept = accept)

  return(output)
}
#-------------------------------------------------------------------------------
#' @title single_like_calc
#'
#' @description Function which calculates and outputs likelihood
#'
#' @details TBA
#'
#' @param type Type of parameter set (FOI only, FOI+R0, FOI and/or R0 coefficients associated with environmental
#'   covariates); choose from "FOI","FOI+R0","FOI enviro","FOI+R0 enviro"
#' @param param_prop Log values of proposed parameters
#' @param pars_min Lower limits of parameter values if specified
#' @param pars_max Upper limits of parameter values if specified
#' @param input_data List of population and vaccination data for multiple regions (created using data input creation
#'   code and usually loaded from RDS file)
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param obs_outbreak_data Outbreak Y/N data for comparison, by region and year, in format 0 = no outbreaks,
#'   1 = 1 or more outbreak(s)
#' @param n_reps Number of times to repeat calculations to get average likelihood at each step
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start = 0, only vaccinated individuals
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity
#' @param prior_type Text indicating which type of calculation to use for prior probability
#'  If prior_type = "zero", prior probability is always zero
#'  If prior_type = "flat", prior probability is zero if FOI/R0 in designated ranges, -Inf otherwise
#'  If prior_type = "exp", prior probability is given by dexp calculation on FOI/R0 values
#'  If prior_type = "norm", prior probability is given by dnorm calculation on parameter values
#' @param dt time increment in days (must be 1 or 5)
#' @param enviro_data Values of environmental variables (if in use)
#' @param R0_fixed_values Values of R0 to use if not being fitted
#' @param vaccine_efficacy Vaccine efficacy (set to NULL if being varied as a parameter)
#' @param p_obs_severe Probability of observation of severe infection (set to NULL if being varied as a parameter)
#' @param p_obs_death Probability of observation of death (set to NULL if being varied as a parameter)
single_like_calc <- function(type=NULL,param_prop=c(),pars_min=NULL,pars_max=NULL,input_data=list(),
                             obs_sero_data=NULL,obs_case_data=NULL,obs_outbreak_data=NULL,
                             n_reps=1,mode_start=0,prior_type="zero",dt=1.0,enviro_data=NULL,R0_fixed_values=c(),
                             vaccine_efficacy=NULL,p_obs_severe=NULL,p_obs_death=NULL) {

  regions=names(table(input_data$region_labels))
  n_regions=length(regions)
  p_severe=0.12
  p_death_severe=0.47

  #Get vaccine efficacy and calculate associated prior
  if(is.numeric(vaccine_efficacy)==FALSE){
    vaccine_efficacy=exp(param_prop[names(param_prop)=="vacc_eff"])
    prior_vacc=log(dtrunc(vaccine_efficacy,"norm",a=0,b=1,mean=0.975,sd=0.05))
  } else {prior_vacc=0}
  prior_report=0
  if(is.numeric(p_obs_severe)==FALSE){
    p_obs_severe=as.numeric(exp(param_prop[names(param_prop)=="p_obs_severe"]))
    if(p_obs_severe<exp(pars_min[names(pars_min)=="p_obs_severe"])){prior_report=-Inf}
    if(p_obs_severe>exp(pars_max[names(pars_max)=="p_obs_severe"])){prior_report=-Inf}
    }
  if(is.numeric(p_obs_death)==FALSE){
    p_obs_death=as.numeric(exp(param_prop[names(param_prop)=="p_obs_death"]))
    if(p_obs_death<exp(pars_min[names(pars_min)=="p_obs_death"])){prior_report=-Inf}
    if(p_obs_death>exp(pars_max[names(pars_max)=="p_obs_death"])){prior_report=-Inf}
    }
  #if(min(p_obs_severe,p_obs_death)<0){prior_report=-Inf}
  #if(max(p_obs_severe,p_obs_death)>1){prior_report=-Inf}

  #Get FOI and R0 values and calculate associated prior
  FOI_R0_data=mcmc_FOI_R0_setup(type,prior_type,regions,param_prop,enviro_data,R0_fixed_values,pars_min,pars_max)
  FOI_values=FOI_R0_data$FOI_values
  R0_values=FOI_R0_data$R0_values
  prior_prop=prior_vacc+prior_report+FOI_R0_data$prior+sum(dnorm(log(c(p_obs_severe,p_obs_death)),
                                                                 mean = 0,sd = 30,log = TRUE))
  #if(min(FOI_values)<1.0e-8){prior_prop=-Inf}
  #if(type %in% c("FOI+R0","FOI+R0 enviro")){if(min(R0_values)<0.1){prior_prop=-Inf}}

  ### if prior finite, evaluate likelihood ###
  if (is.finite(prior_prop)) {
    sero_like_values=cases_like_values=deaths_like_values=outbreak_like_values=rep(NA,n_regions)
    years_sero=years_outbreak=1940
    flag_skip=0

    for(n_region in 1:n_regions){
      if(flag_skip==0){
        #Set up model inputs
        flag_sero=flag_cases=flag_outbreak_risk=0
        years_outbreak=NA
        region=regions[n_region]

        if(is.null(obs_sero_data)==FALSE){
          obs_sero_data_selected=subset(obs_sero_data,obs_sero_data$gadm36==region)
          if(length(obs_sero_data_selected$year)>0){
            flag_sero=1
            years_sero=obs_sero_data_selected$year
          }
        }
        if(is.null(obs_case_data)==FALSE){
          obs_case_data_selected=obs_case_data[obs_case_data$region==region,]
          years_outbreak=obs_case_data_selected$year
          n_years_outbreak=length(years_outbreak)
          if(n_years_outbreak>0){flag_cases=1}
        }
        if(is.null(obs_outbreak_data)==FALSE){
          obs_outbreak_data_selected=obs_outbreak_data[obs_outbreak_data$region==region,]
          years_outbreak_risk=obs_outbreak_data_selected$year
          n_years_outbreak_risk=length(years_outbreak_risk)
          if(n_years_outbreak_risk>0){
            if(flag_cases==1){ assert_that(n_years_outbreak_risk==n_years_outbreak)} else {
              years_outbreak=years_outbreak_risk
              n_years_outbreak=n_years_outbreak_risk
            }
            flag_outbreak_risk=1
          }
        }

        if(flag_sero+flag_cases+flag_outbreak_risk>0){
          vacc_data=input_data$vacc_data[n_region,,]
          pop_data=input_data$pop_data[n_region,,]
          year_end=max(years_sero,years_outbreak,na.rm=TRUE)+1
          year_data_begin=min(years_sero,years_outbreak,na.rm=TRUE)

          #Run model
          if(flag_sero==1){
            model_output=Full_Model_Run_OD(FOI_values[n_region],R0_values[n_region],vacc_data,pop_data,
                                           year0=min(input_data$years_labels),mode_start,n_particles=n_reps,
                                           n_threads=n_reps,year_end,year_data_begin,vaccine_efficacy,dt=dt)
          } else {
            model_output=case_data_generate(FOI_values[n_region],R0_values[n_region],vacc_data,pop_data,
                                            year0=min(input_data$years_labels),mode_start,n_reps,
                                            year_end,year_data_begin,vaccine_efficacy,dt)
          }

          #Compile outbreak/case data and calculate likelihood, if outbreak/case data available
          if(max(flag_cases,flag_outbreak_risk)==1){
            blank=array(data=rep(0,n_reps*n_years_outbreak),dim=c(n_reps,n_years_outbreak))
            annual_data=list(obs_cases=blank,obs_deaths=blank)
            for(i in 1:n_reps){
              for(n_year in 1:n_years_outbreak){
                year=years_outbreak[n_year]
                if(flag_sero==1){
                  cases=sum(model_output$C[,i,model_output$year[1,]==year])
                } else {
                  cases=model_output$C[i,model_output$year==year]
                }
                severe_cases=rbinom(1,floor(cases),p_severe)
                deaths=rbinom(1,floor(severe_cases),p_death_severe)
                annual_data$obs_deaths[i,n_year]=rbinom(1,floor(deaths),p_obs_death)
                annual_data$obs_cases[i,n_year]=annual_data$obs_deaths[i,n_year]+rbinom(1,floor(severe_cases-deaths),
                                                                                        p_obs_severe)
              }
            }

            if(flag_cases==1){
              cases_like_values[n_region]=cases_compare(annual_data,obs_case_data_selected)/n_reps
              if(is.infinite(cases_like_values[n_region])==TRUE){flag_skip=1}
              deaths_like_values[n_region]=deaths_compare(annual_data,obs_case_data_selected)/n_reps
              if(is.infinite(deaths_like_values[n_region])==TRUE){flag_skip=1}
            }

            if(flag_outbreak_risk==1){
              frac=1.0/n_reps
              outbreak_risk=rep(0,n_years_outbreak)
              for(i in 1:n_reps){
                for(n_year in 1:n_years_outbreak){
                  if(annual_data$obs_cases[i,n_year]>0){outbreak_risk[n_year]=outbreak_risk[n_year]+frac}
                }
              }
              for(n_year in 1:n_years_outbreak){
                if(outbreak_risk[n_year]<1.0e-4){outbreak_risk[n_year]=1.0e-4}
                if(outbreak_risk[n_year]>0.9999){outbreak_risk[n_year]=0.9999}
              }
              outbreak_like_values[n_region]=outbreak_risk_compare(model_outbreak_risk=outbreak_risk,
                                                                   obs_data=obs_outbreak_data_selected$outbreak_yn)
              if(is.infinite(outbreak_like_values[n_region])==TRUE){flag_skip=1}
            }
          }

          #Get seroprevalence likelihood if serological data available
          if(flag_sero==1){
            sero_like_values_all=sero_compare_multiparticle(model_output,obs_sero_data_selected)
            sero_like_values[n_region]=mean(sero_like_values_all,na.rm=TRUE)
            if(is.infinite(sero_like_values[n_region])==TRUE){flag_skip=1}
          }
          model_output<-NULL
        }
      }
    }

    if(flag_skip==0){
      likelihood=sum(c(prior_prop,mean(sero_like_values,na.rm=TRUE),mean(cases_like_values,na.rm=TRUE),
                       mean(deaths_like_values,na.rm=TRUE),mean(outbreak_like_values,na.rm=TRUE)),na.rm=TRUE)
    } else {
      likelihood=-Inf
    }

  } else {
    likelihood=-Inf
  }

  return(likelihood)
}
#-------------------------------------------------------------------------------
#' @title mcmc_checks
#'
#' @description Perform checks on MCMC inputs
#'
#' @details TBA
#'
#' @param type Type of parameter set (FOI only, FOI+R0, FOI and/or R0 coefficients associated with environmental
#'   covariates); choose from "FOI","FOI+R0","FOI enviro","FOI+R0 enviro"
#' @param pars_ini = Initial parameter values
#' @param n_params = Number of parameters (equal to length of pars_ini)
#' @param prior_type = Text indicating which type of calculation to use for prior probability
#'  If prior_type = "zero", prior probability is always zero
#'  If prior_type = "flat", prior probability is zero if FOI/R0 in designated ranges, -Inf otherwise
#'  If prior_type = "exp", prior probability is given by dexp calculation on FOI/R0 values
#'  If prior_type = "norm", prior probability is given by dnorm calculation on parameter values
#' @param n_regions = Number of regions
#' @param enviro_data = Values of environmental variables (if in use)
#' @param R0_fixed_values = Values of R0 to use if not being fitted
#' @param vaccine_efficacy = Input value of vaccine efficacy if fixed
#' @param p_obs_severe = Input value of severe case reporting probability if fixed
#' @param p_obs_death = Input value of fatal case reporting probability if fixed
mcmc_checks <- function(type=NULL,pars_ini=c(),n_params=1,prior_type=NULL,n_regions=1,enviro_data=NULL,
                        R0_fixed_values=NULL,vaccine_efficacy=NULL,p_obs_severe=NULL,p_obs_death=NULL){

  param_names=names(pars_ini)
  assert_that(is.null(param_names)==FALSE)
  assert_that(type %in% c("FOI+R0","FOI","FOI+R0 enviro","FOI enviro"))
  assert_that(prior_type %in% c("zero","flat","exp","norm"))
  n_params=length(pars_ini)
  if(is.numeric(vaccine_efficacy)==TRUE){flag_vacc_eff=0} else {
    flag_vacc_eff=1
    assert_that("vacc_eff" %in% param_names)
  }
  if(is.numeric(p_obs_severe)==TRUE){flag_severe=0} else {
    flag_severe=1
    assert_that("p_obs_severe" %in% param_names)
  }
  if(is.numeric(p_obs_death)==TRUE){flag_death=0} else {
    flag_death=1
    assert_that("p_obs_death" %in% param_names)
  }
  if(is.null(enviro_data)==FALSE){
    vars=names(enviro_data[c(2:ncol(enviro_data))])
    n_vars=length(vars)
  }

  if(type=="FOI+R0"){
    assert_that(n_params==(2*n_regions)+flag_vacc_eff+flag_severe+flag_death)
  }
  if(type=="FOI"){
    assert_that(is.null(R0_fixed_values)==FALSE)
    assert_that(length(R0_fixed_values)==n_regions)
    assert_that(n_params==n_regions+flag_vacc_eff+flag_severe+flag_death)
  }
  if(type=="FOI+R0 enviro"){
    assert_that(is.null(enviro_data)==FALSE)
    assert_that(n_params==(2*n_vars)+flag_vacc_eff+flag_severe+flag_death)
  }
  if(type=="FOI enviro"){
    assert_that(is.null(enviro_data)==FALSE)
    assert_that(is.null(R0_fixed_values)==FALSE)
    assert_that(length(R0_fixed_values)==n_regions)
    assert_that(n_params==n_vars+flag_vacc_eff+flag_severe+flag_death)
  }

  return(NULL)
}
#-------------------------------------------------------------------------------
#' @title mcmc_prelim_fit
#'
#' @description Test multiple sets of parameters randomly drawn from range between maximum and minimum
#' values in order to find approximate values giving maximum likelihood
#' @param n_iterations = Number of times to run and adjust maximum/minimum
#' @param n_param_sets = Number of parameter sets to run in each iteration
#' @param n_bounds = Number of parameter sets (with highest likelihood values) to take at each iteration to create new
#' maximum/minimum values
#' @param type Type of parameter set (FOI only, FOI+R0, FOI and/or R0 coefficients associated with environmental
#'   covariates); choose from "FOI","FOI+R0","FOI enviro","FOI+R0 enviro"
#' @param pars_min Initial lower limits of parameter values
#' @param pars_max Initial upper limits of parameter values
#' @param input_data List of population and vaccination data for multiple regions (created using data input creation
#'   code and usually loaded from RDS file)
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param obs_outbreak_data Outbreak Y/N data for comparison, by region and year, in format 0 = no outbreaks,
#'   1 = 1 or more outbreak(s)
#' @param n_reps Number of times to repeat calculations to get average likelihood at each step
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start = 0, only vaccinated individuals
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity
#' @param prior_type Text indicating which type of calculation to use for prior probability
#'  If prior_type = "zero", prior probability is always zero
#'  If prior_type = "flat", prior probability is zero if FOI/R0 in designated ranges, -Inf otherwise
#'  If prior_type = "exp", prior probability is given by dexp calculation on FOI/R0 values
#'  If prior_type = "norm", prior probability is given by dnorm calculation on parameter values
#' @param dt time increment in days (must be 1 or 5)
#' @param enviro_data Values of environmental variables (if in use)
#' @param R0_fixed_values Values of R0 to use if not being fitted
#' @param vaccine_efficacy Vaccine efficacy (set to NULL if being varied as a parameter)
#' @param p_obs_severe Probability of observation of severe infection (set to NULL if being varied as a parameter)
#' @param p_obs_death Probability of observation of death (set to NULL if being varied as a parameter)
mcmc_prelim_fit <- function(n_iterations=1,n_param_sets=1,n_bounds=1,
                            type=NULL,pars_min=NULL,pars_max=NULL,input_data=list(),
                            obs_sero_data=list(),obs_case_data=list(),obs_outbreak_data=list(),
                            n_reps=1,mode_start=0,prior_type="zero",dt=1.0,enviro_data=NULL,R0_fixed_values=c(),
                            vaccine_efficacy=NULL,p_obs_severe=NULL,p_obs_death=NULL){
  #TODO - Add assertthat functions
  assert_that(length(pars_min)==length(pars_max))
  assert_that(type %in% c("FOI+R0","FOI","FOI+R0 enviro","FOI enviro"))
  assert_that(prior_type %in% c("zero","flat","exp","norm"))

  best_fit_results=list()
  n_params=length(pars_min)
  extra_params=c()
  if(is.null(vaccine_efficacy)==TRUE){extra_params=append(extra_params,"vacc_eff")}
  if(is.null(p_obs_severe)==TRUE){extra_params=append(extra_params,"p_obs_severe")}
  if(is.null(p_obs_death)==TRUE){extra_params=append(extra_params,"p_obs_death")}
  param_names=create_param_labels(type,input_data,enviro_data,extra_params)
  assert_that(length(param_names)==n_params)
  names(pars_min)=names(pars_max)=param_names
  xlabels=param_names
  for(i in 1:n_params){xlabels[i]=substr(xlabels[i],1,15)}
  ylabels=10^c(-8,-6,-4,-3,-2,-1,0,1)
  par(mar=c(6,2,1,1))
  ylim=c(min(pars_min),max(pars_max))

  for(iteration in 1:n_iterations){
    cat("\nIteration: ",iteration,"\n",sep="")
    all_param_sets <- lhs(n=n_param_sets,rect=cbind(pars_min,pars_max))
    results=data.frame()

    for(set in 1:n_param_sets){
      cat(set,"\t",sep="")
      param_prop=all_param_sets[set,]
      names(param_prop)=param_names
      like_prop=single_like_calc(type,param_prop,pars_min,pars_max,input_data,obs_sero_data,obs_case_data,
                                 obs_outbreak_data,n_reps,mode_start,prior_type,dt,enviro_data,R0_fixed_values,
                                 vaccine_efficacy,p_obs_severe,p_obs_death)
      results<-rbind(results,c(set,exp(param_prop),like_prop))
      if(set==1){colnames(results)=c("set",param_names,"LogLikelihood")}

    }
    results<-results[order(results$LogLikelihood,decreasing=TRUE), ]
    best_fit_results[[iteration]]=results

    pars_min_new=pars_max_new=rep(0,n_params)
    for(i in 1:n_params){
      pars_min_new[i]=min(log(results[c(1:n_bounds),i+1]))
      pars_max_new[i]=max(log(results[c(1:n_bounds),i+1]))
    }
    names(pars_min_new)=names(pars_max_new)=param_names

    matplot(x=c(1:n_params),y=log(t(results[c(1:n_bounds),c(1:n_params)+1])),type="p",pch=16,col=1,
            xaxt="n",yaxt="n",xlab="",ylab="",ylim=ylim)
    axis(side=1,at=c(1:n_params),labels=xlabels,las=2,cex.axis=0.7)
    axis(side=2,at=log(ylabels),labels=ylabels)
    matplot(x=c(1:n_params),y=pars_min,type="l",col=1,lty=2,add=TRUE)
    matplot(x=c(1:n_params),y=pars_max,type="l",col=1,lty=2,add=TRUE)
    matplot(x=c(1:n_params),y=pars_min_new,type="l",col=2,add=TRUE)
    matplot(x=c(1:n_params),y=pars_max_new,type="l",col=2,add=TRUE)

    pars_min=pars_min_new
    pars_max=pars_max_new
  }

  return(best_fit_results)
}
