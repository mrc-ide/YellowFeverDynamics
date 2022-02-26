# R file for functions used for Markov Chain Monte Carlo fitting (and preliminary maximum-likelihood fitting) in
# YellowFeverDynamics package
#-------------------------------------------------------------------------------
#' @title MCMC2
#'
#' @description Combined MCMC Multi-Region - series of MCMC steps for one or more regions
#'
#' @details This is the master function for running a Markov chain to optimize the parameters of the yellow fever model
#' based on the calculated likelihood of observing supplied data given a particular set of parameters. [TBC]
#'
#' @param pars_ini Initial log values of parameters
#' @param input_data List of population and vaccination data for multiple regions (created using data input creation
#'   code and usually loaded from RDS file)
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no.
#'   cases/no. deaths
#' @param obs_outbreak_data Outbreak Y/N data for comparison, by region and year, in format 0 = no outbreaks,
#'   1 = 1 or more outbreak(s)
#' @param Niter Total number of steps to run
#' @param type Type of parameter set (FOI only, FOI+R0, FOI and/or R0 coefficients associated with environmental
#'   covariates); choose from "FOI","FOI+R0","FOI enviro","FOI+R0 enviro"
#' @param pars_min Lower limits of parameter values if specified
#' @param pars_max Upper limits of parameter values if specified
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
#' @param filename_prefix Prefix of names for output files
#' '
#' @export
#'
MCMC2 <- function(pars_ini=c(),input_data=list(),obs_sero_data=NULL,obs_case_data=NULL,obs_outbreak_data=NULL,Niter=1,
                  type=NULL,pars_min=NULL,pars_max=NULL,n_reps=1,mode_start=0,prior_type="zero",dt=1.0,
                  enviro_data=NULL,R0_fixed_values=NULL,vaccine_efficacy=NULL,p_obs_severe=NULL,p_obs_death=NULL,
                  filename_prefix="Chain"){

  regions=names(table(input_data$region_labels))
  n_regions=length(regions)
  n_params=length(pars_ini)

  checks<-mcmc_checks(type,pars_ini,n_params,prior_type,n_regions,enviro_data,R0_fixed_values,
                      vaccine_efficacy,p_obs_severe,p_obs_death)
  const_list=list(type=type,pars_min=pars_min,pars_max=pars_max,n_reps=n_reps,mode_start=mode_start,
                  prior_type=prior_type,dt=dt,enviro_data=enviro_data,R0_fixed_values=R0_fixed_values,
                  vaccine_efficacy=vaccine_efficacy,p_obs_severe=p_obs_severe,p_obs_death=p_obs_death)

  ### find posterior probability at start ###
  out = MCMC_step2(param=pars_ini,input_data,obs_sero_data,obs_case_data,obs_outbreak_data,
                   chain_cov=1,adapt=0,like_current=-Inf,const_list)

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
    out = MCMC_step2(param,input_data,obs_sero_data,obs_case_data,obs_outbreak_data,chain_cov,adapt,like_current,
                     const_list)
  }

  param_out=exp(out$param)
  names(param_out)=names(pars_ini)

  return(param_out)
}
#-------------------------------------------------------------------------------
#' @title MCMC_step2
#'
#' @description Single MCMC step - one or more regions
#'
#' @details This function runs a single step in a Markov chain set up using the function mcmc(). It proposes a
#' set of parameters using the param_prop_setup() function, calculates the likelihood of observing the observed data
#' based on that proposed parameter set, accepts or rejects the proposed parameter set based on the calculated
#' likelihood and existing chain information, then returns the next line of information for the chain to mcmc(). [TBC]
#'
#' @param param Log values of parameters
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
#' @param const_list = List of constant parameters/flags/etc. loaded to mcmc2() (type,pars_min,pars_max,n_reps,
#' mode_start,prior_type,dt=dt,enviro_data,R0_fixed_values,vaccine_efficacy,p_obs_severe,p_obs_death)
#'
#' @export
#'
MCMC_step2 <- function(param=c(),input_data=list(),obs_sero_data=NULL,obs_case_data=NULL,obs_outbreak_data=NULL,
                       chain_cov=1,adapt=0,like_current=-Inf,const_list=list()) {

  #Propose new parameter values
  param_prop=param_prop_setup(param,chain_cov,adapt)

  #Calculate likelihood using single_like_calc function
  like_prop=single_like_calc2(param_prop,input_data,obs_sero_data,obs_case_data,obs_outbreak_data,const_list)

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
#' @title single_like_calc2
#'
#' @description Function which calculates and outputs likelihood
#'
#' @details This function calculates the total likelihood of observing a set of observations (across multiple regions
#' and data types) for a given proposed parameter set. [TBC]
#'
#' @param param_prop Log values of proposed parameters
#' @param input_data List of population and vaccination data for multiple regions (created using data input creation
#'   code and usually loaded from RDS file)
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param obs_outbreak_data Outbreak Y/N data for comparison, by region and year, in format 0 = no outbreaks,
#'   1 = 1 or more outbreak(s)
#' @param const_list = List of constant parameters/flags/etc. loaded to mcmc2() (type,pars_min,pars_max,n_reps,
#' mode_start,prior_type,dt=dt,enviro_data,R0_fixed_values,vaccine_efficacy,p_obs_severe,p_obs_death)
#'
#' @export
#'
single_like_calc2 <- function(param_prop=c(),input_data=list(),obs_sero_data=NULL,obs_case_data=NULL,
                              obs_outbreak_data=NULL,const_list=list()) {

  regions=names(table(input_data$region_labels))
  n_regions=length(regions)
  p_severe=0.12
  p_death_severe=0.47
  frac=1.0/const_list$n_reps

  #Get vaccine efficacy and calculate associated prior
  if(is.numeric(const_list$vaccine_efficacy)==FALSE){
    vaccine_efficacy=exp(param_prop[names(param_prop)=="vacc_eff"])
    prior_vacc=log(dtrunc(vaccine_efficacy,"norm",a=0,b=1,mean=0.975,sd=0.05))
  } else {prior_vacc=0}
  prior_report=0
  if(is.numeric(const_list$p_obs_severe)==FALSE){
    p_obs_severe=as.numeric(exp(param_prop[names(param_prop)=="p_obs_severe"]))
    if(p_obs_severe<exp(pars_min[names(pars_min)=="p_obs_severe"])){prior_report=-Inf}
    if(p_obs_severe>exp(pars_max[names(pars_max)=="p_obs_severe"])){prior_report=-Inf}
    }
  if(is.numeric(const_list$p_obs_death)==FALSE){
    p_obs_death=as.numeric(exp(param_prop[names(param_prop)=="p_obs_death"]))
    if(p_obs_death<exp(pars_min[names(pars_min)=="p_obs_death"])){prior_report=-Inf}
    if(p_obs_death>exp(pars_max[names(pars_max)=="p_obs_death"])){prior_report=-Inf}
    }

  #Get FOI and R0 values and calculate associated prior
  FOI_R0_data=mcmc_FOI_R0_setup(const_list$type,const_list$prior_type,regions,param_prop,const_list$enviro_data,
                                const_list$R0_fixed_values,const_list$pars_min,const_list$pars_max)
  FOI_values=FOI_R0_data$FOI_values
  R0_values=FOI_R0_data$R0_values
  prior_prop=prior_vacc+prior_report+FOI_R0_data$prior+sum(dnorm(log(c(p_obs_severe,p_obs_death)),
                                                                 mean = 0,sd = 30,log = TRUE))

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
            model_output=Full_Model_Run(FOI_values[n_region],R0_values[n_region],vacc_data,pop_data,
                                        year0=min(input_data$years_labels),const_list$mode_start,
                                        n_particles=const_list$n_reps,n_threads=const_list$n_reps,year_end,
                                        year_data_begin,vaccine_efficacy,dt=const_list$dt)
          } else {
            model_output=case_data_generate(FOI_values[n_region],R0_values[n_region],vacc_data,pop_data,
                                            year0=min(input_data$years_labels),const_list$mode_start,const_list$n_reps,
                                            year_end,year_data_begin,vaccine_efficacy,const_list$dt)
          }

          #Compile outbreak/case data and calculate likelihood, if outbreak/case data available
          if(max(flag_cases,flag_outbreak_risk)==1){
            blank=array(data=rep(0,const_list$n_reps*n_years_outbreak),dim=c(const_list$n_reps,n_years_outbreak))
            annual_data=list(obs_cases=blank,obs_deaths=blank)
            for(i in 1:const_list$n_reps){
              for(n_year in 1:n_years_outbreak){
                year=years_outbreak[n_year]
                if(flag_sero==1){
                  cases=sum(model_output$C[,i,model_output$year[1,]==year])
                } else {
                  cases=model_output$C[i,model_output$year==year]
                }
                severe_cases=rbinom(1,floor(cases),p_severe)
                deaths=rbinom(1,severe_cases,p_death_severe)
                annual_data$obs_deaths[i,n_year]=rbinom(1,deaths,p_obs_death)
                annual_data$obs_cases[i,n_year]=annual_data$obs_deaths[i,n_year]+rbinom(1,severe_cases-deaths,
                                                                                        p_obs_severe)
              }
            }

            if(flag_cases==1){
              cases_like_values[n_region]=cases_compare(annual_data,obs_case_data_selected)/const_list$n_reps
              if(is.infinite(cases_like_values[n_region])==TRUE){flag_skip=1}
              deaths_like_values[n_region]=deaths_compare(annual_data,obs_case_data_selected)/const_list$n_reps
              if(is.infinite(deaths_like_values[n_region])==TRUE){flag_skip=1}
            }

            if(flag_outbreak_risk==1){
              outbreak_risk=rep(0,n_years_outbreak)
              for(i in 1:const_list$n_reps){
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
