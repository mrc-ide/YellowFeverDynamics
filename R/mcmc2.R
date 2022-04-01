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
#' @param pars_ini Initial values of parameters to be fitted. These should always be the log() values of the actual
#'   epidemiological parameters, ordered as follows:
#'   1) Parameters controlling the value of spillover force of infection FOI, either a) a number of FOI values equal to
#'   the total number of regions to be considered or b) a number of environmental coefficients used to calculate FOI
#'   values from environmental covariates equal to the number of environmental covariates listed in the enviro_data
#'   frame. Values should be in alphabetical order by region in case (a) or in the order of the columns in the
#'   environmental data frame in case (b).
#'   2) If the basic reproduction number for human-human transmission R0 is to be fitted (i.e. type is set to
#'   "FOI+R0" or "FOI+R0 enviro"), parameters controlling the value of R0, either a) a number of R0 values equal to
#'   the total number of regions to be considered or b) a number of environmental coefficients used to calculate R0
#'   values from environmental covariates equal to the number of environmental covariates listed in the enviro_data
#'   frame. Values should be in alphabetical order by region in case (a) or in the order of the columns in the
#'   environmental data frame in case (b).
#'   3) Values of the additional parameters (vaccine efficacy vaccine_efficacy, severe case reporting probability
#'   p_rep_severe and fatal case reporting probability p_rep_death) if these are to be fitted, in the order
#'   vaccine_efficacy->p_rep_severe->p_rep_death. If these parameters are to be fitted, the values separately supplied
#'   to this function (see below) should be set to NULL, the default.
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
#' @param enviro_data Data frame containing values of environmental covariates; set to NULL if not in use
#' @param R0_fixed_values Values of R0 to use if only FOI is subject to fitting (i.e. type set to "FOI" or "FOI
#'   enviro"); set to NULL if not in use
#' @param vaccine_efficacy Vaccine efficacy (set to NULL if being varied as a parameter)
#' @param p_rep_severe Probability of observation of severe infection (set to NULL if being varied as a parameter)
#' @param p_rep_death Probability of observation of death (set to NULL if being varied as a parameter)
#' @param filename_prefix Prefix of names for output files
#' '
#' @export
#'
MCMC2 <- function(pars_ini=c(),input_data=list(),obs_sero_data=NULL,obs_case_data=NULL,obs_outbreak_data=NULL,Niter=1,
                  type=NULL,pars_min=NULL,pars_max=NULL,n_reps=1,mode_start=0,prior_type="zero",dt=1.0,
                  enviro_data=NULL,R0_fixed_values=NULL,vaccine_efficacy=NULL,p_rep_severe=NULL,p_rep_death=NULL,
                  filename_prefix="Chain"){

  #Check that initial, minimum and maximum parameters are in vectors of same sizes
  n_params=length(pars_ini)
  assert_that(length(pars_min)==n_params)
  assert_that(length(pars_max)==n_params)

  #TODO - Assert that environmental variables must be in alphabetical order?

  extra_params=c()
  if(is.null(vaccine_efficacy)==TRUE){extra_params=append(extra_params,"vaccine_efficacy")}
  if(is.null(p_rep_severe)==TRUE){extra_params=append(extra_params,"p_rep_severe")}
  if(is.null(p_rep_death)==TRUE){extra_params=append(extra_params,"p_rep_death")}

  #Process input data to check that all regions with sero, case and/or outbreak data supplied are present, remove
  #regions without any supplied data, and add cross-referencing tables for use when calculating likelihood
  input_data=input_data_process2(input_data,obs_sero_data,obs_case_data,obs_outbreak_data)
  regions=names(table(input_data$region_labels)) #Regions in new processed input data list
  n_regions=length(regions)

  #Label parameters according to order and fitting type
  param_names=create_param_labels(type,input_data,enviro_data,extra_params)
  names(pars_ini)=names(pars_min)=names(pars_max)=param_names

  #Run checks to ensure that number of parameters is correct for fitting type and number of regions/environmental covariates
  checks<-mcmc_checks2(pars_ini,n_regions,type,pars_min,pars_max,prior_type,enviro_data,
                       R0_fixed_values,vaccine_efficacy,p_rep_severe,p_rep_death)

  #Set up list of invariant parameter values to supply to other functions
  const_list=list(type=type,pars_min=pars_min,pars_max=pars_max,n_reps=n_reps,mode_start=mode_start,
                  prior_type=prior_type,dt=dt,enviro_data=enviro_data,R0_fixed_values=R0_fixed_values,
                  vaccine_efficacy=vaccine_efficacy,p_rep_severe=p_rep_severe,p_rep_death=p_rep_death)

  ### find posterior probability at start ###
  out = MCMC_step2(param=pars_ini,input_data,obs_sero_data,obs_case_data,obs_outbreak_data,
                   chain_cov=1,adapt=0,like_current=-Inf,const_list)

  #MCMC setup
  chain=chain_prop=posterior_current=posterior_prop=flag_accept=chain_cov_all=NULL
  burnin = min(2*n_params, Niter)
  fileIndex = 0
  chain_cov=1

  #Iterative fitting
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

    #Set output headings
    if(iter==1){
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

    #Output chain to file every 10 iterations; start new file every 10,000 iterations
    if (iter %% 10 == 0){
      if (iter %% 10000 == 0){fileIndex  = iter/10000}

      filename=paste(filename_prefix,fileIndex,".csv",sep="")
      lines=min((fileIndex * 10000+1),iter):iter
      cat("\nIteration ",iter,sep="")
      data_out<-cbind(posterior_current,posterior_prop,exp(chain),flag_accept,exp(chain_prop),chain_cov_all)[lines,]
      write.csv(data_out,filename,row.names=FALSE)
    }

    #Decide whether next iteration will be adaptive or use [TBA]
    if (iter>burnin & runif(1)<0.9){ #adapt
      adapt = 1
      chain_cov  = cov(chain[max(nrow(chain)-10000, 1):nrow(chain),])
    } else {
      adapt = 0
      chain_cov = 1
    }

    #Next iteration in chain
    out = MCMC_step2(param,input_data,obs_sero_data,obs_case_data,obs_outbreak_data,chain_cov,adapt,like_current,
                     const_list)
  }

  #Get final parameter values
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
#' mode_start,prior_type,dt=dt,enviro_data,R0_fixed_values,vaccine_efficacy,p_rep_severe,p_rep_death)
#'
#' @export
#'
MCMC_step2 <- function(param=c(),input_data=list(),obs_sero_data=NULL,obs_case_data=NULL,obs_outbreak_data=NULL,
                       chain_cov=1,adapt=0,like_current=-Inf,const_list=list()) {

  #Propose new parameter values
  param_prop=param_prop_setup(param,chain_cov,adapt)

  #Calculate likelihood using single_like_calc2 function
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
#'   code and usually loaded from RDS file), with cross-reference tables added using input_data_process2 in MCMC2
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param obs_outbreak_data Outbreak Y/N data for comparison, by region and year, in format 0 = no outbreaks,
#'   1 = 1 or more outbreak(s)
#' @param const_list = List of constant parameters/flags/etc. loaded to mcmc2() (type,pars_min,pars_max,n_reps,
#' mode_start,prior_type,dt=dt,enviro_data,R0_fixed_values,vaccine_efficacy,p_rep_severe,p_rep_death)
#'
#' @export
#'
single_like_calc2 <- function(param_prop=c(),input_data=list(),obs_sero_data=NULL,obs_case_data=NULL,
                              obs_outbreak_data=NULL,const_list=list()) {

  regions=input_data$region_labels
  n_regions=length(regions)
  # p_severe_inf=0.12
  # p_death_severe_inf=0.47
  frac=1.0/const_list$n_reps

  #Get vaccine efficacy and calculate associated prior
  if(is.numeric(const_list$vaccine_efficacy)==FALSE){
    vaccine_efficacy=exp(param_prop[names(param_prop)=="vaccine_efficacy"])
    prior_vacc=log(dtrunc(vaccine_efficacy,"norm",a=0,b=1,mean=0.975,sd=0.05))
  } else {
    vaccine_efficacy=const_list$vaccine_efficacy
    prior_vacc=0
  }

  #Get reporting probabilities and check they are within specified bounds
  prior_report=0
  if(is.numeric(const_list$p_rep_severe)==FALSE){
    p_rep_severe=as.numeric(exp(param_prop[names(param_prop)=="p_rep_severe"]))
    if(p_rep_severe<exp(const_list$pars_min[names(const_list$pars_min)=="p_rep_severe"])){prior_report=-Inf}
    if(p_rep_severe>exp(const_list$pars_max[names(const_list$pars_max)=="p_rep_severe"])){prior_report=-Inf}
  } else {
    p_rep_severe=const_list$p_rep_severe
  }
  if(is.numeric(const_list$p_rep_death)==FALSE){
    p_rep_death=as.numeric(exp(param_prop[names(param_prop)=="p_rep_death"]))
    if(p_rep_death<exp(const_list$pars_min[names(const_list$pars_min)=="p_rep_death"])){prior_report=-Inf}
    if(p_rep_death>exp(const_list$pars_max[names(const_list$pars_max)=="p_rep_death"])){prior_report=-Inf}
  } else {
    p_rep_death=const_list$p_rep_death
  }

  #Get FOI and R0 values and calculate associated prior
  FOI_R0_data=mcmc_FOI_R0_setup(const_list$type,const_list$prior_type,regions,param_prop,const_list$enviro_data,
                                const_list$R0_fixed_values,const_list$pars_min,const_list$pars_max)
  FOI_values=FOI_R0_data$FOI_values
  R0_values=FOI_R0_data$R0_values
  prior_prop=prior_vacc+prior_report+FOI_R0_data$prior+sum(dnorm(log(c(p_rep_severe,p_rep_death)),
                                                                 mean = 0,sd = 30,log = TRUE))

  ### If prior finite, evaluate likelihood ###
  if (is.finite(prior_prop)) {

    #Set up data structures to take modelled data corresponding to observed data and likelihood values
    if(is.null(obs_sero_data)==FALSE){
      sero_like_values=model_sero_values=rep(0,nrow(obs_sero_data))
      # model_sero_data=list()
      # blank1=rep(0,nrow(obs_sero_data))
      # blank2=data.frame(samples=blank1,positives=blank1,sero=blank1)
      # for(rep in 1:const_list$n_reps){
      #   model_sero_data[[rep]]=blank2
      # }
    } else {sero_like_values=model_sero_values=NA}
    if(is.null(obs_case_data)==FALSE){
      cases_like_values=deaths_like_values=model_case_values=model_death_values=rep(0,nrow(obs_case_data))
      # model_case_data=list()
      # blank1=rep(0,nrow(obs_case_data))
      # blank2=data.frame(rep_cases=blank1,rep_deaths=blank1)
      # for(rep in 1:const_list$n_reps){
      #   model_case_data[[rep]]=blank2
      # }
    } else {cases_like_values=deaths_like_values=model_case_values=model_death_values=NA}
    if(is.null(obs_outbreak_data)==FALSE){
      outbreak_like_values=model_outbreak_risk_values=rep(0,nrow(obs_outbreak_data))
      # model_outbreak_data=list()
      # blank1=rep(0,nrow(obs_outbreak_data))
      # for(rep in 1:const_list$n_reps){
      #   model_outbreak_data[[rep]]=blank1
      # }
    } else {outbreak_like_values=model_outbreak_risk_values=NA}

    #Model all regions and save relevant output data
    for(n_region in 1:n_regions){

      #Get information on which observed data types are available for considered region
      flag_sero=input_data$flag_sero[n_region]
      flag_case=input_data$flag_case[n_region]
      flag_outbreak=input_data$flag_outbreak[n_region]

      #Get input data on region
      region=regions[n_region]
      vacc_data=input_data$vacc_data[n_region,,]
      pop_data=input_data$pop_data[n_region,,]
      year_end=input_data$year_end[n_region]
      year_data_begin=input_data$year_data_begin[n_region]

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

      #Compile outbreak/case data if needed
      if(max(flag_case,flag_outbreak)==1){
        if(flag_case==1){
          case_line_list=input_data$case_line_list[[n_region]]
          years_outbreak=obs_case_data$year[case_line_list]
          n_years_outbreak=length(case_line_list)
        } else {
          years_outbreak=obs_outbreak_data$year[input_data$outbreak_line_list[[n_region]]]
          n_years_outbreak=length(input_data$outbreak_line_list[[n_region]])
        }
        blank=array(data=rep(0,const_list$n_reps*n_years_outbreak),dim=c(const_list$n_reps,n_years_outbreak))
        annual_data=list(rep_cases=blank,rep_deaths=blank)
        for(i in 1:const_list$n_reps){
          for(n_year in 1:n_years_outbreak){
            year=years_outbreak[n_year]
            if(flag_sero==1){
              infs=sum(model_output$C[,i,model_output$year[1,]==year])
            } else {
              infs=model_output$C[i,model_output$year==year]
            }
            severe_infs=rbinom(1,floor(infs),p_severe_inf)
            deaths=rbinom(1,severe_infs,p_death_severe_inf)
            annual_data$rep_deaths[i,n_year]=rbinom(1,deaths,p_rep_death)
            annual_data$rep_cases[i,n_year]=annual_data$rep_deaths[i,n_year]+rbinom(1,severe_infs-deaths,
                                                                                    p_rep_severe)
          }
        }

        if(flag_case==1){
          for(rep in 1:const_list$n_reps){
            model_case_values[case_line_list]=model_case_values[case_line_list]+annual_data$rep_cases[rep,]
            model_death_values[case_line_list]=model_death_values[case_line_list]+annual_data$rep_deaths[rep,]
          }
        }

        if(flag_outbreak==1){
          outbreak_risk=rep(0,n_years_outbreak)
          for(i in 1:const_list$n_reps){
            for(n_year in 1:n_years_outbreak){
              if(annual_data$rep_cases[i,n_year]>0){outbreak_risk[n_year]=outbreak_risk[n_year]+frac}
            }
          }
          for(n_year in 1:n_years_outbreak){
            if(outbreak_risk[n_year]<1.0e-4){outbreak_risk[n_year]=1.0e-4}
            if(outbreak_risk[n_year]>0.9999){outbreak_risk[n_year]=0.9999}
          }
          for(i in 1:const_list$n_reps){
            model_outbreak_data[input_data$outbreak_line_list[[n_region]]]=outbreak_risk
          }
        }
      }

      #Compile seroprevalence data if necessary
      if(flag_sero==1){
        sero_line_list=input_data$sero_line_list[[n_region]]
        for(i in 1:const_list$n_reps){
          sero_results=sero_calculate2(obs_sero_data[sero_line_list,],model_data=list(day=model_output$day[i,],
                                                                                      year=model_output$year[i,],S=t(model_output$S[,i,]),E=t(model_output$E[,i,]),
                                                                                      I=t(model_output$I[,i,]),R=t(model_output$R[,i,]),V=t(model_output$V[,i,])))
          model_sero_values[sero_line_list]=model_sero_values[sero_line_list]+(sero_results$positives/sero_results$samples)
        }
      }
      model_output<-NULL
    }

    #Likelihood of observing serological data
    if(is.null(obs_sero_data)==FALSE){
      model_sero_values=model_sero_values*frac
      sero_like_values=sero_like_values+lgamma(obs_sero_data$samples+1)-
        lgamma(obs_sero_data$positives+1)-lgamma(obs_sero_data$samples-obs_sero_data$positives+1)+
        obs_sero_data$positives*log(model_sero_values)+
        (obs_sero_data$samples-obs_sero_data$positives)*log(1.0-model_sero_values)
    }
    #Likelihood of observing annual case/death data
    if(is.null(obs_case_data)==FALSE){
      for(i in 1:length(model_case_values)){
        model_case_values[i]=max(model_case_values[i]*frac,0.1)
        model_death_values[i]=max(model_death_values[i]*frac,0.1)
      }
      cases_like_values=dnbinom(x=obs_case_data$cases,mu=model_case_values,
                                size=rep(1,length(obs_case_data$cases)),log=TRUE)
      deaths_like_values=dnbinom(x=obs_case_data$deaths,mu=model_death_values,
                                 size=rep(1,length(obs_case_data$deaths)),log=TRUE)
    }
    #Likelihood of observing annual outbreak Y/N data
    if(is.null(obs_outbreak_data)==FALSE){
      outbreak_like_values=outbreak_risk_compare(model_outbreak_risk=model_outbreak_data,
                                                 obs_data=obs_outbreak_data$outbreak_yn)
    }

    # likelihood=sum(c(prior_prop,mean(sero_like_values,na.rm=TRUE),mean(cases_like_values,na.rm=TRUE),
    #                  mean(deaths_like_values,na.rm=TRUE),mean(outbreak_like_values,na.rm=TRUE)),na.rm=TRUE)
    # likelihood=sum(c(prior_prop,sum(sero_like_values,na.rm=TRUE),sum(cases_like_values,na.rm=TRUE),
    #                  sum(deaths_like_values,na.rm=TRUE),sum(outbreak_like_values,na.rm=TRUE)),na.rm=TRUE)
    likelihood=prior_prop+mean(c(sum(sero_like_values,na.rm=TRUE),sum(cases_like_values,na.rm=TRUE),
                     sum(deaths_like_values,na.rm=TRUE),sum(outbreak_like_values,na.rm=TRUE)),na.rm=TRUE)

  } else {
    likelihood=-Inf
  }

  return(likelihood)
}
#-------------------------------------------------------------------------------
#' @title input_data_process2
#'
#' @description TBA
#'
#' @details TBA
#'
#' @param input_data List of population and vaccination data for multiple regions (created using data input creation
#'   code and usually loaded from RDS file)
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param obs_outbreak_data Outbreak Y/N data for comparison, by region and year, in format 0 = no outbreaks,
#'   1 = 1 or more outbreak(s)
#'
#' @export
#'
input_data_process2 <- function(input_data=list(),obs_sero_data=NULL,obs_case_data=NULL,obs_outbreak_data=NULL){

  regions_input_data=input_data$region_labels
  #TODO - Make sure regions always in alphabetical order?
  #if(table(sort(regions_input_data)==regions_input_data)[[TRUE]]==length(regions_input_data)){}
  regions_sero_com=names(table(obs_sero_data$adm1))
  regions_case_com=names(table(obs_case_data$adm1))
  regions_outbreak_com=names(table(obs_outbreak_data$adm1))
  regions_sero_unc=regions_case_unc=regions_outbreak_unc=c()

  for(region in regions_sero_com){regions_sero_unc=append(regions_sero_unc,strsplit(region,",")[[1]])}
  for(region in regions_case_com){regions_case_unc=append(regions_case_unc,strsplit(region,",")[[1]])}
  for(region in regions_outbreak_com){regions_outbreak_unc=append(regions_outbreak_unc,strsplit(region,",")[[1]])}
  regions_sero_unc=names(table(regions_sero_unc))
  regions_case_unc=names(table(regions_case_unc))
  regions_outbreak_unc=names(table(regions_outbreak_unc))
  regions_all_data_unc=names(table(c(regions_sero_unc,regions_case_unc,regions_outbreak_unc)))

  for(region in regions_all_data_unc){
    if(region %in% regions_input_data==FALSE){
      cat("\nInput data error - ",region," does not appear in input data\n",sep="")
      stop()
    }
  }

  input_regions_check=regions_input_data %in% regions_all_data_unc
  regions_input_data_new=regions_input_data[input_regions_check]
  n_regions_input_data=length(regions_input_data_new)

  #TODO - Skip steps below if the input data already has the same set of regions and the cross-reference tables

  blank=rep(0,n_regions_input_data)
  flag_sero=flag_case=flag_outbreak=year_end=blank
  year_data_begin=rep(Inf,n_regions_input_data)
  sero_line_list=case_line_list=outbreak_line_list=list()
  for(i in 1:n_regions_input_data){
    region=regions_input_data_new[i]
    if(region %in% regions_sero_unc){
      flag_sero[i]=1
      sero_line_list[[i]]=c(0)
      for(j in 1:nrow(obs_sero_data)){
        if(grepl(region,obs_sero_data$adm1[j])==TRUE){
          sero_line_list[[i]]=append(sero_line_list[[i]],j)
          year_data_begin[i]=min(obs_sero_data$year[j],year_data_begin[i])
          year_end[i]=max(obs_sero_data$year[j]+1,year_end[i])
        }
      }
      sero_line_list[[i]]=sero_line_list[[i]][c(2:length(sero_line_list[[i]]))]
    }
    if(region %in% regions_case_unc){
      flag_case[i]=1
      case_line_list[[i]]=c(0)
      for(j in 1:nrow(obs_case_data)){
        if(grepl(region,obs_case_data$adm1[j])==TRUE){
          case_line_list[[i]]=append(case_line_list[[i]],j)
          year_data_begin[i]=min(obs_case_data$year[j],year_data_begin[i])
          year_end[i]=max(obs_case_data$year[j]+1,year_end[i])
        }
      }
      case_line_list[[i]]=case_line_list[[i]][c(2:length(case_line_list[[i]]))]
    }
    if(region %in% regions_outbreak_unc){
      flag_outbreak[i]=1
      outbreak_line_list[[i]]=c(0)
      for(j in 1:nrow(obs_outbreak_data)){
        if(grepl(region,obs_outbreak_data$adm1[j])==TRUE){
          outbreak_line_list[[i]]=append(outbreak_line_list[[i]],j)
          year_data_begin[i]=min(obs_outbreak_data$year[j],year_data_begin[i])
          year_end[i]=max(obs_outbreak_data$year[j]+1,year_end[i])
        }
      }
      outbreak_line_list[[i]]=outbreak_line_list[[i]][c(2:length(outbreak_line_list[[i]]))]
    }
  }

  input_data_new=list(region_labels=input_data$region_labels[input_regions_check],
                      years_labels=input_data$years_labels,age_labels=input_data$age_labels,
                      vacc_data=input_data$vacc_data[input_regions_check,,],
                      pop_data=input_data$pop_data[input_regions_check,,],year_data_begin=year_data_begin,
                      year_end=year_end,flag_sero=flag_sero,flag_case=flag_case,flag_outbreak=flag_outbreak,
                    sero_line_list=sero_line_list,case_line_list=case_line_list,outbreak_line_list=outbreak_line_list)


  return(input_data_new)
}
#-------------------------------------------------------------------------------
#' @title mcmc_prelim_fit2
#'
#' @description Test multiple sets of parameters randomly drawn from range between maximum and minimum
#' values in order to find approximate values giving maximum likelihood
#'
#' @details This function is used to estimate the model parameter values giving maximum likelihood; it is primarily
#' intended to be used to generate initial parameter values for Markov Chain Monte Carlo fitting (using the mcmc()
#' function). [TBC]
#'
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
#' @param p_rep_severe Probability of observation of severe infection (set to NULL if being varied as a parameter)
#' @param p_rep_death Probability of observation of death (set to NULL if being varied as a parameter)
#' '
#' @export
#'
mcmc_prelim_fit2 <- function(n_iterations=1,n_param_sets=1,n_bounds=1,
                            type=NULL,pars_min=NULL,pars_max=NULL,input_data=list(),
                            obs_sero_data=list(),obs_case_data=list(),obs_outbreak_data=list(),
                            n_reps=1,mode_start=0,prior_type="zero",dt=1.0,enviro_data=NULL,R0_fixed_values=c(),
                            vaccine_efficacy=NULL,p_rep_severe=NULL,p_rep_death=NULL){

  #TODO - Add assertthat functions
  assert_that(length(pars_min)==length(pars_max))
  assert_that(type %in% c("FOI+R0","FOI","FOI+R0 enviro","FOI enviro"))
  assert_that(prior_type %in% c("zero","flat","exp","norm"))

  best_fit_results=list()
  n_params=length(pars_min)
  extra_params=c()
  if(is.null(vaccine_efficacy)==TRUE){extra_params=append(extra_params,"vaccine_efficacy")}
  if(is.null(p_rep_severe)==TRUE){extra_params=append(extra_params,"p_rep_severe")}
  if(is.null(p_rep_death)==TRUE){extra_params=append(extra_params,"p_rep_death")}
  param_names=create_param_labels(type,input_data,enviro_data,extra_params)
  #TODO - Additional assert_that checks
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
    const_list=list(type=type,pars_min=pars_min,pars_max=pars_max,n_reps=n_reps,mode_start=mode_start,
                    prior_type=prior_type,dt=dt,enviro_data=enviro_data,R0_fixed_values=R0_fixed_values,
                    vaccine_efficacy=vaccine_efficacy,p_rep_severe=p_rep_severe,p_rep_death=p_rep_death)

    for(set in 1:n_param_sets){
      if(set %% 10 == 0){cat("\n")}
      cat(set,"\t",sep="")
      param_prop=all_param_sets[set,]
      names(param_prop)=param_names
      like_prop=single_like_calc2(param_prop,input_data,obs_sero_data,obs_case_data,obs_outbreak_data,const_list)
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
#-------------------------------------------------------------------------------
#' @title mcmc_checks2
#'
#' @description Perform checks on MCMC inputs
#'
#' @details This function, which is called by MCMC2(), performs a number of checks on data to be used in fitting to
#' ensure proper functionality. It verifies that the number of parameters being fitted is consistent with other
#' settings, that certain values are not outwith sensible boundaries (e.g. probabilities must be between 0 and 1) and
#' [TBA]
#'
#' @param pars_ini = Initial parameter values
#' @param n_regions = Number of regions
#' @param type Type of parameter set (FOI only, FOI+R0, FOI and/or R0 coefficients associated with environmental
#'   covariates); choose from "FOI","FOI+R0","FOI enviro","FOI+R0 enviro"
#' @param pars_min Lower limits of parameter values if specified
#' @param pars_max Upper limits of parameter values if specified
#' @param prior_type Text indicating which type of calculation to use for prior probability
#'  If prior_type = "zero", prior probability is always zero
#'  If prior_type = "flat", prior probability is zero if FOI/R0 in designated ranges, -Inf otherwise
#'  If prior_type = "exp", prior probability is given by dexp calculation on FOI/R0 values
#'  If prior_type = "norm", prior probability is given by dnorm calculation on parameter values
#' @param enviro_data Values of environmental covariates (if in use)
#' @param R0_fixed_values Values of R0 to use if not being fitted (set to NULL if R0 is being fitted)
#' @param vaccine_efficacy Vaccine efficacy (set to NULL if being varied as a parameter)
#' @param p_rep_severe Probability of observation of severe infection (set to NULL if being varied as a parameter)
#' @param p_rep_death Probability of observation of death (set to NULL if being varied as a parameter)
#'
#' @export
#'
mcmc_checks2 <- function(pars_ini=c(),n_regions=1,type=NULL,pars_min=c(),pars_max=c(),prior_type=NULL,enviro_data=NULL,
                        R0_fixed_values=NULL,vaccine_efficacy=NULL,p_rep_severe=NULL,p_rep_death=NULL){

  param_names=names(pars_ini)
  n_params=length(pars_ini)
  assert_that(is.null(param_names)==FALSE) #Check that parameters have been named (should always be done in MCMC2())
  assert_that(type %in% c("FOI+R0","FOI","FOI+R0 enviro","FOI enviro")) #Check that type is one of those allowed
  assert_that(prior_type %in% c("zero","flat","exp","norm")) #Check that prior type is one of those allowed

  #If vaccine efficacy, severe case reporting probability and/or fatal case reporting probability have been set to NULL,
  # check that they appear among the parameters (should always have been set up in MCMC2() if create_param_labels()
  #ran correctly )
  if(is.numeric(vaccine_efficacy)==TRUE){
    flag_vacc_eff=0
    assert_that(0.0<=vaccine_efficacy && vaccine_efficacy<=1.0)
    } else {
    flag_vacc_eff=1
    assert_that("vaccine_efficacy" %in% param_names)
  }
  if(is.numeric(p_rep_severe)==TRUE){
    flag_severe=0
    assert_that(0.0<=p_rep_severe && p_rep_severe<=1.0)
    } else {
    flag_severe=1
    assert_that("p_rep_severe" %in% param_names)
  }
  if(is.numeric(p_rep_death)==TRUE){
    flag_death=0
    assert_that(0.0<=p_rep_death && p_rep_death<=1.0)
    } else {
    flag_death=1
    assert_that("p_rep_death" %in% param_names)
  }

  #If environmental data has been supplied, get names of variables
  if(is.null(enviro_data)==FALSE){
    env_vars=names(enviro_data[c(2:ncol(enviro_data))])
    n_env_vars=length(env_vars)
  }

  #Check that total number of parameters is correct based on "type" and on number of additional parameters (vaccine
  #efficacy, reporting probabilities)
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
    assert_that(n_params==(2*n_env_vars)+flag_vacc_eff+flag_severe+flag_death)
  }
  if(type=="FOI enviro"){
    assert_that(is.null(enviro_data)==FALSE)
    assert_that(is.null(R0_fixed_values)==FALSE)
    assert_that(length(R0_fixed_values)==n_regions)
    assert_that(n_params==n_env_vars+flag_vacc_eff+flag_severe+flag_death)
  }

  #TBA

  return(NULL)
}
