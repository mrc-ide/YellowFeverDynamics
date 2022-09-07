# R file for functions used for Markov Chain Monte Carlo fitting (and preliminary maximum-likelihood fitting) in
# YellowFeverDynamics package
#-------------------------------------------------------------------------------
#' @title MCMC
#'
#' @description Combined MCMC Multi-Region - series of MCMC steps for one or more regions
#'
#' @details This is the master function for running a Markov chain to optimize the parameters of the yellow fever model
#' based on the calculated likelihood of observing supplied data given a particular set of parameters.
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
MCMC <- function(pars_ini=c(),input_data=list(),obs_sero_data=NULL,obs_case_data=NULL,obs_outbreak_data=NULL,Niter=1,
                  type=NULL,pars_min=NULL,pars_max=NULL,n_reps=1,mode_start=0,prior_type="zero",dt=1.0,
                  enviro_data=NULL,R0_fixed_values=NULL,vaccine_efficacy=NULL,p_rep_severe=NULL,p_rep_death=NULL,
                  filename_prefix="Chain"){

  #Check that initial, minimum and maximum parameters are in vectors of same sizes
  n_params=length(pars_ini)
  assert_that(length(pars_min)==n_params)
  assert_that(length(pars_max)==n_params)

  extra_params=c()
  if(is.null(vaccine_efficacy)==TRUE){extra_params=append(extra_params,"vaccine_efficacy")}
  if(is.null(p_rep_severe)==TRUE){extra_params=append(extra_params,"p_rep_severe")}
  if(is.null(p_rep_death)==TRUE){extra_params=append(extra_params,"p_rep_death")}

  #Process input data to check that all regions with sero, case and/or outbreak data supplied are present, remove
  #regions without any supplied data, and add cross-referencing tables for use when calculating likelihood. Take
  #subset of environmental data (if used) and check that environmental data available for all regions
  input_data=input_data_process(input_data,obs_sero_data,obs_case_data,obs_outbreak_data)
  regions=names(table(input_data$region_labels)) #Regions in new processed input data list
  n_regions=length(regions)
  if(is.null(enviro_data)==FALSE){
    for(region in regions){assert_that(region %in% enviro_data$adm1)}
    enviro_data=subset(enviro_data,enviro_data$adm1 %in% regions)
    }

  #Label parameters according to order and fitting type
  param_names=create_param_labels(type,input_data,enviro_data,extra_params)
  names(pars_ini)=names(pars_min)=names(pars_max)=param_names

  #Run checks to ensure that number of parameters is correct for fitting type and number of regions/environmental
  #covariates
  checks<-mcmc_checks(pars_ini,n_regions,type,pars_min,pars_max,prior_type,enviro_data,
                       R0_fixed_values,vaccine_efficacy,p_rep_severe,p_rep_death)

  #Set up list of invariant parameter values to supply to other functions
  const_list=list(type=type,pars_min=pars_min,pars_max=pars_max,n_reps=n_reps,mode_start=mode_start,
                  prior_type=prior_type,dt=dt,enviro_data=enviro_data,R0_fixed_values=R0_fixed_values,
                  vaccine_efficacy=vaccine_efficacy,p_rep_severe=p_rep_severe,p_rep_death=p_rep_death)

  ### find posterior probability at start ###
  out = MCMC_step(param=pars_ini,input_data,obs_sero_data,obs_case_data,obs_outbreak_data,
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
      colnames(chain)=colnames(chain_prop)=names(pars_ini)
      for(i in 1:n_params){colnames(chain_prop)[i]=paste("Test_",colnames(chain_prop)[i],sep="")}
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

    #Decide whether next iteration will be adaptive
    if (iter>burnin & runif(1)<0.9){ #adapt
      adapt = 1
      chain_cov  = cov(chain[max(nrow(chain)-10000, 1):nrow(chain),])
    } else {
      adapt = 0
      chain_cov = 1
    }

    #Next iteration in chain
    out = MCMC_step(param,input_data,obs_sero_data,obs_case_data,obs_outbreak_data,chain_cov,adapt,like_current,
                     const_list)
  }

  #Get final parameter values
  param_out=exp(out$param)
  names(param_out)=names(pars_ini)

  return(param_out)
}
#-------------------------------------------------------------------------------
#' @title MCMC_step
#'
#' @description Single MCMC step - one or more regions
#'
#' @details This function runs a single step in a Markov chain set up using the function mcmc(). It proposes a
#' set of parameters using the param_prop_setup() function, calculates the likelihood of observing the observed data
#' based on that proposed parameter set, accepts or rejects the proposed parameter set based on the calculated
#' likelihood and existing chain information, then returns the next line of information for the chain to mcmc().
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
#' @param const_list = List of constant parameters/flags/etc. loaded to mcmc() (type,pars_min,pars_max,n_reps,
#' mode_start,prior_type,dt=dt,enviro_data,R0_fixed_values,vaccine_efficacy,p_rep_severe,p_rep_death)
#'
#' @export
#'
MCMC_step <- function(param=c(),input_data=list(),obs_sero_data=NULL,obs_case_data=NULL,obs_outbreak_data=NULL,
                       chain_cov=1,adapt=0,like_current=-Inf,const_list=list()) {

  #Propose new parameter values
  param_prop=param_prop_setup(param,chain_cov,adapt)

  #Calculate likelihood using single_like_calc function
  like_prop=single_like_calc(param_prop,input_data,obs_sero_data,obs_case_data,obs_outbreak_data,const_list)

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

  return(list(param=param,param_prop=param_prop,like_current=like_current,like_prop=like_prop,accept=accept))
}
#-------------------------------------------------------------------------------
#' @title single_like_calc
#'
#' @description Function which calculates and outputs likelihood of observing simulated data
#'
#' @details This function calculates the total likelihood of observing a set of observations (across multiple regions
#' and data types) for a given proposed parameter set.
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
#' @param const_list = List of constant parameters/flags/etc. loaded to mcmc() (type,pars_min,pars_max,n_reps,
#' mode_start,prior_type,dt=dt,enviro_data,R0_fixed_values,vaccine_efficacy,p_rep_severe,p_rep_death)
#'
#' @export
#'
single_like_calc <- function(param_prop=c(),input_data=list(),obs_sero_data=NULL,obs_case_data=NULL,
                              obs_outbreak_data=NULL,const_list=list()) {

  regions=input_data$region_labels
  n_regions=length(regions)
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

  if(is.null(obs_sero_data)){sero_like_values=NA}
  if(is.null(obs_case_data)){cases_like_values=deaths_like_values=NA}
  if(is.null(obs_outbreak_data)){outbreak_like_values=NA}

  ### If prior finite, evaluate likelihood ###
  if (is.finite(prior_prop)) {

    #Generate modelled data over all regions
    dataset <- Generate_Dataset(input_data,const_list$enviro_data,FOI_values,R0_values,
                                obs_sero_data,obs_case_data,obs_outbreak_data,
                                vaccine_efficacy,p_rep_severe,p_rep_death,
                                const_list$mode_start,const_list$n_reps,const_list$dt)

    #Likelihood of observing serological data
    if(is.null(obs_sero_data)==FALSE){
      sero_like_values=lgamma(obs_sero_data$samples+1)-lgamma(obs_sero_data$positives+1)-
        lgamma(obs_sero_data$samples-obs_sero_data$positives+1)+obs_sero_data$positives*log(dataset$model_sero_values)+
        (obs_sero_data$samples-obs_sero_data$positives)*log(1.0-dataset$model_sero_values)
    }
    #Likelihood of observing annual case/death data
    if(is.null(obs_case_data)==FALSE){
      model_case_values=dataset$model_case_values
      model_death_values=dataset$model_death_values
      for(i in 1:length(model_case_values)){
        model_case_values[i]=max(model_case_values[i],0.1)
        model_death_values[i]=max(model_death_values[i],0.1)
      }
      cases_like_values=dnbinom(x=obs_case_data$cases,mu=model_case_values,
                                size=rep(1,length(obs_case_data$cases)),log=TRUE)
      deaths_like_values=dnbinom(x=obs_case_data$deaths,mu=model_death_values,
                                 size=rep(1,length(obs_case_data$deaths)),log=TRUE)
    }
    #Likelihood of observing annual outbreak Y/N data
    if(is.null(obs_outbreak_data)==FALSE){
      outbreak_like_values=outbreak_risk_compare(model_outbreak_risk=dataset$model_outbreak_risk_values,
                                                 obs_data=obs_outbreak_data$outbreak_yn)
    }

    likelihood=prior_prop+mean(c(sum(sero_like_values,na.rm=TRUE),sum(cases_like_values,na.rm=TRUE),
                     sum(deaths_like_values,na.rm=TRUE),sum(outbreak_like_values,na.rm=TRUE)),na.rm=TRUE)

  } else {likelihood=-Inf}

  return(likelihood)
}
#-------------------------------------------------------------------------------
#' @title mcmc_prelim_fit
#'
#' @description Test multiple sets of parameters randomly drawn from range between maximum and minimum
#' values in order to find approximate values giving maximum likelihood
#'
#' @details This function is used to estimate the model parameter values giving maximum likelihood; it is primarily
#' intended to be used to generate initial parameter values for Markov Chain Monte Carlo fitting (using the mcmc()
#' function).
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
mcmc_prelim_fit <- function(n_iterations=1,n_param_sets=1,n_bounds=1,
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
      like_prop=single_like_calc(param_prop,input_data,obs_sero_data,obs_case_data,obs_outbreak_data,const_list)
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
#' @title mcmc_checks
#'
#' @description Perform checks on MCMC inputs
#'
#' @details This function, which is called by MCMC2(), performs a number of checks on data to be used in fitting to
#' ensure proper functionality. It verifies that the number of parameters being fitted is consistent with other
#' settings and that certain values are not outwith sensible boundaries (e.g. probabilities must be between 0 and 1).
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
mcmc_checks <- function(pars_ini=c(),n_regions=1,type=NULL,pars_min=c(),pars_max=c(),prior_type=NULL,enviro_data=NULL,
                        R0_fixed_values=NULL,vaccine_efficacy=NULL,p_rep_severe=NULL,p_rep_death=NULL){

  param_names=names(pars_ini)
  n_params=length(pars_ini)
  assert_that(is.null(param_names)==FALSE) #Check that parameters have been named (should always be done in MCMC2())
  assert_that(type %in% c("FOI+R0","FOI","FOI+R0 enviro","FOI enviro")) #Check that type is one of those allowed
  assert_that(prior_type %in% c("zero","flat","exp","norm")) #Check that prior type is one of those allowed

  #If vaccine efficacy, severe case reporting probability and/or fatal case reporting probability have been set to
  #NULL, check that they appear among the parameters (should always have been set up in MCMC() if
  #create_param_labels() ran correctly )
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

  return(NULL)
}
#-------------------------------------------------------------------------------
#' @title mcmc_FOI_R0_setup
#'
#' @description Set up FOI and R0 values and calculate some prior probability values for MCMC calculation
#'
#' @details Takes in parameter values used for Markov Chain Monte Carlo fitting, calculates spillover force of
#' infection and (optionally) reproduction number values either directly or from environmental covariates. Also
#' calculates related components of prior probability.
#'
#' @param type Type of parameter set (FOI only, FOI+R0, FOI and/or R0 coefficients associated with environmental
#'   covariates); choose from "FOI","FOI+R0","FOI enviro","FOI+R0 enviro"
#' @param prior_type Text indicating which type of calculation to use for prior probability
#'  If prior_type = "zero", prior probability is always zero
#'  If prior_type = "flat", prior probability is zero if FOI/R0 in designated ranges, -Inf otherwise
#'  If prior_type = "exp", prior probability is given by dexp calculation on FOI/R0 values
#'  If prior_type = "norm", prior probability is given by dnorm calculation on parameter values
#' @param regions Vector of region names
#' @param param_prop Proposed parameter values
#' @param enviro_data Environmental data frame, containing only relevant environmental variables
#' @param R0_fixed_values Values of R0 to use if not being fitted
#' @param pars_min Lower limits of parameter values if specified
#' @param pars_max Upper limits of parameter values if specified
#' '
#' @export
#'
mcmc_FOI_R0_setup <- function(type="",prior_type="",regions="",param_prop=c(),enviro_data=list(),R0_fixed_values=c(),
                              pars_min=c(),pars_max=c()){

  n_params=length(param_prop)
  n_regions=length(regions)
  if(type %in% c("FOI+R0 enviro","FOI enviro")){n_env_vars=ncol(enviro_data)-1}
  FOI_values=R0_values=rep(0,n_regions)

  if(type %in% c("FOI+R0 enviro","FOI enviro")){
    for(i in 1:n_regions){
      model_params=param_calc_enviro(param=param_prop,enviro_data=enviro_data[enviro_data$adm1==regions[i],])
      FOI_values[i]=model_params$FOI
      if(type=="FOI+R0 enviro"){R0_values[i]=model_params$R0} else {R0_values[i]=R0_fixed_values[i]}
    }
  }
  if(type %in% c("FOI+R0","FOI")){
    FOI_values=exp(param_prop[c(1:n_regions)])
    if(type=="FOI+R0"){R0_values=exp(param_prop[c((n_regions+1):(2*n_regions))])
    } else {R0_values=R0_fixed_values}
  }

  prior=0
  if(prior_type=="exp"){
    prior_FOI=dexp(FOI_values,rate=1,log=TRUE)
    if(type %in% c("FOI+R0","FOI+R0 enviro")){prior_R0=dexp(R0_values,rate=1,log=TRUE)} else {prior_R0=0}
    prior = prior+sum(prior_FOI)+sum(prior_R0)
  }
  if(prior_type=="flat"){
    if(is.null(pars_min)==FALSE){
      for(i in 1:n_params){
        if(param_prop[i]<pars_min[i]){prior=-Inf}
      }
    }
    if(is.null(pars_max)==FALSE){
      for(i in 1:n_params){
        if(param_prop[i]>pars_max[i]){prior=-Inf}
      }
    }
  }
  if(prior_type=="norm"){
    if(type=="FOI"){n_params_check=n_regions}
    if(type=="FOI+R0"){n_params_check=2*n_regions}
    if(type=="FOI enviro"){n_params_check=n_env_vars}
    if(type=="FOI+R0 enviro"){n_params_check=2*n_env_vars}
    prior=sum(dnorm(param_prop[c(1:n_params_check)],mean = 0,sd = 30,log = TRUE))
  }

  output=list(FOI_values=FOI_values,R0_values=R0_values,prior=prior)
  return(output)
}
#-------------------------------------------------------------------------------
#' @title param_prop_setup
#'
#' @description Set up proposed new parameter values for next step in chain
#'
#' @details Takes in current values of parameter set used for Markov Chain Monte Carlo fitting and proposes new values
#' from multivariate normal distribution where the existing values form the mean and the standard deviation is
#' based on the chain covariance or (if the flag "adapt" is set to 1) a flat value based on the number of parameters.
#'
#' @param param Previous parameter values used as input
#' @param chain_cov Covariance calculated from previous steps in chain
#' @param adapt 0/1 flag indicating which type of calculation to use for proposition value
#' '
#' @export
#'
param_prop_setup <- function(param=c(),chain_cov=1,adapt=0){

  n_params = length(param)
  if (adapt==1) {
    sigma = (2.38 ^ 2) * chain_cov / n_params #'optimal' scaling of chain covariance
    param_prop_a = rmvnorm(n = 1, mean = param, sigma = sigma)
  } else {
    sigma = ((1e-2) ^ 2) * diag(n_params) / n_params #this is an inital proposal covariance, see [Mckinley et al 2014]
    param_prop_a = rmvnorm(n = 1, mean = param, sigma = sigma)
  }
  param_prop = param_prop_a[1,]
  names(param_prop)=names(param)

  return(param_prop)
}
#-------------------------------------------------------------------------------
#' @title calc_like0
#'
#' @description Calculate "best case" likelihood
#'
#' @details Calculate the likelihood value which would be obtained (with no prior) if modelled data values were
#'              all equal to observed values (for use as comparison with MCMC likelihood values)
#'
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param obs_outbreak_data Outbreak Y/N data for comparison, by region and year, in format 0 = no outbreaks,
#'   1 = 1 or more outbreak(s)
#'
#' @export
#'
calc_like0 <- function(obs_sero_data=NULL,obs_case_data=NULL,obs_outbreak_data=NULL) {

  if(is.null(obs_sero_data)){sero_like_values=NA}
  if(is.null(obs_case_data)){cases_like_values=deaths_like_values=NA}
  if(is.null(obs_outbreak_data)){outbreak_like_values=NA}

  #Likelihood of observing serological data
  if(is.null(obs_sero_data)==FALSE){
    model_sero_values=obs_sero_data$positives/obs_sero_data$samples
    sero_like_values=lgamma(obs_sero_data$samples+1)-lgamma(obs_sero_data$positives+1)-
      lgamma(obs_sero_data$samples-obs_sero_data$positives+1)+obs_sero_data$positives*log(model_sero_values)+
      (obs_sero_data$samples-obs_sero_data$positives)*log(1.0-model_sero_values)
  }
  #Likelihood of observing annual case/death data
  if(is.null(obs_case_data)==FALSE){
    model_case_values=obs_case_data$cases
    model_death_values=obs_case_data$deaths
    for(i in 1:length(model_case_values)){
      model_case_values[i]=max(model_case_values[i],0.1)
      model_death_values[i]=max(model_death_values[i],0.1)
    }
    cases_like_values=dnbinom(x=obs_case_data$cases,mu=model_case_values,
                              size=rep(1,length(obs_case_data$cases)),log=TRUE)
    deaths_like_values=dnbinom(x=obs_case_data$deaths,mu=model_death_values,
                               size=rep(1,length(obs_case_data$deaths)),log=TRUE)
  }
  #Likelihood of observing annual outbreak Y/N data
  if(is.null(obs_outbreak_data)==FALSE){
    outbreak_like_values=outbreak_risk_compare(model_outbreak_risk=obs_outbreak_data$outbreak_yn,
                                               obs_data=obs_outbreak_data$outbreak_yn)
  }

  likelihood=mean(c(sum(sero_like_values,na.rm=TRUE),sum(cases_like_values,na.rm=TRUE),
                    sum(deaths_like_values,na.rm=TRUE),sum(outbreak_like_values,na.rm=TRUE)),na.rm=TRUE)

  return(likelihood)
}
