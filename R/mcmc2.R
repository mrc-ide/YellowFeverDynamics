# R file for functions used for Markov Chain Monte Carlo fitting (and preliminary maximum-likelihood fitting) - additional functionality
#-------------------------------------------------------------------------------
#' @title MCMC2
#'
#' @description Combined MCMC Multi-Region - series of MCMC steps for one or more regions
#'
#' @details This is the master function for running a Markov chain to optimize the parameters of the yellow fever
#' model based on the calculated likelihood of observing supplied data given a particular set of parameters.
#'
#' @param log_params_ini Initial values of parameters to be estimated. These should always be the log() values of the
#'   actual parameters, ordered as follows:
#'   1) Parameters controlling the value of spillover force of infection FOI, either a) a number of FOI values equal
#'   to the total number of regions to be considered or b) a number of environmental coefficients used to calculate
#'   FOI values from environmental covariates equal to the number of environmental covariates listed in the
#'   enviro_data frame. Values should be in alphabetical order by region in case (a) or in the order of the columns
#'   in the environmental data frame in case (b).
#'   2) If type is set to "FOI+R0" or "FOI+R0 enviro", parameters controlling the value of R0, either a) a number of R0 values equal to
#'   the total number of regions to be considered or b) a number of environmental coefficients used to calculate R0
#'   values from environmental covariates equal to the number of environmental covariates listed in the enviro_data
#'   frame. Values should be in alphabetical order by region in case (a) or in the order of the columns in the
#'   environmental data frame in case (b).
#'   3) Values of the additional parameters (vaccine efficacy vaccine_efficacy, severe case reporting probability
#'   p_rep_severe and fatal case reporting probability p_rep_death) if these are to be estimated, in the order
#'   vaccine_efficacy->p_rep_severe->p_rep_death. If these parameters are to be estimated, the values separately
#'   supplied to this function (see below) should be set to NULL, the default.
#' @param input_data List of population and vaccination data for multiple regions (created using data input creation
#'   code and usually loaded from RDS file)
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no.
#'   cases/no. deaths
#' @param filename_prefix Prefix of names for output files
#' @param Niter Total number of steps to run
#' @param type Type of parameter set (FOI only, FOI+R0, FOI and/or R0 coefficients associated with environmental
#'   covariates); choose from "FOI","FOI+R0","FOI enviro","FOI+R0 enviro"
#' @param prior_settings List containing settings for priors: must contain text named "type":
#'  If type = "zero", prior probability is always zero
#'  If type = "flat", prior probability is zero if log parameter values in designated ranges log_params_min and log_params_max,
#'   -Inf otherwise; log_params_min and log_params_max included in prior_settings as vectors of same length as log_params_ini
#'  If type = "exp", prior probability is given by dexp calculation on FOI/R0 values
#'  If type = "norm", prior probability is given by dnorm calculation on parameter values with settings based on vectors of values
#'   in prior_settings; norm_params_mean and norm_params_sd (mean and standard deviation values applied to log FOI/R0
#'   parameters and to actual values of additional parameters) + R0_mean + R0_sd (mean + standard deviation of computed R0, single values)
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start = 0, only vaccinated individuals
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (uniform by age, R0 based only)
#'  If mode_start = 3, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age)
#' @param dt time increment in days (must be 1 or 5)
#' @param n_reps Number of times to repeat calculations to get average likelihood at each step
#' @param enviro_data Data frame containing values of environmental covariates; set to NULL if not in use
#' @param R0_fixed_values Values of R0 to use if only FOI is subject to fitting (i.e. type set to "FOI" or "FOI
#'   enviro"); set to NULL if not in use
#' @param p_severe_inf Probability of an infection being severe
#' @param p_death_severe_inf Probability of a severe infection resulting in death
#' @param add_values List of parameters in addition to those governing FOI/R0, either giving a fixed value or giving NA to
#'   indicate that they are part of the fitted  parameter set
#'  vaccine_efficacy Vaccine efficacy (proportion of reported vaccinations causing immunity) (must be present)
#'  p_rep_severe Probability of observation of severe infection
#'  p_rep_death Probability of observation of death
#'  m_FOI_Brazil Multiplier of spillover FOI for Brazil regions (only relevant if regions in Brazil to be considered)
#'  p_rep_severe_af Probability of observation of severe infection in African regions
#'  p_rep_death_af Probability of observation of death in African regions
#'  p_rep_severe_sa Probability of observation of severe infection in South American regions
#'  p_rep_death_sa Probability of observation of death in South American regions
#'  NOTE: Either p_rep_severe and p_rep_death OR p_rep_severe_af, p_rep_death_af, p_rep_severe_sa and p_rep_death_sa must be used
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' @param mode_parallel TRUE/FALSE - indicate whether to use parallel processing on supplied cluster for speed
#' @param cluster Cluster of threads to use if mode_parallel = TRUE
#' '
#' @export
#'
MCMC2 <- function(log_params_ini = c(),input_data = list(),obs_sero_data = NULL,obs_case_data = NULL,filename_prefix = "Chain",
                  Niter = 1,type = NULL,mode_start = 0,prior_settings = list(type = "zero"),dt = 1.0,n_reps = 1,enviro_data = NULL,
                  R0_fixed_values = NULL,p_severe_inf = 0.12,p_death_severe_inf = 0.39,
                  add_values = list(vaccine_efficacy = 1.0,p_rep_severe = 1.0,p_rep_death = 1.0,m_FOI_Brazil = 1.0),
                  deterministic = FALSE,mode_parallel = FALSE,cluster = NULL){

  assert_that(is.logical(deterministic))
  assert_that(mode_start %in% c(0,1,3),msg = "mode_start must have value 0, 1 or 3")
  n_params = length(log_params_ini)

  extra_estimated_params = c()
  for(var_name in names(add_values)){
    if(is.null(add_values[[var_name]]) ==TRUE){extra_estimated_params = append(extra_estimated_params,var_name)}
  }

  #Process input data to check that all regions with sero and/or case data supplied are present, remove
  #regions without any supplied data, and add cross-referencing tables for use when calculating likelihood. Take
  #subset of environmental data (if used) and check that environmental data available for all regions
  input_data = input_data_process(input_data,obs_sero_data,obs_case_data)
  regions = names(table(input_data$region_labels)) #Regions in new processed input data list
  n_regions = length(regions)
  if(is.null(enviro_data) ==FALSE){
    for(region in regions){assert_that(region %in% enviro_data$region)}
    enviro_data = subset(enviro_data,enviro_data$region %in% regions)
  }

  #Label parameters according to order and fitting type
  param_names = create_param_labels2(type,input_data,enviro_data,extra_estimated_params)
  names(log_params_ini) = param_names

  #Run checks on inputs
  checks<-mcmc_checks2(log_params_ini,n_regions,type,prior_settings,enviro_data,R0_fixed_values,add_values,
                       extra_estimated_params)
  if(prior_settings$type =="flat"){names(prior_settings$log_params_min) = names(prior_settings$log_params_max) = param_names}

  #Set up list of invariant parameter values to supply to other functions
  consts = list(type = type,mode_start = mode_start,prior_settings = prior_settings,dt = dt,n_reps = n_reps,enviro_data = enviro_data,
                R0_fixed_values = R0_fixed_values,p_severe_inf = p_severe_inf,p_death_severe_inf = p_death_severe_inf,
                add_values = add_values,extra_estimated_params = extra_estimated_params,
                deterministic = deterministic,mode_parallel = mode_parallel,cluster = cluster)

  #MCMC setup
  chain = chain_prop = posterior_current = posterior_prop = flag_accept = chain_cov_all = NULL
  burnin = min(2*n_params, Niter)
  fileIndex = 0
  log_params = log_params_ini
  chain_cov = 1
  adapt = 0
  posterior_value_current = -Inf

  #Iterative estimation
  for (iter in 1:Niter){

    #Propose new parameter values
    log_params_prop = param_prop_setup(log_params,chain_cov,adapt)

    #Calculate likelihood using single_posterior_calc2 function
    posterior_value_prop = single_posterior_calc2(log_params_prop,input_data,obs_sero_data,obs_case_data,consts)
    gc()

    if(is.finite(posterior_value_prop) ==FALSE) {
      p_accept = -Inf
    } else {
      p_accept = posterior_value_prop - posterior_value_current
      if(is.na(p_accept) ){ p_accept = -Inf}
    }

    ## accept/reject step:
    tmp = runif(1)
    if(tmp<min(exp(p_accept),1)) {
      log_params = log_params_prop
      posterior_value_current = posterior_value_prop
      accept = 1
    } else {accept = 0}

    #save current step
    chain = rbind(chain, log_params)
    chain_prop = rbind(chain_prop,log_params_prop)
    posterior_current = rbind(posterior_current,posterior_value_current)
    posterior_prop = rbind(posterior_prop,posterior_value_prop)
    flag_accept = rbind(flag_accept, accept)
    chain_cov_all = rbind(chain_cov_all,max(chain_cov))

    #Set output headings
    if(iter ==1){
      colnames(chain) = colnames(chain_prop) = names(log_params_ini)
      for(i in 1:n_params){colnames(chain_prop)[i] = paste("Test_",colnames(chain_prop)[i],sep = "")}
      colnames(posterior_current) = "posterior_current"
      colnames(posterior_prop) = "posterior_prop"
      colnames(flag_accept) = "flag_accept"
      colnames(chain_cov_all) = "chain_cov_all"
    }

    #Output chain to file every 10 iterations; start new file every 10,000 iterations
    if (iter %% 10 == 0){
      if (iter %% 10000 == 10){fileIndex = (iter-10)/10000}

      if(fileIndex >= 10){fn = paste0(filename_prefix, fileIndex, ".csv")} else {fn = paste0(filename_prefix, "0", fileIndex, ".csv")}
      if(file.exists(fn) == FALSE){file.create(fn)}
      lines = min(((fileIndex*10000)+1), iter):iter

      data_out<-cbind(posterior_current, posterior_prop, exp(chain), flag_accept, exp(chain_prop), chain_cov_all)[lines, ]
      if(fileAccess(fn, 2) == 0){write.csv(data_out, fn, row.names = FALSE)}
    }

    #Decide whether next iteration will be adaptive
    if (iter>burnin & runif(1)<0.9){ #adapt
      adapt = 1
      chain_cov  = cov(chain[max(nrow(chain)-10000, 1):nrow(chain),])
    } else {
      adapt = 0
      chain_cov = 1
    }
  }

  #Get final parameter values
  param_out = exp(log_params)
  names(param_out) = names(log_params_ini)

  return(param_out)
}
#-------------------------------------------------------------------------------
#' @title single_posterior_calc2
#'
#' @description Function which calculates and outputs posterior likelihood of observing simulated data
#'
#' @details This function calculates the posterior likelihood of observing a set of observations (across multiple
#' regions and data types) for a given proposed parameter set. [TBA]
#'
#' @param log_params_prop Proposed values of parameters (natural logarithm of actual parameters)
#' @param input_data List of population and vaccination data for multiple regions (created using data input
#'   creation code and usually loaded from RDS file), with cross-reference tables added using input_data_process
#'   in MCMC
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param consts = List of constant parameters/flags/etc. loaded to mcmc() (type,
#'   mode_start,prior_settings,dt,n_reps,enviro_data,R0_fixed_values,p_severe_inf,
#'   p_death_severe_inf,add_values list,extra_estimated_params,deterministic, mode_parallel, cluster)
#'
#' @export
#'
single_posterior_calc2 <- function(log_params_prop = c(),input_data = list(),obs_sero_data = NULL,obs_case_data = NULL,
                                   consts = list()){

  #Get additional values and calculate associated priors
  vaccine_efficacy = p_rep_severe = p_rep_death = p_rep_severe_af = p_rep_death_af = p_rep_severe_sa = p_rep_death_sa = m_FOI_Brazil = 1.0
  prior_add = 0
  for(var_name in names(consts$add_values)){
    if(is.numeric(consts$add_values[[var_name]]) ==FALSE){
      i = match(var_name,names(log_params_prop))
      value = exp(as.numeric(log_params_prop[i]))
      assign(var_name,value)
      if(consts$prior_settings$type =="norm"){
        prior_add = prior_add+log(dtrunc(value,"norm",a = 0,b = 1,mean = consts$prior_settings$norm_params_mean[i],
                                         sd = consts$prior_settings$norm_params_sd[i]))
      } else {
        if(consts$prior_settings$type =="flat"){
          if(value<consts$prior_settings$log_params_min[i] || value>consts$prior_settings$log_params_max[i]){prior_add = -Inf}
        }
      }
    } else {assign(var_name,consts$add_values[[var_name]])}
  }

  #If additional values give finite prior, get FOI and R0 values and calculate associated prior
  if(is.finite(prior_add)){
    regions = input_data$region_labels
    n_regions = length(regions)

    FOI_R0_data = mcmc_FOI_R0_setup(consts$type,consts$prior_settings,regions,log_params_prop,
                                    consts$enviro_data,consts$R0_fixed_values)
    FOI_values = FOI_R0_data$FOI_values

    for(n_region in 1:n_regions){if(substr(regions[n_region],1,3) =="BRA"){FOI_values[n_region] = FOI_values[n_region]*m_FOI_Brazil}}
    R0_values = FOI_R0_data$R0_values
    if(consts$prior_settings$type =="norm"){
      prior_prop = FOI_R0_data$prior + prior_add +
        sum(log(dtrunc(R0_values,"norm",a = 0,b = Inf,mean = consts$prior_settings$R0_mean,sd = consts$prior_settings$R0_sd))) +
        sum(log(dtrunc(FOI_values,"norm",a = 0,b = 1,mean = consts$prior_settings$FOI_mean,sd = consts$prior_settings$FOI_sd)))
    } else {
      prior_prop = FOI_R0_data$prior+prior_add
    }
  } else {prior_prop = -Inf}

  ### If prior finite, evaluate likelihood ###
  if (is.finite(prior_prop)) {

    #Generate modelled data over all regions
    if("p_rep_severe" %in% names(consts$add_values)){
      dataset <- Generate_Dataset(input_data,FOI_values,R0_values,obs_sero_data,obs_case_data,vaccine_efficacy,
                                  consts$p_severe_inf,consts$p_death_severe_inf,p_rep_severe,p_rep_death,
                                  consts$mode_start,start_SEIRV = NULL,consts$dt,consts$n_reps,consts$deterministic,
                                  consts$mode_parallel,consts$cluster)
    } else {
      dataset <- Generate_Dataset2(input_data,FOI_values,R0_values,obs_sero_data,obs_case_data,vaccine_efficacy,
                                   consts$p_severe_inf,consts$p_death_severe_inf,
                                   p_rep_severe_af,p_rep_death_af,p_rep_severe_sa,p_rep_death_sa,
                                   consts$mode_start,start_SEIRV = NULL,consts$dt,consts$n_reps,consts$deterministic,
                                   consts$mode_parallel,consts$cluster)

    }

    #Likelihood of observing serological data
    if(is.null(obs_sero_data) ==FALSE){
      sero_like_values = sero_data_compare(dataset$model_sero_values,obs_sero_data)
    } else {sero_like_values = 0}

    #Likelihood of observing annual case/death data
    if(is.null(obs_case_data) ==FALSE){
      cases_like_values = case_data_compare(dataset$model_case_values,obs_case_data$cases)
      if(is.null(obs_case_data$deaths) ==FALSE){
        deaths_like_values = case_data_compare(dataset$model_death_values,obs_case_data$deaths)
      } else {deaths_like_values = 0}
    } else {cases_like_values = deaths_like_values = 0}

    # posterior = prior_prop+mean(c(sum(sero_like_values,na.rm = TRUE),sum(cases_like_values,na.rm = TRUE),
    #                             sum(deaths_like_values,na.rm = TRUE)),na.rm = TRUE)
    posterior = prior_prop+sum(sero_like_values,na.rm = TRUE)+sum(cases_like_values,na.rm = TRUE)+sum(deaths_like_values,na.rm = TRUE)

  } else {posterior = -Inf}

  return(posterior)
}
#-------------------------------------------------------------------------------
#' @title mcmc_checks2
#'
#' @description Perform checks on MCMC inputs
#'
#' @details This function, which is called by MCMC(), performs a number of checks on data to be used in fitting to
#' ensure proper functionality. It verifies that the number of parameters being estimated is consistent with other
#' settings and that certain values are not outwith sensible boundaries (e.g. probabilities must be between 0 and 1).
#'
#' @param log_params_ini Initial values of parameters (natural logarithm of actual parameters)
#' @param n_regions Number of regions
#' @param type Type of parameter set (FOI only, FOI+R0, FOI and/or R0 coefficients associated with environmental
#'   covariates); choose from "FOI","FOI+R0","FOI enviro","FOI+R0 enviro"
#' @param prior_settings TBA
#' @param enviro_data Values of environmental covariates (if in use)
#' @param R0_fixed_values Values of R0 to use if type set to "FOI" or "FOI enviro"
#' @param add_values TBA
#' @param extra_estimated_params TBA
#'
#' @export
#'
mcmc_checks2 <- function(log_params_ini = c(),n_regions = 1,type = NULL,prior_settings = list(type = "zero"),enviro_data = NULL,
                         R0_fixed_values = NULL,
                         add_values = list(vaccine_efficacy = 1.0,p_rep_severe = 1.0,p_rep_death = 1.0,m_FOI_Brazil = 1.0),
                         extra_estimated_params = list()){

  param_names = names(log_params_ini)
  n_params = length(log_params_ini)
  assert_that(is.null(param_names) ==FALSE,msg = "Parameters should be named using create_param_labels2")
  assert_that(type %in% c("FOI+R0","FOI","FOI+R0 enviro","FOI enviro"),msg = "Check type parameter")
  assert_that(prior_settings$type %in% c("zero","flat","exp","norm"),msg = "Check prior_settings$type")
  if(prior_settings$type =="flat"){
    assert_that(length(prior_settings$log_params_min) ==n_params)
    assert_that(length(prior_settings$log_params_max) ==n_params)
  }
  if(prior_settings$type =="norm"){
    assert_that(length(prior_settings$norm_params_mean) ==n_params)
    assert_that(length(prior_settings$norm_params_sd) ==n_params)
    assert_that(is.numeric(prior_settings$R0_mean))
    assert_that(is.numeric(prior_settings$R0_sd))
    assert_that(is.numeric(prior_settings$FOI_mean))
    assert_that(is.numeric(prior_settings$FOI_sd))
  }

  # Check additional values
  add_value_names = names(add_values) #TODO: formalize flexibility
  # assert_that(all(add_value_names ==c("vaccine_efficacy","p_rep_severe","p_rep_death","m_FOI_Brazil")) ||
  #               all(add_value_names ==c("vaccine_efficacy","p_rep_severe_af","p_rep_death_af",
  #                                      "p_rep_severe_sa","p_rep_death_sa","m_FOI_Brazil")))
  assert_that(all(extra_estimated_params %in% add_value_names))
  for(var_name in add_value_names){
    if(var_name %in% extra_estimated_params){
      assert_that(is.null(add_values[[var_name]]))
    } else {assert_that(add_values[[var_name]] <= 1.0 && add_values[[var_name]] >= 0.0)}
  }

  # If environmental data has been supplied, get names of variables
  if(is.null(enviro_data) ==FALSE){
    env_vars = names(enviro_data[c(2:ncol(enviro_data))])
    n_env_vars = length(env_vars)
  } else {assert_that(type %in% c("FOI+R0","FOI"))}

  # Check that total number of parameters is correct based on "type" and on number of additional parameters (vaccine
  # efficacy, reporting probabilities); check parameters named in correct order (TBA); all should be correct if parameter
  # names created using create_param_labels2
  if(type =="FOI+R0"){
    assert_that(n_params ==(2*n_regions)+length(extra_estimated_params))
  }
  if(type =="FOI"){
    assert_that(n_params ==n_regions+length(extra_estimated_params))
    assert_that(is.null(R0_fixed_values) ==FALSE)
    assert_that(length(R0_fixed_values) ==n_regions)
  }
  if(type =="FOI+R0 enviro"){
    assert_that(is.null(enviro_data) ==FALSE)
    assert_that(n_params ==(2*n_env_vars)+length(extra_estimated_params))
    for(i in 1:n_env_vars){
      assert_that(param_names[i] ==paste("FOI_",env_vars[i],sep = ""))
      assert_that(param_names[i+n_env_vars] ==paste("R0_",env_vars[i],sep = ""))
    }
  }
  if(type =="FOI enviro"){
    assert_that(is.null(enviro_data) ==FALSE)
    assert_that(n_params ==n_env_vars+length(extra_estimated_params))
    assert_that(is.null(R0_fixed_values) ==FALSE)
    assert_that(length(R0_fixed_values) ==n_regions)
    for(i in 1:n_env_vars){
      assert_that(param_names[i] ==paste("FOI_",env_vars[i],sep = ""))
    }
  }

  return(NULL)
}
#-------------------------------------------------------------------------------
#' @title param_prop_setup
#'
#' @description Set up proposed new log parameter values for next step in chain
#'
#' @details Takes in current values of parameter set used for Markov Chain Monte Carlo fitting and proposes new values
#' from multivariate normal distribution where the existing values form the mean and the standard deviation is
#' based on the chain covariance or (if the flag "adapt" is set to 1) a flat value based on the number of parameters.
#'
#' @param log_params Previous log parameter values used as input
#' @param chain_cov Covariance calculated from previous steps in chain
#' @param adapt 0/1 flag indicating which type of calculation to use for proposition value
#' '
#' @export
#'
param_prop_setup <- function(log_params = c(),chain_cov = 1,adapt = 0){

  n_params = length(log_params)
  if (adapt ==1) {
    sigma = (2.38 ^ 2) * chain_cov / n_params #'optimal' scaling of chain covariance
    log_params_prop_a = rmvnorm(n = 1, mean = log_params, sigma = sigma)
  } else {
    sigma = ((1e-2) ^ 2) * diag(n_params) / n_params #this is an inital proposal covariance, see [Mckinley et al 2014]
    log_params_prop_a = rmvnorm(n = 1, mean = log_params, sigma = sigma)
  }
  log_params_prop = log_params_prop_a[1,]

  return(log_params_prop)
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
#' @param prior_settings TBA
#' @param regions Vector of region names
#' @param log_params_prop Proposed values of parameters (natural logarithm of actual parameters)
#' @param enviro_data Environmental data frame, containing only relevant environmental variables
#' @param R0_fixed_values Values of R0 to use if type set to "FOI" or "FOI enviro"
#' '
#' @export
#'
mcmc_FOI_R0_setup <- function(type = "",prior_settings = list(type = "zero"),regions = "",log_params_prop = c(),enviro_data = list(),
                              R0_fixed_values = c()){

  n_regions = length(regions)
  FOI_values = R0_values = rep(0,n_regions)

  if(type %in% c("FOI+R0 enviro","FOI enviro")){
    n_env_vars = ncol(enviro_data)-1
    if(type =="FOI+R0 enviro"){n_values = 2*n_env_vars} else {n_values = n_env_vars}
    enviro_coeffs = exp(log_params_prop[c(1:n_values)])

    for(i in 1:n_regions){
      model_params = param_calc_enviro(enviro_coeffs,
                                       as.numeric(enviro_data[enviro_data$region ==regions[i],1+c(1:n_env_vars)]))
      FOI_values[i] = model_params$FOI
      if(type =="FOI+R0 enviro"){R0_values[i] = model_params$R0} else {R0_values[i] = R0_fixed_values[i]}
    }
  } else {
    FOI_values = exp(log_params_prop[c(1:n_regions)])
    if(type =="FOI+R0"){
      n_values = 2*n_regions
      R0_values = exp(log_params_prop[c((n_regions+1):n_values)])
    } else {
      n_values = n_regions
      R0_values = R0_fixed_values
    }
  }

  prior = 0
  if(prior_settings$type =="norm"){
    prior = sum(dnorm(log_params_prop[c(1:n_values)],mean = prior_settings$norm_params_mean[c(1:n_values)],
                      sd = prior_settings$norm_params_sd[c(1:n_values)],log = TRUE))
  } else {
    if(prior_settings$type =="flat"){
      for(i in 1:n_values){
        if(log_params_prop[i]<prior_settings$log_params_min[i]){prior = -Inf}
        if(log_params_prop[i]>prior_settings$log_params_max[i]){prior = -Inf}
      }
    } else {
      if(prior_settings$type =="exp"){
        prior_FOI = dexp(FOI_values,rate = 1,log = TRUE)
        if(type %in% c("FOI+R0","FOI+R0 enviro")){prior_R0 = dexp(R0_values,rate = 1,log = TRUE)} else {prior_R0 = 0}
        prior = prior+sum(prior_FOI)+sum(prior_R0)
      }
    }
  }

  return(list(FOI_values = FOI_values,R0_values = R0_values,prior = prior))
}
#-------------------------------------------------------------------------------
#' @title mcmc_prelim_fit2
#'
#' @description Test multiple sets of parameters randomly drawn from range between maximum and minimum
#' values in order to find approximate values giving maximum posterior likelihood
#'
#' @details This function is used to estimate the model parameter values giving maximum posterior likelihood; it is
#' primarily intended to be used to generate initial parameter values for Markov Chain Monte Carlo fitting (using
#' the mcmc() function).
#'
#' @param n_iterations = Number of times to run and adjust maximum/minimum
#' @param n_param_sets = Number of parameter sets to run in each iteration
#' @param n_bounds = Number of parameter sets (with highest likelihood values) to take at each iteration to create new
#' maximum/minimum values
#' @param type Type of parameter set (FOI only, FOI+R0, FOI and/or R0 coefficients associated with environmental
#'   covariates); choose from "FOI","FOI+R0","FOI enviro","FOI+R0 enviro"
#' @param log_params_min Initial lower limits of values of parameters being estimated (natural logarithm of actual limits)
#' @param log_params_max Initial upper limits of values of parameters being estimated (natural logarithm of actual limits)
#' @param input_data List of population and vaccination data for multiple regions (created using data input creation
#'   code and usually loaded from RDS file)
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start = 0, only vaccinated individuals
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (uniform by age, R0 based only)
#'  If mode_start = 3, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age)
#' @param prior_settings TBA
#' @param dt time increment in days (must be 1 or 5)
#' @param n_reps Number of repetitions
#' @param enviro_data Values of environmental variables (if in use)
#' @param R0_fixed_values Values of R0 to use if type set to "FOI" or "FOI enviro"
#' @param p_severe_inf TBA
#' @param p_death_severe_inf TBA
#' @param add_values List of parameters in addition to those governing FOI/R0, either giving a fixed value or giving NA to
#'   indicate that they are part of the fitted  parameter set
#'  vaccine_efficacy Vaccine efficacy (proportion of reported vaccinations causing immunity) (must be present)
#'  p_rep_severe Probability of observation of severe infection
#'  p_rep_death Probability of observation of death
#'  m_FOI_Brazil Multiplier of spillover FOI for Brazil regions (only relevant if regions in Brazil to be considered)
#'  p_rep_severe_af Probability of observation of severe infection in African regions
#'  p_rep_death_af Probability of observation of death in African regions
#'  p_rep_severe_sa Probability of observation of severe infection in South American regions
#'  p_rep_death_sa Probability of observation of death in South American regions
#'  NOTE: Either p_rep_severe and p_rep_death OR p_rep_severe_af, p_rep_death_af, p_rep_severe_sa and p_rep_death_sa must be used
#' @param deterministic TBA
#' @param mode_parallel TRUE/FALSE - indicate whether to use parallel processing on supplied cluster for speed
#' @param cluster Cluster of threads to use if multithreading to be used; set to NULL otherwise
#' '
#' @export
#'
mcmc_prelim_fit2 <- function(n_iterations = 1,n_param_sets = 1,n_bounds = 1,type = NULL,log_params_min = NULL,
                             log_params_max = NULL,input_data = list(),obs_sero_data = list(),obs_case_data = list(),
                             mode_start = 0,prior_settings = list(type = "zero"),dt = 1.0,n_reps = 1,enviro_data = NULL,R0_fixed_values = c(),
                             p_severe_inf = 0.12, p_death_severe_inf = 0.39,
                             add_values = list(vaccine_efficacy = 1.0,p_rep_severe = 1.0,p_rep_death = 1.0,m_FOI_Brazil = 1.0),
                             deterministic = TRUE,mode_parallel = FALSE,cluster = NULL){

  #TODO - Add assertthat functions
  assert_that(mode_start %in% c(0,1,3),msg = "mode_start must have value 0, 1 or 3")
  assert_that(length(log_params_min) ==length(log_params_max),msg = "Parameter limit vectors must have same lengths")
  assert_that(type %in% c("FOI+R0","FOI","FOI+R0 enviro","FOI enviro"))
  assert_that(prior_settings$type %in% c("zero","exp","norm"),msg = "Prior type must be 'zero', 'exp' or 'norm'")

  best_fit_results = list()
  n_params = length(log_params_min)

  #Get additional values
  extra_estimated_params = c()
  add_value_names = names(add_values)
  for(var_name in add_value_names){
    if(is.null(add_values[[var_name]]) ==TRUE){extra_estimated_params = append(extra_estimated_params,var_name)}
  }
  param_names = create_param_labels2(type,input_data,enviro_data,extra_estimated_params)

  #TODO - Additional assert_that checks
  assert_that(length(param_names) ==n_params)
  names(log_params_min) = names(log_params_max) = param_names
  xlabels = param_names
  for(i in 1:n_params){xlabels[i] = substr(xlabels[i],1,15)}
  ylabels = 10^c(-8,-6,-4,-3,-2,-1,0,1)
  par(mar = c(6,2,1,1))
  ylim = c(min(log_params_min),max(log_params_max))

  for(iteration in 1:n_iterations){
    cat("\nIteration: ",iteration,"\n",sep = "")
    all_param_sets <- lhs(n = n_param_sets,rect = cbind(log_params_min,log_params_max))
    results = data.frame()
    consts = list(type = type,mode_start = mode_start,prior_settings = prior_settings,
                  dt = dt,n_reps = n_reps,enviro_data = enviro_data,R0_fixed_values = R0_fixed_values,
                  p_severe_inf = p_severe_inf, p_death_severe_inf = p_death_severe_inf,add_values = add_values,
                  deterministic = deterministic,mode_parallel = mode_parallel,cluster = cluster)

    for(set in 1:n_param_sets){
      cat("\n\tSet: ",set,sep = "")
      log_params_prop = all_param_sets[set,]

      cat("\n\tParams: ",signif(log_params_prop,3))

      names(log_params_prop) = param_names
      posterior_value = single_posterior_calc2(log_params_prop,input_data,obs_sero_data,obs_case_data,consts)
      gc()
      results<-rbind(results,c(set,exp(log_params_prop),posterior_value))
      if(set ==1){colnames(results) = c("set",param_names,"posterior")}

      cat("\n\tPosterior likelihood = ",posterior_value,sep = "")

    }
    results<-results[order(results$posterior,decreasing = TRUE), ]
    best_fit_results[[iteration]] = results

    log_params_min_new = log_params_max_new = rep(0,n_params)
    for(i in 1:n_params){
      log_params_min_new[i] = min(log(results[c(1:n_bounds),i+1]))
      log_params_max_new[i] = max(log(results[c(1:n_bounds),i+1]))
    }
    names(log_params_min_new) = names(log_params_max_new) = param_names

    matplot(x = c(1:n_params),y = log(t(results[c(1:n_bounds),c(1:n_params)+1])),type = "p",pch = 16,col = 1,
            xaxt = "n",yaxt = "n",xlab = "",ylab = "",ylim = ylim)
    axis(side = 1,at = c(1:n_params),labels = xlabels,las = 2,cex.axis = 0.7)
    axis(side = 2,at = log(ylabels),labels = ylabels)
    matplot(x = c(1:n_params),y = log_params_min,type = "l",col = 1,lty = 2,add = TRUE)
    matplot(x = c(1:n_params),y = log_params_max,type = "l",col = 1,lty = 2,add = TRUE)
    matplot(x = c(1:n_params),y = log_params_min_new,type = "l",col = 2,add = TRUE)
    matplot(x = c(1:n_params),y = log_params_max_new,type = "l",col = 2,add = TRUE)

    log_params_min = log_params_min_new
    log_params_max = log_params_max_new
  }

  return(best_fit_results)
}
