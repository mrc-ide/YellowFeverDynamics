# R file for functions used for Markov Chain Monte Carlo fitting (and preliminary maximum-likelihood fitting) in
# YEP package
#-------------------------------------------------------------------------------
#' @title MCMC_VarFR
#'
#' @description Combined MCMC Multi-Region - series of MCMC iterations for one or more regions
#'
#' @details This is the master function for running a Markov chain to optimize the parameters of the yellow fever
#' model based on the calculated likelihood of observing supplied data given a particular set of parameters.
#'
#' @param log_params_ini Initial values of parameters to be estimated. These should always be the log() values of the
#'   actual parameters, ordered as follows: \cr
#'   1) A number of environmental coefficients used to calculate spillover force of infection values from environmental
#'   covariates equal to the total number of environmental covariates. Values should be in the
#'   order of the columns in enviro_data_const and then the covariates in enviro_data_var. \cr
#'   2) A number of environmental coefficients used to calculate basic reproduction number
#'   values from environmental covariates equal to the total number of environmental covariates.
#'   Values should be in the order of the columns in enviro_data_const and then the covariates in enviro_data_var. \cr
#'   3) Values of the additional parameters (reported vaccination effectiveness vaccine_efficacy, severe case reporting
#'   probability p_rep_severe, and fatal case reporting probability p_rep_death, and Brazil spillover FOI multiplier m_FOI_Brazil);
#'   if these are to be estimated, in the order vaccine_efficacy->p_rep_severe->p_rep_death->m_FOI_Brazil. If these parameters
#'   are to be estimated, the values separately supplied to this function (see add_values below) should be set to NA.
#' @param input_data List of population and vaccination data for multiple regions (created using data input creation
#'   code and usually loaded from RDS file)
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no.
#'   cases/no. deaths
#' @param filename_prefix Prefix of names for output files; function outputs a CSV file every 10,000 iterations with a name in the
#'  format: "(filename_prefix)XX.csv", e.g. Chain00.csv
#' @param Niter Total number of iterations to run
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination \cr
#'  If mode_start = 0, only vaccinated individuals \cr
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (uniform by age, R0 based only) \cr
#'  If mode_start = 3, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age)
#' @param prior_settings List containing settings for priors: must contain text named "type":
#'  If type = "zero", prior probability is always zero \cr
#'  If type = "flat", prior probability is zero if log parameter values in designated ranges param_min_limits and param_max_limits,
#'   -Inf otherwise; param_min_limits and param_max_limits included in prior_settings as vectors of same length as log_params_ini \cr
#'  If type = "norm", prior probability is given by truncated normal distribution calculation on parameter values with settings based
#'  on vectors of values in prior_settings: \cr
#'   norm_params_mean and norm_params_sd (vectors of mean and standard deviation values applied to log FOI/R0
#'   parameters and to actual values of additional parameters) \cr
#'   + FOI_mean + FOI_sd (mean + standard deviation of computed FOI, single values)  \cr
#'   + R0_mean + R0_sd (mean + standard deviation of computed R0, single values) \cr
#'   + param_min_limits and param_max_limits (lower and upper limits applied to truncated normal distributions)
#' @param dt time increment in days (must be 1 or 5)
#' @param n_reps Number of times to repeat calculations to get average likelihood at each iteration
#' @param enviro_data_const Data frame of values of constant environmental covariates (columns) by region (rows)
#' @param enviro_data_var List containing values of time-varying environmental covariates (TBA)
#' @param p_severe_inf Probability of an infection being severe
#' @param p_death_severe_inf Probability of a severe infection resulting in death
#' @param add_values List of parameters in addition to those governing FOI/R0, either giving a fixed value or giving NA to
#'   indicate that they are part of the fitted  parameter set \cr
#'  vaccine_efficacy Vaccine efficacy (proportion of reported vaccinations causing immunity) (must be present) \cr
#'  p_rep_severe Probability of observation of severe infection \cr
#'  p_rep_death Probability of observation of death \cr
#'  m_FOI_Brazil Multiplier of spillover FOI for Brazil regions (only relevant if regions in Brazil to be considered)
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' @param mode_time TBA
#' @param mode_parallel TRUE/FALSE - indicate whether to use parallel processing on supplied cluster for speed
#' @param cluster Cluster of threads to use if mode_parallel = TRUE
#' '
#' @export
#'
MCMC_VarFR <- function(log_params_ini = c(), input_data = list(), obs_sero_data = NULL, obs_case_data = NULL, filename_prefix = "Chain",
                       Niter = 1, mode_start = 0, prior_settings = list(type = "zero"), dt = 1.0, n_reps = 1,
                       enviro_data_const = list(), enviro_data_var = list(), p_severe_inf = 0.12, p_death_severe_inf = 0.39,
                       add_values = list(vaccine_efficacy = 1.0, p_rep_severe = 1.0, p_rep_death = 1.0, m_FOI_Brazil = 1.0),
                       deterministic = FALSE, mode_time = 1, mode_parallel = FALSE, cluster = NULL){

  assert_that(is.logical(deterministic))
  assert_that(mode_start %in% c(0, 1, 3), msg = "mode_start must have value 0, 1 or 3")
  if(is.null(obs_case_data)==FALSE){assert_that(all(c("p_rep_severe","p_rep_death") %in% names(add_values)),
                                                msg = "Reporting probabilities required for case data")}

  #Process input data to check that all regions with sero and/or case data supplied are present, remove
  #regions without any supplied data, and add cross-referencing tables for use when calculating likelihood. Take
  #subset of environmental data and check that environmental data available for all regions
  input_data = input_data_process(input_data, obs_sero_data, obs_case_data)
  regions = names(table(input_data$region_labels)) #Regions in new processed input data list
  n_regions = length(regions)
  assert_that(all(regions %in% enviro_data_const$region),msg="Environmental data must be available for all regions in observed data")
  enviro_data_const = subset(enviro_data_const, enviro_data_const$region %in% regions)
  assert_that(all(regions==enviro_data_var$regions),msg="Environmental data must be available for all regions in observed data")
  #enviro_data_var = subset(enviro_data_var, enviro_data_var$region %in% regions) #TBA

  #Get names of additional parameters to be estimated
  extra_estimated_params = c()
  for(var_name in names(add_values)){
    if(is.na(add_values[[var_name]]) == TRUE){extra_estimated_params = append(extra_estimated_params, var_name)}
  }

  #Label parameters according to order and fitting type
  param_names = create_param_labels_VarFR(enviro_data_const, enviro_data_var, extra_estimated_params)
  names(log_params_ini) = param_names

  #Run checks on inputs
  checks <- mcmc_checks_VarFR(log_params_ini, n_regions, prior_settings, enviro_data_const, enviro_data_var, add_values, extra_estimated_params)
  if(prior_settings$type == "flat"){names(prior_settings$param_min_limits) = names(prior_settings$param_max_limits) = param_names}

  #Designate constant and variable covariates
  const_covars=colnames(enviro_data_const)[c(2:ncol(enviro_data_const))]
  var_covars=enviro_data_var$env_vars
  covar_names=c(const_covars,var_covars)
  n_env_vars=length(covar_names)
  i_FOI_const = c(1:n_env_vars)[covar_names %in% const_covars]
  i_FOI_var = c(1:n_env_vars)[covar_names %in% var_covars]
  i_R0_const = i_FOI_const + n_env_vars
  i_R0_var = i_FOI_var + n_env_vars

  #MCMC setup
  chain = chain_prop = posterior_current = posterior_prop = flag_accept = chain_cov_all = NULL
  n_params = length(log_params_ini)
  burnin = min(2*n_params, Niter)
  fileIndex = 0
  log_params = log_params_ini
  chain_cov = 1
  adapt = 0
  posterior_value_current = -Inf

  #Iterative estimation
  for (iter in 1:Niter){

    #Propose new parameter values
    log_params_prop = param_prop_setup(log_params, chain_cov, adapt)

    #Calculate likelihood using single_posterior_calc function
    posterior_value_prop = single_posterior_calc_VarFR(log_params_prop, input_data, obs_sero_data, obs_case_data,
                                                       mode_start=mode_start,prior_settings=prior_settings,dt=dt,n_reps=n_reps,
                                                       enviro_data_const=enviro_data_const, enviro_data_var=enviro_data_var,
                                                       p_severe_inf = p_severe_inf, p_death_severe_inf=p_death_severe_inf,
                                                       add_values = add_values, extra_estimated_params = extra_estimated_params,
                                                       deterministic = deterministic, mode_time = mode_time,
                                                       mode_parallel = mode_parallel, cluster = cluster,
                                                       i_FOI_const = i_FOI_const, i_FOI_var = i_FOI_var,
                                                       i_R0_const = i_R0_const, i_R0_var = i_R0_var)
    gc() #Clear garbage to prevent memory creep

    if(is.finite(posterior_value_prop) == FALSE) {
      p_accept = -Inf
    } else {
      p_accept = posterior_value_prop - posterior_value_current
      if(is.na(p_accept) ){ p_accept = -Inf}
    }

    ## accept/reject iteration:
    tmp = runif(1)
    if(tmp<min(exp(p_accept), 1)) {
      log_params = log_params_prop
      posterior_value_current = posterior_value_prop
      accept = 1
    } else {accept = 0}

    #save current iteration
    chain = rbind(chain, log_params)
    chain_prop = rbind(chain_prop, log_params_prop)
    posterior_current = rbind(posterior_current, posterior_value_current)
    posterior_prop = rbind(posterior_prop, posterior_value_prop)
    flag_accept = rbind(flag_accept, accept)
    chain_cov_all = rbind(chain_cov_all, max(chain_cov))

    #Set output headings
    if(iter == 1){
      colnames(chain) = colnames(chain_prop) = names(log_params_ini)
      for(i in 1:n_params){colnames(chain_prop)[i] = paste("Test_", colnames(chain_prop)[i], sep = "")}
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
      lines = min(((fileIndex*10000) + 1), iter):iter

      data_out <- cbind(posterior_current, posterior_prop, exp(chain), flag_accept, exp(chain_prop), chain_cov_all)[lines, ]
      if(fileAccess(fn, 2) == 0){write.csv(data_out, fn, row.names = FALSE)}
    }

    #Decide whether next iteration will be adaptive
    if (iter>burnin & runif(1)<0.9){ #adapt
      adapt = 1
      chain_cov = cov(chain[max(nrow(chain)-10000, 1):nrow(chain), ])
    } else {
      adapt = 0
      chain_cov = 1
    }
  }

  return(NULL)
}
#-------------------------------------------------------------------------------
#' @title single_posterior_calc_VarFR
#'
#' @description Function which calculates and outputs posterior likelihood of observing simulated data
#'
#' @details This function calculates the posterior likelihood of observing a set of observations (across multiple
#' regions and data types) for a given proposed parameter set. [TBA]
#'
#' @param log_params_prop Proposed values of parameters to be estimated (natural logarithm of actual parameters)
#' @param input_data List of population and vaccination data for multiple regions (created using data input
#'   creation code and usually loaded from RDS file), with cross-reference tables added using input_data_process
#'   in MCMC
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param ... = Constant parameters/flags/etc. loaded to or determined by mcmc() and mcmc_prelim_fit, including mode_start,
#'  prior_settings, dt, n_reps, enviro_data, p_severe_inf, p_death_severe_inf, add_values list, extra_estimated_params,
#'  deterministic, mode_time, mode_parallel, cluster, i_FOI_const, i_FOI_var, i_R0_const, i_R0_var
#'
#' @export
#'
single_posterior_calc_VarFR <- function(log_params_prop = c(), input_data = list(), obs_sero_data = NULL, obs_case_data = NULL,...){

  consts=list(...)

  #Check values for flat prior
  prior_like = 0
  if(consts$prior_settings$type == "flat"){
    if(any(log_params_prop<consts$prior_settings$param_min_limits) || any(log_params_prop>consts$prior_settings$param_max_limits)){
      prior_like = -Inf}
  }

  #Get additional values, calculate associated normal-distribution prior values if relevant
  if(is.finite(prior_like)){
    vaccine_efficacy = p_rep_severe = p_rep_death = m_FOI_Brazil = 1.0
    for(var_name in names(consts$add_values)){
      if(var_name %in% consts$extra_estimated_params){
        i = match(var_name, names(log_params_prop))
        value = exp(as.numeric(log_params_prop[i]))
        assign(var_name, value)
        if(consts$prior_settings$type == "norm"){
          prior_like = prior_like + log(dtrunc(value, "norm", a = consts$prior_settings$param_min_limits[i],
                                               b = consts$prior_settings$param_max_limits[i],
                                               mean = consts$prior_settings$norm_params_mean[i],
                                               sd = consts$prior_settings$norm_params_sd[i]))
        }
      } else {assign(var_name, consts$add_values[[var_name]])}
    }
  }

  #If prior is finite so far, get normal-distribution prior values for environmental coefficients if relevant
  if(is.finite(prior_like) && consts$prior_settings$type == "norm"){
    # values=c(1:(2*(ncol(consts$enviro_data)-1)))
    # prior_like = prior_like + sum(dnorm(log_params_prop[values], mean = consts$prior_settings$norm_params_mean[values],
    #                                     sd = consts$prior_settings$norm_params_sd[values], log = TRUE))
    for(i in c(1:(2*(ncol(consts$enviro_data)-1)))){
      prior_like = prior_like + log(dtrunc(log_params_prop[i], "norm", a = consts$prior_settings$param_min_limits[i],
                                           b = consts$prior_settings$param_max_limits[i],
                                           mean = consts$prior_settings$norm_params_mean[i],
                                           sd = consts$prior_settings$norm_params_sd[i]))
    }
  }

  #If prior is finite so far, get FOI and R0 values and calculate any associated prior
  if(is.finite(prior_like)){
    regions = input_data$region_labels
    n_regions = length(regions)

    #CALCULATE FOI AND R0 VALUES FROM CONSTANT AND VARIABLE ENVIRONMENTAL DATA
    {
      FOI_values=calc_var_epi(coeffs_const=exp(log_params_prop[consts$i_FOI_const]),coeffs_var=exp(log_params_prop[consts$i_FOI_var]),
                              enviro_data_const=consts$enviro_data_const,enviro_data_var=consts$enviro_data_var)
      R0_values=calc_var_epi(coeffs_const=exp(log_params_prop[consts$i_R0_const]),coeffs_var=exp(log_params_prop[consts$i_R0_var]),
                             enviro_data_const=consts$enviro_data_const,enviro_data_var=consts$enviro_data_var)
    }

    for(n_region in 1:n_regions){if(substr(regions[n_region], 1, 3) == "BRA"){FOI_values[n_region] = FOI_values[n_region]*m_FOI_Brazil}}
    # if(consts$prior_settings$type == "norm"){ #TBC - Apply prior to mean FOI/R0 over time for each region?
    #   prior_like = prior_like  +
    #     sum(log(dtrunc(R0_values, "norm", a = 0, b = Inf, mean = consts$prior_settings$R0_mean, sd = consts$prior_settings$R0_sd)))  +
    #     sum(log(dtrunc(FOI_values, "norm", a = 0, b = 1, mean = consts$prior_settings$FOI_mean, sd = consts$prior_settings$FOI_sd)))
    # }
  }

  ### If prior finite, evaluate likelihood ###
  if (is.finite(prior_like)) {

    #Generate modelled data over all regions
    dataset <- Generate_Dataset_VarFR(input_data, FOI_values, R0_values, obs_sero_data, obs_case_data, vaccine_efficacy,
                                      consts$p_severe_inf, consts$p_death_severe_inf, p_rep_severe, p_rep_death,
                                      consts$mode_start, start_SEIRV = NULL, consts$dt, consts$n_reps, consts$deterministic,
                                      consts$mode_time, consts$mode_parallel, consts$cluster)

    #Likelihood of observing serological data
    if(is.null(obs_sero_data) == FALSE){
      sero_like_values = sero_data_compare(dataset$model_sero_values, obs_sero_data)
    } else {sero_like_values = 0}

    #Likelihood of observing annual case/death data
    if(is.null(obs_case_data) == FALSE){
      cases_like_values = case_data_compare(dataset$model_case_values, obs_case_data$cases)
      if(is.null(obs_case_data$deaths) == FALSE){
        deaths_like_values = case_data_compare(dataset$model_death_values, obs_case_data$deaths)
      } else {deaths_like_values = 0}
    } else {cases_like_values = deaths_like_values = 0}

    posterior = prior_like + sum(sero_like_values, na.rm = TRUE) + sum(cases_like_values, na.rm = TRUE) + sum(deaths_like_values, na.rm = TRUE)

  } else {posterior = -Inf}

  return(posterior)
}
#-------------------------------------------------------------------------------
#' @title mcmc_checks_VarFR
#'
#' @description Perform checks on MCMC inputs
#'
#' @details This function, which is called by MCMC(), performs a number of checks on data to be used in fitting to
#' ensure proper functionality. It verifies that the number of parameters being estimated is consistent with other
#' settings and that certain values are not outwith sensible boundaries (e.g. probabilities must be between 0 and 1).
#'
#' @param log_params_ini Initial values of parameters to be estimated (natural logarithm of actual parameters; see
#'  documentation for MCMC() function for more details)
#' @param n_regions Number of regions
#' @param prior_settings List containing settings for priors; see documentation for MCMC() function for more details)
#' @param enviro_data_const Data frame of values of constant environmental covariates (columns) by region (rows)
#' @param enviro_data_var List containing values of time-varying environmental covariates (TBA)
#' @param add_values List of parameters in addition to those governing FOI/R0, either giving a fixed value or giving NA to
#'   indicate that they are part of the parameter set to be estimated \cr
#'  vaccine_efficacy Vaccine efficacy (proportion of reported vaccinations causing immunity) (must be present) \cr
#'  p_rep_severe Probability of observation of severe infection \cr
#'  p_rep_death Probability of observation of death \cr
#'  m_FOI_Brazil Multiplier of spillover FOI for Brazil regions (only relevant if regions in Brazil to be considered) \cr
#' @param extra_estimated_params Vector of names of parameters to be estimated in addition to those governing FOI and R0; see add_values
#'
#' @export
#'
mcmc_checks_VarFR <- function(log_params_ini = c(), n_regions = 1, prior_settings = list(type = "zero"),
                              enviro_data_const = list(), enviro_data_var = list(),
                              add_values = list(vaccine_efficacy = 1.0, p_rep_severe = 1.0, p_rep_death = 1.0, m_FOI_Brazil = 1.0),
                              extra_estimated_params = list()){

  param_names = names(log_params_ini)
  n_params = length(log_params_ini)
  assert_that(is.null(param_names) == FALSE, msg = "Parameters should be named using create_param_labels")
  assert_that(prior_settings$type %in% c("zero", "flat", "norm"), msg = "Prior settings type must be 'zero', 'flat' or 'norm'")
  if(prior_settings$type %in% c("flat", "norm")){
    assert_that(length(prior_settings$param_min_limits) == n_params, msg = "Check prior_settings$param_min_limits")
    assert_that(length(prior_settings$param_max_limits) == n_params, msg = "Check prior_settings$param_max_limits")
  }
  if(prior_settings$type == "norm"){
    assert_that(length(prior_settings$norm_params_mean) == n_params, msg = "Check prior_settings$norm_params_mean")
    assert_that(length(prior_settings$norm_params_sd) == n_params, msg = "Check prior_settings$norm_params_sd")
    assert_that(is.numeric(prior_settings$R0_mean), msg = "Check prior_settings$R0_mean")
    assert_that(is.numeric(prior_settings$R0_sd), msg = "Check prior_settings$R0_sd")
    assert_that(is.numeric(prior_settings$FOI_mean), msg = "Check prior_settings$FOI_mean")
    assert_that(is.numeric(prior_settings$FOI_sd), msg = "Check prior_settings$FOI_sd")
  }

  # Check additional values
  add_value_names = names(add_values)
  assert_that("vaccine_efficacy" %in% add_value_names,
              msg="Reported vaccination effectiveness vaccine_efficacy must be included in add_values")
  assert_that(all(extra_estimated_params %in% add_value_names),
              sg="Additional parameters to be estimated must be included in add_values")
  for(var_name in add_value_names){
    if(var_name %in% extra_estimated_params){
      assert_that(is.na(add_values[[var_name]]),msg="Additional parameters to be estimated must be set to NA in add_values")
    } else {assert_that(add_values[[var_name]] <= 1.0 && add_values[[var_name]] >= 0.0,
                        msg="Fixed additional parameters must be in range 0-1")}
  }

  # Get names of environmental covariates
  assert_that(is.null(enviro_data_const) == FALSE, msg="Constant environmental data required")
  assert_that(is.null(enviro_data_var) == FALSE, msg="Variable environmental data required")
  env_vars = c(names(enviro_data_const[c(2:ncol(enviro_data_const))]),enviro_data_var$env_vars) #TBC
  n_env_vars = length(env_vars)

  # Check that total number of parameters is correct; check parameters named in correct order (TBA); all should be correct if
  # parameter names created using create_param_labels
  if(is.null(extra_estimated_params)){n_extra_params=0}else{n_extra_params=length(extra_estimated_params)}
  assert_that(n_params == (2*n_env_vars) + n_extra_params,
              msg="Length of initial parameter vector must equal twice number of environmental covariates +
              number of additional estimated parameters")
  for(i in 1:n_env_vars){
    assert_that(param_names[i] == paste0("FOI_", env_vars[i]), msg="Initial parameter vector must start with FOI coefficients")
    assert_that(param_names[i + n_env_vars] == paste0("R0_", env_vars[i]),
                msg="R0 coefficients must follow FOI coefficients in initial parameter vector")
  }
  if(length(extra_estimated_params)>0){
    assert_that(all(param_names[(2*n_env_vars)+c(1:n_extra_params)] == extra_estimated_params),
                msg="Initial parameter vector must end with additional estimated parameters")
  }

  return(NULL)
}
#-------------------------------------------------------------------------------
#' @title mcmc_prelim_fit_VarFR
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
#' @param log_params_min Initial lower limits of estimated parameter values (natural logarithm of actual limits)
#' @param log_params_max Initial upper limits of estimated parameter values (natural logarithm of actual limits)
#' @param input_data List of population and vaccination data for multiple regions (created using data input creation
#'   code and usually loaded from RDS file)
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination \cr
#'  If mode_start = 0, only vaccinated individuals \cr
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (uniform by age, R0 based only) \cr
#'  If mode_start = 3, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age) \cr
#' @param prior_settings List containing settings for priors: must contain text named "type":
#'  If type = "zero", prior probability is always zero \cr
#'  If type = "norm", prior probability is given by dnorm calculation on parameter values with settings based on vectors of values
#'   in prior_settings: \cr
#'   norm_params_mean and norm_params_sd (vectors of mean and standard deviation values applied to log FOI/R0
#'   parameters and to actual values of additional parameters) \cr
#'   + FOI_mean + FOI_sd (mean + standard deviation of computed FOI, single values)  \cr
#'   + R0_mean + R0_sd (mean + standard deviation of computed R0, single values) \cr
#' @param dt time increment in days (must be 1 or 5)
#' @param n_reps Number of repetitions
#' @param enviro_data_const Data frame of values of constant environmental covariates (columns) by region (rows)
#' @param enviro_data_var List containing values of time-varying environmental covariates (TBA)
#' @param p_severe_inf Probability of an infection being severe
#' @param p_death_severe_inf Probability of a severe infection resulting in death
#' @param add_values List of parameters in addition to those governing FOI/R0, either giving a fixed value or giving NA to
#'   indicate that they are part of the fitted  parameter set \cr
#'  vaccine_efficacy Vaccine efficacy (proportion of reported vaccinations causing immunity) (must be present) \cr
#'  p_rep_severe Probability of observation of severe infection \cr
#'  p_rep_death Probability of observation of death \cr
#'  m_FOI_Brazil Multiplier of spillover FOI for Brazil regions (only relevant if regions in Brazil to be considered)
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' @param mode_time TBA
#' @param mode_parallel TRUE/FALSE - indicate whether to use parallel processing on supplied cluster for speed
#' @param cluster Cluster of threads to use if mode_parallel = TRUE
#' @param plot_graphs TRUE/FALSE - plot graphs of evolving parameter space
#' '
#' @export
#'
mcmc_prelim_fit_VarFR <- function(n_iterations = 1, n_param_sets = 1, n_bounds = 1, log_params_min = NULL,
                                  log_params_max = NULL, input_data = list(), obs_sero_data = list(), obs_case_data = list(),
                                  mode_start = 0, prior_settings = list(type = "zero"), dt = 1.0, n_reps = 1,
                                  enviro_data_const = list(), enviro_data_var=list(),
                                  p_severe_inf = 0.12, p_death_severe_inf = 0.39,
                                  add_values = list(vaccine_efficacy = 1.0,p_rep_severe = 1.0,p_rep_death = 1.0,m_FOI_Brazil = 1.0),
                                  deterministic = TRUE, mode_time=0, mode_parallel = FALSE, cluster = NULL, plot_graphs = FALSE){

  #TODO - Add assertthat functions
  assert_that(mode_start %in% c(0, 1, 3), msg = "mode_start must have value 0, 1 or 3")
  assert_that(length(log_params_min) == length(log_params_max), msg = "Parameter limit vectors must have same lengths")
  assert_that(prior_settings$type %in% c("zero", "norm"), msg = "Prior type must be 'zero' or 'norm'")

  best_fit_results = list()
  n_params = length(log_params_min)

  #Get additional values
  extra_estimated_params = c()
  add_value_names = names(add_values)
  assert_that("vaccine_efficacy" %in% add_value_names)
  for(var_name in add_value_names){
    if(is.na(add_values[[var_name]]) == TRUE){extra_estimated_params = append(extra_estimated_params, var_name)}
  }
  param_names = create_param_labels_VarFR(enviro_data_const, enviro_data_var, extra_estimated_params)

  #TODO - Additional assert_that checks

  assert_that(length(param_names) == n_params)
  names(log_params_min) = names(log_params_max) = param_names

  #Designate constant and variable covariates
  const_covars=colnames(enviro_data_const)[c(2:ncol(enviro_data_const))]
  var_covars=enviro_data_var$env_vars
  covar_names=c(const_covars,var_covars)
  n_env_vars=length(covar_names)
  i_FOI_const = c(1:n_env_vars)[covar_names %in% const_covars]
  i_FOI_var = c(1:n_env_vars)[covar_names %in% var_covars]
  i_R0_const = i_FOI_const + n_env_vars
  i_R0_var = i_FOI_var + n_env_vars

  if(plot_graphs){
    xlabels = param_names
    for(i in 1:n_params){xlabels[i] = substr(xlabels[i], 1, 15)}
    ylabels = 10^c(-8, -6, -4, -3, -2, -1, 0, 1)
    par(mar = c(6, 2, 1, 1))
    ylim = c(min(log_params_min), max(log_params_max))
  }

  for(iteration in 1:n_iterations){
    cat("\nIteration: ", iteration, "\n", sep = "")
    all_param_sets <- lhs(n = n_param_sets, rect = cbind(log_params_min, log_params_max))
    results = data.frame()

    for(set in 1:n_param_sets){
      cat("\n\tSet: ", set, sep = "")
      log_params_prop = all_param_sets[set, ]

      cat("\n\tParams: ", signif(log_params_prop, 3))

      names(log_params_prop) = param_names
      posterior_value = single_posterior_calc_VarFR(log_params_prop, input_data, obs_sero_data, obs_case_data,
                                                    mode_start=mode_start,prior_settings=prior_settings,dt=dt,n_reps=n_reps,
                                                    enviro_data_const=enviro_data_const, enviro_data_var=enviro_data_var,
                                                    p_severe_inf = p_severe_inf, p_death_severe_inf=p_death_severe_inf,
                                                    add_values = add_values, extra_estimated_params = extra_estimated_params,
                                                    deterministic = deterministic, mode_time = mode_time,
                                                    mode_parallel = mode_parallel, cluster = cluster,
                                                    i_FOI_const = i_FOI_const, i_FOI_var = i_FOI_var,
                                                    i_R0_const = i_R0_const, i_R0_var = i_R0_var)
      gc() #Clear garbage to prevent memory creep
      results <- rbind(results, c(set, exp(log_params_prop), posterior_value))
      if(set == 1){colnames(results) = c("set", param_names, "posterior")}

      cat("\n\tPosterior likelihood = ", posterior_value, sep = "")

    }
    results <- results[order(results$posterior, decreasing = TRUE), ]
    best_fit_results[[iteration]] = results

    log_params_min_new = log_params_max_new = rep(0, n_params)
    for(i in 1:n_params){
      log_params_min_new[i] = min(log(results[c(1:n_bounds), i + 1]))
      log_params_max_new[i] = max(log(results[c(1:n_bounds), i + 1]))
    }
    names(log_params_min_new) = names(log_params_max_new) = param_names

    if(plot_graphs){
      matplot(x = c(1:n_params), y = log(t(results[c(1:n_bounds), c(1:n_params) + 1])), type = "p", pch = 16, col = 1,
              xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = ylim)
      axis(side = 1, at = c(1:n_params), labels = xlabels, las = 2, cex.axis = 0.7)
      axis(side = 2, at = log(ylabels), labels = ylabels)
      matplot(x = c(1:n_params), y = log_params_min, type = "l", col = 1, lty = 2, add = TRUE)
      matplot(x = c(1:n_params), y = log_params_max, type = "l", col = 1, lty = 2, add = TRUE)
      matplot(x = c(1:n_params), y = log_params_min_new, type = "l", col = 2, add = TRUE)
      matplot(x = c(1:n_params), y = log_params_max_new, type = "l", col = 2, add = TRUE)
    }

    log_params_min = log_params_min_new
    log_params_max = log_params_max_new
  }

  return(best_fit_results)
}
#-------------------------------------------------------------------------------
#' @title create_param_labels_VarFR
#'
#' @description Apply names to the parameters in a set used for data matching and parameter fitting
#'
#' @details Takes in environmental covariate data along with names of additional parameters (vaccine efficacy
#' and reporting probabilities) and generates list of names for parameter set to use as input for fitting functions
#'
#' @param enviro_data_const TBA
#' @param enviro_data_var TBA
#' @param extra_estimated_params Vector of strings listing variable parameters besides ones determining FOI/R0 (may include
#' vaccine efficacy and/or infection/death reporting probabilities and/or Brazil FOI adjustment factor)
#'
#' @export
#'
create_param_labels_VarFR <- function(enviro_data_const = NULL, enviro_data_var=list(), extra_estimated_params = c("vacc_eff")){

  #TODO - Assert_that functions

  if(is.null(extra_estimated_params)){n_extra=0}else{n_extra = length(extra_estimated_params)}
  env_vars = c(names(enviro_data_const[c(2:ncol(enviro_data_const))]),enviro_data_var$env_vars) #TBC
  n_env_vars = length(env_vars)
  n_params = (2*n_env_vars)+n_extra
  param_names = rep("", n_params)
  for(i in 1:n_env_vars){
    param_names[i] = paste("FOI_", env_vars[i], sep = "")
    param_names[i+n_env_vars] = paste("R0_", env_vars[i], sep = "")
  }
  if(n_extra>0){param_names[(n_params-n_extra+1):n_params] = extra_estimated_params}

  return(param_names)
}
#-------------------------------------------------------------------------------
#' @title calc_var_epi
#'
#' @description Calculate time-varying FOI_spillover or R0 values from environmental covariates and coefficients
#'
#' @details TBA
#'
#' @param coeffs_const TBA
#' @param coeffs_var TBA
#' @param enviro_data_const TBA
#' @param enviro_data_var TBA
#' '
#' @export
#'
calc_var_epi <- function(coeffs_const = c(), coeffs_var = c(), enviro_data_const = data.frame(), enviro_data_var = list()){
  n_pts=dim(enviro_data_var$values)[3]
  #TODO - Add assertthat checks
  #TODO - Ensure function works if only variable or only constant covariates

  base_output_values=as.numeric(colSums(coeffs_const*t(enviro_data_const[,c(2:ncol(enviro_data_const))])))
  var_output_values=colSums(coeffs_var*enviro_data_var$values)

  total_output_values=array(NA,dim=dim(var_output_values))
  for(i in 1:n_pts){
    total_output_values[,i]=var_output_values[,i]+base_output_values
  }

  return(total_output_values)
}
