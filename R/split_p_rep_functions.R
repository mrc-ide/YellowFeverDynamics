#-------------------------------------------------------------------------------
#TODO: Harmonize order of variables between model running, dataset generation, and MCMC functions (see YEP)
#' @title MCMC2
#'
#' @description Combined MCMC Multi-Region - series of MCMC iterations for one or more regions (split reporting probabilities)
#'
#' @details This is a version of the master function for running a Markov chain to optimize the parameters of the yellow fever
#' model based on the calculated likelihood of observing supplied data given a particular set of parameters, using different
#' values of the infection reporting probabilities for African and South American countries.
#'
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
#'   probability p_rep_severe_af/p_rep_severe_sa, and fatal case reporting probability p_rep_death_af/p_rep_death_sa, and
#'   Brazil spillover FOI multiplier m_FOI_Brazil); if these are to be estimated, in the order
#'   vaccine_efficacy->p_rep_severe_af->p_rep_severe_sa->p_rep_death_af->p_rep_death_sa->m_FOI_Brazil.
#'   If these parameters are to be estimated, the values separately supplied to this function (see add_values below) should
#'   be set to NA.
#' @param input_data List of population and vaccination data for multiple regions (created using data input creation
#'   code and usually loaded from RDS file)
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no.
#'   cases/no. deaths
#' @param filename_prefix Prefix of names for output files; function outputs a CSV file every 10, 000 iterations with a
#'   name in the format: "(filename_prefix)XX.csv", e.g. Chain00.csv
#' @param Niter Total number of iterations to run
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination \cr
#'   If mode_start = 0, only vaccinated individuals \cr
#'   If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (uniform by age,
#'   R0 based only) \cr
#'   If mode_start = 3, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age)
#' @param prior_settings List containing settings for priors: must contain text named "type":
#'   If type = "zero", prior probability is always zero \cr
#'   If type = "flat", prior probability is zero if log parameter values in designated ranges param_min_limits and
#'   param_max_limits, -Inf otherwise; param_min_limits and param_max_limits included in prior_settings as vectors of
#'   same length as log_params_ini \cr
#'   If type = "norm", prior probability is given by truncated normal distribution with settings based on vectors of values
#'   in prior_settings: \cr
#'   norm_params_mean and norm_params_sd (vectors of mean and standard deviation values applied to log FOI/R0
#'   parameters and to actual values of additional parameters) \cr
#'   + FOI_mean + FOI_sd (mean + standard deviation of computed FOI, single values)  \cr
#'   + R0_mean + R0_sd (mean + standard deviation of computed R0, single values) \cr
#'   + param_min_limits and param_max_limits (lower and upper limits applied to truncated normal distributions)
#' @param time_inc time increment in days (must be 1 or 5)
#' @param n_reps Number of times to repeat calculations to get average likelihood at each iteration
#' @param enviro_data_const Data frame of values of constant environmental covariates (columns) by region (rows)
#' @param enviro_data_var List containing time-varying environmental covariate data:\cr
#'   regions: Vector of region labels\cr
#'   env_vars: Vector of covariate names\cr
#'   values: Array of covariate values with dimensions (number of covariates, number of regions, number of time points).
#'   Number of time points must be correct for mode_time setting.\cr
#' @param p_severe_inf Probability of an infection being severe
#' @param p_death_severe_inf Probability of a severe infection resulting in death
#' @param add_values List of parameters in addition to those governing FOI/R0, either giving a fixed value or giving NA to
#'   indicate that they are part of the fitted  parameter set \cr
#'   vaccine_efficacy Vaccine efficacy (proportion of reported vaccinations causing immunity) (must be present) \cr
#'   p_rep_severe_af Probability of observation of severe infection in African regions
#'   p_rep_severe_sa Probability of observation of severe infection in South American regions
#'   p_rep_death_af Probability of observation of death in African regions
#'   p_rep_death_sa Probability of observation of death in South American regions
#'   m_FOI_Brazil Multiplier of spillover FOI for Brazil regions (only relevant if regions in Brazil to be considered)
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' @param mode_time Type of time dependence of FOI_spillover and R0 to be used: \cr
#'   If mode_time = 0, no time variation (constant values)\cr
#'   If mode_time = 1, FOI/R0 vary annually without seasonality (number of values = number of years to consider) \cr
#'   If mode_time = 2, FOI/R0 vary with monthly seasonality without inter-annual variation (number of values = 12) \cr
#'   If mode_time = 3, FOI/R0 vary with daily seasonality without inter-annual variation (number of values = 365/dt) \cr
#'   If mode_time = 4, FOI/R0 vary annually with monthly seasonality (number of values = 12*number of years to consider) \cr
#'   If mode_time = 5, FOI/R0 vary annually with daily seasonality (number of values = (365/dt)*number of years to consider)
#' @param mode_parallel TRUE/FALSE - indicate whether to use parallel processing on supplied cluster for speed
#' @param cluster Cluster of threads to use if mode_parallel = TRUE
#' '
#' @export
#'
MCMC2 <- function(log_params_ini = c(), input_data = list(), obs_sero_data = NULL, obs_case_data = NULL,
                  filename_prefix = "Chain", Niter = 1, mode_start = 0, prior_settings = list(type = "zero"),
                  time_inc = 1.0, n_reps = 1, enviro_data_const = list(), enviro_data_var = list(),
                  p_severe_inf = 0.12, p_death_severe_inf = 0.39,
                  add_values = list(vaccine_efficacy = 1.0, p_rep_severe_af = 1.0, p_rep_severe_sa = 1.0,
                                    p_rep_death_af = 1.0, p_rep_death_sa = 1.0, m_FOI_Brazil = 1.0),
                  deterministic = FALSE, mode_time = 1, mode_parallel = FALSE, cluster = NULL){

  assert_that(is.logical(deterministic))
  assert_that(mode_start %in% c(0, 1, 3), msg = "mode_start must have value 0, 1 or 3")
  assert_that(all(c("p_rep_severe_af", "p_rep_severe_sa", "p_rep_death_af", "p_rep_death_sa") %in% names(add_values)),
              msg = "Reporting probabilities required for case data")

  #Process input data to check that all regions with sero and/or case data supplied are present, remove
  #regions without any supplied data, and add cross-referencing tables for use when calculating likelihood. Take
  #subset of environmental data and check that environmental data available for all regions
  input_data = input_data_process(input_data, obs_sero_data, obs_case_data)
  regions = names(table(input_data$region_labels)) #Regions in new processed input data list
  n_regions = length(regions)
  assert_that(all(regions %in% enviro_data_const$region),
              msg = "Time-invariant environmental data must be available for all regions in observed data")
  enviro_data_const = subset(enviro_data_const, enviro_data_const$region %in% regions)
  assert_that(enviro_data_var_check(enviro_data_var))
  assert_that(all(regions %in% enviro_data_var$regions),
              msg = "Time-variant environmental data must be available for all regions in observed data")
  enviro_data_var = enviro_data_var_truncate(enviro_data_var, regions)

  #Get names of additional parameters to be estimated
  extra_estimated_params = c()
  for(var_name in names(add_values)){
    if(is.na(add_values[[var_name]]) == TRUE){extra_estimated_params = append(extra_estimated_params, var_name)}
  }

  #Label parameters according to order and fitting type
  param_names = create_param_labels(enviro_data_const, enviro_data_var, extra_estimated_params)
  names(log_params_ini) = param_names

  #Run checks on inputs
  checks <- mcmc_checks(log_params_ini, n_regions, prior_settings, enviro_data_const, enviro_data_var, add_values,
                        extra_estimated_params)
  if(prior_settings$type == "flat"){
    names(prior_settings$param_min_limits) = names(prior_settings$param_max_limits) = param_names
    }

  #Designate constant and variable covariates
  const_covars = colnames(enviro_data_const)[c(2:ncol(enviro_data_const))]
  var_covars = enviro_data_var$env_vars
  covar_names = c(const_covars, var_covars)
  n_env_vars = length(covar_names)
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

    #Calculate likelihood using single_posterior_calc2 function
    posterior_value_prop = single_posterior_calc2(log_params_prop, input_data, obs_sero_data, obs_case_data,
                                                 mode_start = mode_start, prior_settings = prior_settings, time_inc = time_inc,
                                                 n_reps = n_reps, enviro_data_const = enviro_data_const,
                                                 enviro_data_var = enviro_data_var, p_severe_inf = p_severe_inf,
                                                 p_death_severe_inf = p_death_severe_inf, add_values = add_values,
                                                 extra_estimated_params = extra_estimated_params,
                                                 deterministic = deterministic, mode_time = mode_time,
                                                 mode_parallel = mode_parallel, cluster = cluster,
                                                 i_FOI_const = i_FOI_const, i_FOI_var = i_FOI_var,
                                                 i_R0_const = i_R0_const, i_R0_var = i_R0_var, n_env_vars = n_env_vars)
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

    #Output chain to file every 10 iterations; start new file every 10, 000 iterations
    if (iter %% 10 == 0){
      if (iter %% 10000 == 10){fileIndex = (iter-10)/10000}
      if(fileIndex >=  10){fn = paste0(filename_prefix, fileIndex, ".csv")}else{fn = paste0(filename_prefix, "0", fileIndex, ".csv")}
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
#' @title single_posterior_calc2
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
#' @param ... = Constant parameters/flags/etc. loaded to or determined by MCMC2() and mcmc_prelim_fit2, including mode_start,
#' prior_settings, time_inc, n_reps, enviro_data, p_severe_inf, p_death_severe_inf, add_values, extra_estimated_params,
#' deterministic, mode_time, mode_parallel, cluster, i_FOI_const, i_FOI_var, i_R0_const, i_R0_var
#'
#' @export
#'
single_posterior_calc2 <- function(log_params_prop = c(), input_data = list(), obs_sero_data = NULL, obs_case_data = NULL, ...){

  consts = list(...)

  #Check values for flat prior
  prior_like = 0
  if(consts$prior_settings$type == "flat"){
    if(any(log_params_prop<consts$prior_settings$param_min_limits) ||
       any(log_params_prop>consts$prior_settings$param_max_limits)){
      prior_like = -Inf}
  }

  #Get additional values, calculate associated normal-distribution prior values if relevant
  if(is.finite(prior_like)){
    vaccine_efficacy = p_rep_severe_af = p_rep_severe_sa = p_rep_death_af = p_rep_death_sa = m_FOI_Brazil = 1.0
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
    for(i in 1:(2*consts$n_env_vars)){
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

    FOI_values = epi_param_calc(coeffs_const = exp(log_params_prop[consts$i_FOI_const]),
                              coeffs_var = exp(log_params_prop[consts$i_FOI_var]),
                              enviro_data_const = consts$enviro_data_const, enviro_data_var = consts$enviro_data_var)
    R0_values = epi_param_calc(coeffs_const = exp(log_params_prop[consts$i_R0_const]),
                             coeffs_var = exp(log_params_prop[consts$i_R0_var]),
                             enviro_data_const = consts$enviro_data_const, enviro_data_var = consts$enviro_data_var)

    for(n_region in 1:n_regions){
      if(substr(regions[n_region], 1, 3) == "BRA"){FOI_values[n_region] = FOI_values[n_region]*m_FOI_Brazil}
      }
    if(consts$prior_settings$type == "norm"){
      prior_like = prior_like  +
        sum(log(dtrunc(rowMeans(R0_values),
                       "norm", a = 0, b = Inf, mean = consts$prior_settings$R0_mean, sd = consts$prior_settings$R0_sd))) +
        sum(log(dtrunc(rowMeans(FOI_values),
                       "norm", a = 0, b = 1, mean = consts$prior_settings$FOI_mean, sd = consts$prior_settings$FOI_sd)))
    }
  }

  ### If prior finite, evaluate likelihood ###
  if (is.finite(prior_like)) {

    #Generate modelled data over all regions
    dataset <- Generate_Dataset2(FOI_values, R0_values, input_data, obs_sero_data, obs_case_data, vaccine_efficacy,
                                consts$time_inc, consts$mode_start, start_SEIRV = NULL, consts$mode_time, consts$n_reps,
                                consts$deterministic, consts$p_severe_inf, consts$p_death_severe_inf,
                                p_rep_severe_af, p_rep_severe_sa, p_rep_death_af, p_rep_death_sa,
                                consts$mode_parallel, consts$cluster, output_frame = FALSE)

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

    posterior = prior_like + sum(sero_like_values, na.rm = TRUE) + sum(cases_like_values, na.rm = TRUE) +
      sum(deaths_like_values, na.rm = TRUE)

  } else {posterior = -Inf}

  return(posterior)
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
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (uniform by age,
#'  R0 based only) \cr
#'  If mode_start = 3, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age)
#' @param prior_settings List containing settings for priors: must contain text named "type":
#'  If type = "zero", prior probability is always zero \cr
#'  If type = "norm", prior probability is given by truncated normal distribution with settings based on vectors
#'  of values in prior_settings: \cr
#'   norm_params_mean and norm_params_sd (vectors of mean and standard deviation values applied to log FOI/R0
#'   parameters and to actual values of additional parameters) \cr
#'   + FOI_mean + FOI_sd (mean + standard deviation of computed FOI, single values)  \cr
#'   + R0_mean + R0_sd (mean + standard deviation of computed R0, single values) \cr
#' @param time_inc time increment in days (must be 1 or 5)
#' @param n_reps Number of repetitions
#' @param enviro_data_const Data frame of values of constant environmental covariates (columns) by region (rows)
#' @param enviro_data_var List containing values of time-varying environmental covariates (TBA)
#' @param p_severe_inf Probability of an infection being severe
#' @param p_death_severe_inf Probability of a severe infection resulting in death
#' @param add_values List of parameters in addition to those governing FOI/R0, either giving a fixed value or giving NA to
#'   indicate that they are part of the fitted  parameter set \cr
#'   vaccine_efficacy Vaccine efficacy (proportion of reported vaccinations causing immunity) (must be present) \cr
#'   p_rep_severe_af Probability of observation of severe infection in African regions
#'   p_rep_severe_sa Probability of observation of severe infection in South American regions
#'   p_rep_death_af Probability of observation of death in African regions
#'   p_rep_death_sa Probability of observation of death in South American regions
#'   m_FOI_Brazil Multiplier of spillover FOI for Brazil regions (only relevant if regions in Brazil to be considered)
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' @param mode_time TBA
#' @param mode_parallel TRUE/FALSE - indicate whether to use parallel processing on supplied cluster for speed
#' @param cluster Cluster of threads to use if mode_parallel = TRUE
#' @param plot_graphs TRUE/FALSE - plot graphs of evolving parameter space
#' '
#' @export
#'
mcmc_prelim_fit2 <- function(n_iterations = 1, n_param_sets = 1, n_bounds = 1, log_params_min = NULL,
                            log_params_max = NULL, input_data = list(), obs_sero_data = list(), obs_case_data = list(),
                            mode_start = 0, prior_settings = list(type = "zero"), time_inc = 1.0, n_reps = 1,
                            enviro_data_const = list(), enviro_data_var = list(), p_severe_inf = 0.12, p_death_severe_inf = 0.39,
                            add_values = list(vaccine_efficacy = 1.0, p_rep_severe_af = 1.0, p_rep_severe_sa = 1.0,
                                              p_rep_death_af = 1.0, p_rep_death_sa = 1.0, m_FOI_Brazil = 1.0),
                            deterministic = TRUE, mode_time = 0, mode_parallel = FALSE, cluster = NULL, plot_graphs = FALSE){

  #TODO - Add assertthat functions
  assert_that(is.logical(deterministic))
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
  param_names = create_param_labels(enviro_data_const, enviro_data_var, extra_estimated_params)

  #TODO - Additional assert_that checks

  assert_that(length(param_names) == n_params)
  names(log_params_min) = names(log_params_max) = param_names

  #Designate constant and variable covariates
  const_covars = colnames(enviro_data_const)[c(2:ncol(enviro_data_const))]
  var_covars = enviro_data_var$env_vars
  covar_names = c(const_covars, var_covars)
  n_env_vars = length(covar_names)
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
      posterior_value = single_posterior_calc2(log_params_prop, input_data, obs_sero_data, obs_case_data,
                                              mode_start = mode_start, prior_settings = prior_settings, time_inc = time_inc,
                                              n_reps = n_reps,
                                              enviro_data_const = enviro_data_const, enviro_data_var = enviro_data_var,
                                              p_severe_inf = p_severe_inf, p_death_severe_inf = p_death_severe_inf,
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
#' @title create_param_labels
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
create_param_labels <- function(enviro_data_const = NULL, enviro_data_var = list(),
                                extra_estimated_params = c("vacc_eff")){

  #TODO - Assert_that functions

  if(is.null(extra_estimated_params)){n_extra = 0}else{n_extra = length(extra_estimated_params)}
  covar_names = c(names(enviro_data_const[c(2:ncol(enviro_data_const))]), enviro_data_var$env_vars) #TBC
  n_env_vars = length(covar_names)
  n_params = (2*n_env_vars)+n_extra
  param_names = rep("", n_params)
  for(i in 1:n_env_vars){
    param_names[i] = paste("FOI_", covar_names[i], sep = "")
    param_names[i+n_env_vars] = paste("R0_", covar_names[i], sep = "")
  }
  if(n_extra>0){param_names[(n_params-n_extra+1):n_params] = extra_estimated_params}

  return(param_names)
}
#-------------------------------------------------------------------------------
#' @title Generate_Dataset2
#'
#' @description Generate annual serological and/or case/death data with split reporting probabilities
#'
#' @details This function is used to generate annual serological and/or case/death data based on templates;
#' it is normally used by the single_posterior_calc() function. This is a version which uses different
#' infection reporting probabilities
#'
#' [TBA - Explanation of breakdown of regions to model and how to set lengths of FOI_values and R0_values]
#'
#' #TODO: Harmonize order of variables between model running, dataset generation, and MCMC functions (see YEP)
#'
#' @param FOI_values Array of values of force of infection due to spillover from sylvatic reservoir by region + time point
#' @param R0_values Array of values of basic reproduction number for human-human transmission by region and time point
#' @param input_data List of population and vaccination data for multiple regions in standard format [TBA]
#' @param sero_template Seroprevalence data template - data frame with region, year, minimum/maximum age, vc_factor [TBA]
#'   and number of samples
#' @param case_template Annual reported case/death data template - data frame with region and year
#' @param vaccine_efficacy Fractional vaccine efficacy
#' @param time_inc Time increment in days to use in model (should be either 1.0, 2.5 or 5.0 days)
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination \cr
#'  If mode_start = 0, only vaccinated individuals \cr
#'  If mode_start = 1, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age) \cr
#'  If mode_start = 2, use SEIRV input in list from previous run(s)
#' @param start_SEIRV SEIRV data from end of a previous run to use as input (list of datasets, one per region)
#' @param mode_time Type of time dependence of FOI_spillover and R0 to be used: \cr
#'  If mode_time = 0, no time variation (constant values)\cr
#'  If mode_time = 1, FOI/R0 vary annually without seasonality (number of values = number of years to consider) \cr
#'  If mode_time = 2, FOI/R0 vary with monthly seasonality without inter-annual variation (number of values = 12) \cr
#'  If mode_time = 3, FOI/R0 vary with daily seasonality without inter-annual variation (number of values = 365/dt) \cr
#'  If mode_time = 4, FOI/R0 vary annually with monthly seasonality (number of values = 12*number of years to consider) \cr
#'  If mode_time = 5, FOI/R0 vary annually with daily seasonality (number of values = (365/dt)*number of years to consider)
#' @param n_reps number of stochastic repetitions
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' @param p_severe_inf Probability of an infection being severe
#' @param p_death_severe_inf Probability of a severe infection resulting in death
#' @param p_rep_severe_af Probability of reporting of a severe but non-fatal infection (Africa)
#' @param p_rep_severe_sa Probability of reporting of a severe but non-fatal infection (South America)
#' @param p_rep_death_af Probability of reporting of a fatal infection (Africa)
#' @param p_rep_death_sa Probability of reporting of a fatal infection (South America)
#' @param mode_parallel TRUE/FALSE - set model to run in parallel using cluster if TRUE
#' @param cluster Cluster of threads to use if mode_parallel = TRUE
#' @param output_frame TRUE/FALSE - indicate whether to output a complete data frame of results in template format (if TRUE)
#'   or calculated values only (if FALSE)
#' '
#' @export
#'
Generate_Dataset2 <- function(FOI_values = c(), R0_values = c(), input_data = list(), sero_template = NULL, case_template = NULL,
                              vaccine_efficacy = 1.0, time_inc = 1.0, mode_start = 1, start_SEIRV = NULL, mode_time = 0,
                              n_reps = 1, deterministic = FALSE, p_severe_inf = 0.12, p_death_severe_inf = 0.39,
                              p_rep_severe_af = 1.0, p_rep_severe_sa = 1.0, p_rep_death_af = 1.0, p_rep_death_sa = 1.0,
                              mode_parallel = FALSE, cluster = NULL, output_frame = FALSE){

  assert_that(input_data_check(input_data), msg = paste("Input data must be in standard format",
                                                       " (see https://mrc-ide.github.io/YEP/articles/CGuideAInputs.html)"))
  assert_that(any(is.null(sero_template) == FALSE, is.null(case_template) == FALSE), msg = "Need at least one template")
  if(is.null(sero_template) == FALSE){
    assert_that(all(c("region", "year", "age_min", "age_max", "samples", "vc_factor") %in% names(sero_template)))
  }
  if(is.null(case_template) == FALSE){
    assert_that(all(c("region", "year") %in% names(case_template)))
    assert_that(p_severe_inf >= 0.0 && p_severe_inf <= 1.0, msg = "Severe infection rate must be between 0-1")
    assert_that(p_death_severe_inf >= 0.0 && p_death_severe_inf <= 1.0,
                msg = "Fatality rate of severe infections must be between 0-1")
    assert_that(p_rep_severe_af >= 0.0 && p_rep_severe_af <= 1.0, msg = "Severe infection reporting probability must be between 0-1")
    assert_that(p_rep_severe_sa >= 0.0 && p_rep_severe_sa <= 1.0, msg = "Severe infection reporting probability must be between 0-1")
    assert_that(p_rep_death_af >= 0.0 && p_rep_death_af <= 1.0, msg = "Fatal infection reporting probability must be between 0-1")
    assert_that(p_rep_death_sa >= 0.0 && p_rep_death_sa <= 1.0, msg = "Fatal infection reporting probability must be between 0-1")
  }
  assert_that(is.logical(mode_parallel))
  if(mode_parallel){assert_that(is.null(cluster) == FALSE)}

  #Prune input data based on regions
  regions = regions_breakdown(c(sero_template$region, case_template$region))
  input_data = input_data_truncate(input_data, regions)
  n_regions = length(input_data$region_labels)

  #Cross-reference templates with input regions
  if(is.null(sero_template) == FALSE){
    xref_sero = template_region_xref(sero_template, input_data$region_labels)
    sero_line_list = xref_sero$line_list
  } else {
    xref_sero = data.frame(year_data_begin = rep(Inf, n_regions), year_end = rep(-Inf, n_regions))
    sero_line_list = rep(NA, n_regions)
  }
  if(is.null(case_template) == FALSE){
    xref_case = template_region_xref(case_template, input_data$region_labels)
    case_line_list = xref_case$line_list
  } else {
    xref_case = data.frame(year_data_begin = rep(Inf, n_regions), year_end = rep(-Inf, n_regions))
    case_line_list = rep(NA, n_regions)
  }
  year_data_begin = year_end = rep(NA, length(input_data$region_labels))
  for(i in 1:length(year_data_begin)){
    year_data_begin[i] = min(xref_sero$year_data_begin[i], xref_case$year_data_begin[i])
    year_end[i] = max(xref_sero$year_end[i], xref_case$year_end[i])
  }

  assert_that(length(dim(FOI_values)) == 2, msg = "FOI_values must be 2-D array")
  assert_that(length(dim(R0_values)) == 2, msg = "R0_values must be 2-D array")
  assert_that(dim(FOI_values)[1] == n_regions, msg = "1st dimension of FOI_values must match number of regions to be modelled")
  assert_that(dim(R0_values)[1] == n_regions, msg = "1st dimension of R0_values must match number of regions to be modelled")
  if(mode_start == 2){assert_that(length(start_SEIRV) == n_regions,
                                  msg = "Number of start_SEIRV datasets must match number of regions")}

  #Set up data structures to take modelled data corresponding to observed data
  if(is.null(sero_template)){model_sero_data = NULL} else {
    blank = rep(0, nrow(sero_template))
    model_sero_data = data.frame(samples = blank, positives = blank, sero = blank)
  }
  if(is.null(case_template)){
    model_case_values = model_death_values = NA
  } else {
    model_case_values = model_death_values = rep(0, nrow(case_template))
    country_list_af=c("AGO", "BDI", "BEN", "BFA", "CAF", "CIV", "CMR", "COD", "COG", "ERI", "ETH", "GAB", "GHA", "GIN", "GMB", "GNB", "GNQ",
                      "KEN", "LBR", "MLI", "MRT", "NER", "NGA", "RWA", "SDN", "SEN", "SLE", "SOM", "SSD", "TCD", "TGO", "TZA", "UGA", "ZMB")
    country_list_sa=c("ARG", "BOL", "BRA", "COL", "ECU", "GUF", "GUY", "PER", "PRY", "SUR", "VEN")
    countries=substr(regions, 1, 3)
    assert_that(all(countries %in% c(country_list_af, country_list_sa)))
  }

  #Set up vector of output types to get from model
  output_types = rep(NA, n_regions)
  for(n_region in 1:n_regions){
    if(is.na(case_line_list[[n_region]][1]) == FALSE){
      if(is.na(sero_line_list[[n_region]][1]) == FALSE){ output_types[n_region] = "case+sero"}else{output_types[n_region] = "case"}
    } else {
      output_types[n_region] = "sero"
    }
  }

  if(mode_parallel){
    vacc_data_subsets = pop_data_subsets = years_data_sets = list() #TODO - change input data?
    for(n_region in 1:n_regions){
      vacc_data_subsets[[n_region]] = input_data$vacc_data[n_region, , ]
      pop_data_subsets[[n_region]] = input_data$pop_data[n_region, , ]
      years_data_sets[[n_region]] = c(year_data_begin[n_region]:year_end[n_region])
    }
    if(is.null(start_SEIRV)){start_SEIRV = rep(NA, n_regions)}
    model_output_all = clusterMap(cl = cluster, fun = Model_Run, FOI_spillover = FOI_values, R0 = R0_values,
                                  vacc_data = vacc_data_subsets, pop_data = pop_data_subsets,
                                  years_data = years_data_sets, start_SEIRV = start_SEIRV, output_type = output_types,
                                  MoreArgs = list(year0 = input_data$years_labels[1], vaccine_efficacy = vaccine_efficacy,
                                                  time_inc = time_inc, mode_start = mode_start, mode_time = mode_time,
                                                  n_particles = n_reps, n_threads = 1 , deterministic = deterministic))
  }

  #Save relevant output data from each region
  for(n_region in 1:n_regions){

    #Run model if not using parallelization
    if(mode_parallel == FALSE){
      model_output = Model_Run(FOI_spillover = FOI_values[n_region, ], R0 = R0_values[n_region, ],
                               vacc_data = input_data$vacc_data[n_region, , ], pop_data = input_data$pop_data[n_region, , ],
                               years_data = c(year_data_begin[n_region]:year_end[n_region]), year0 = input_data$years_labels[1],
                               vaccine_efficacy = vaccine_efficacy, time_inc = time_inc, output_type = output_types[n_region],
                               mode_start = mode_start, start_SEIRV = start_SEIRV[[n_region]], mode_time = mode_time,
                               n_particles = n_reps, n_threads = n_reps, deterministic = deterministic)
    } else {
      model_output = model_output_all[[n_region]]
    }
    t_pts = length(model_output$year)

    #Compile case data if needed
    if(is.na(case_line_list[[n_region]][1]) == FALSE){
      #Get reporting probabilities based on region
      if(countries[n_region] %in% country_list_af){
        p_rep_severe=p_rep_severe_af
        p_rep_death=p_rep_death_af
      } else {
        p_rep_severe=p_rep_severe_sa
        p_rep_death=p_rep_death_sa
      }
      case_line_list_region = case_line_list[[n_region]]
      years_case = case_template$year[case_line_list_region]
      n_lines = length(case_line_list_region)

      for(n_rep in 1:n_reps){
        rep_cases = rep_deaths = rep(0, n_lines)
        for(n_line in 1:n_lines){
          pts = c(1:t_pts)[model_output$year == years_case[n_line]]
          infs = sum(model_output$C[n_rep, pts])
          if(deterministic){
            severe_infs = floor(infs)*p_severe_inf
            deaths = severe_infs*p_death_severe_inf
            rep_deaths[n_line] = round(deaths*p_rep_death)
            rep_cases[n_line] = rep_deaths[n_line]+round((severe_infs-deaths)*p_rep_severe)
          } else {
            severe_infs = rbinom(1, floor(infs), p_severe_inf)
            deaths = rbinom(1, severe_infs, p_death_severe_inf)
            rep_deaths[n_line] = rbinom(1, deaths, p_rep_death)
            rep_cases[n_line] = rep_deaths[n_line]+rbinom(1, floor(severe_infs-deaths), p_rep_severe)
          }
        }
        model_case_values[case_line_list_region] = model_case_values[case_line_list_region]+rep_cases
        model_death_values[case_line_list_region] = model_death_values[case_line_list_region]+rep_deaths
      }
    }

    #Compile seroprevalence data if necessary
    if(is.na(sero_line_list[[n_region]][1]) == FALSE){
      sero_line_list_region = sero_line_list[[n_region]]
      for(n_rep in 1:n_reps){
        sero_results = sero_calculate2(sero_template[sero_line_list_region, ], model_output, n_rep)
        model_sero_data$samples[sero_line_list_region] = model_sero_data$samples[sero_line_list_region]+sero_results$samples
        model_sero_data$positives[sero_line_list_region] = model_sero_data$positives[sero_line_list_region] +
          sero_results$positives
      }
    }
  }

  if(is.null(sero_template) == FALSE){model_sero_data$sero = model_sero_data$positives/model_sero_data$samples}
  if(is.null(case_template) == FALSE){
    model_case_values = model_case_values/n_reps
    model_death_values = model_death_values/n_reps
  }

  if(output_frame) { #Output complete frames of data
    return(list(model_sero_data = data.frame(region = sero_template$region, year = sero_template$year,
                                             age_min = sero_template$age_min, age_max = sero_template$age_max,
                                             samples = sero_template$samples,
                                             positives = sero_template$samples*model_sero_data$sero),
                model_case_data = data.frame(region = case_template$region, year = case_template$year,
                                             cases = model_case_values, deaths = model_death_values)))
  } else { #Minimal output for MCMC
    return(list(model_sero_values = model_sero_data$sero, model_case_values = model_case_values,
                model_death_values = model_death_values))
  }
}
