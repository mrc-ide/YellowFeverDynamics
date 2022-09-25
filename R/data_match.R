#Functions for generating sets of modelled data to compare with observed data and displaying comparative graphs
#-------------------------------------------------------------------------------
#' @title data_match_single
#'
#' @description Function which runs the model to create simulated data corresponding to supplied observed data for a
#'   single set of parameters with one or more repetitions
#'
#' @details This function takes in observed data and produces a corresponding simulated dataset for comparison using
#'   supplied model parameters (in the param_prop input variable), population and vaccination settings (in the
#'   input_data input variable) and other settings including parameter type, time increment and environmental covariate
#'   values (in the const_list input variable). The simulated dataset will be able to be compared directly to the
#'   observed dataset - individual data values will be produced for the same years and regions (or combinations of
#'   regions), seroprevalence data will be produced for the same age groups, etc.
#'
#' @param param_prop Log values of proposed parameters
#' @param input_data List of population and vaccination data for multiple regions, with tables to cross-reference
#' with observed data, added using input_data_process
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param obs_outbreak_data Outbreak Y/N data for comparison, by region and year, in format 0 = no outbreaks,
#'   1 = 1 or more outbreak(s)
#' @param const_list = List of constant parameters/flags/etc. (type,n_reps,mode_start,dt,enviro_data,R0_fixed_values,
#'   vaccine_efficacy,p_rep_severe,p_rep_death)
#'
#' @export
#'
data_match_single <- function(param_prop=c(),input_data=list(),obs_sero_data=NULL,obs_case_data=NULL,
                              obs_outbreak_data=NULL,const_list=list()) {

  assert_that(is.numeric(param_prop))
  assert_that(input_data_check(input_data))
  assert_that(any(is.null(obs_sero_data)==FALSE,is.null(obs_case_data)==FALSE,is.null(obs_outbreak_data)==FALSE),
              msg="Need at least one of obs_sero_data, obs_case_data or obs_outbreak_data")
  assert_that(is.list(const_list))

  #Process input data to check that all regions with sero, case and/or outbreak data supplied are present, remove
  #regions without any supplied data, and add cross-referencing tables for use when calculating likelihood. Take
  #subset of environmental data (if used) and check that environmental data available for all regions
  input_data=input_data_process(input_data,obs_sero_data,obs_case_data,obs_outbreak_data)
  regions=names(table(input_data$region_labels)) #Regions in new processed input data list
  n_regions=length(regions)
  if(is.null(const_list$enviro_data)==FALSE){
    for(region in regions){assert_that(region %in% const_list$enviro_data$adm1)}
    enviro_data=subset(const_list$enviro_data,const_list$enviro_data$adm1 %in% regions)
  }

  frac=1.0/const_list$n_reps
  n_params=length(param_prop)
  extra_params=c()
  if(is.null(const_list$vaccine_efficacy)==TRUE){extra_params=append(extra_params,"vaccine_efficacy")}
  if(is.null(const_list$p_rep_severe)==TRUE){extra_params=append(extra_params,"p_rep_severe")}
  if(is.null(const_list$p_rep_death)==TRUE){extra_params=append(extra_params,"p_rep_death")}
  names(param_prop)=create_param_labels(const_list$type,input_data,const_list$enviro_data,extra_params)

  #Get vaccine efficacy
  if(is.numeric(const_list$vaccine_efficacy)==FALSE){
    vaccine_efficacy=exp(param_prop[names(param_prop)=="vaccine_efficacy"])
  } else {
    vaccine_efficacy=const_list$vaccine_efficacy
  }

  #Get reporting probabilities
  if(is.numeric(const_list$p_rep_severe)==FALSE){
    p_rep_severe=as.numeric(exp(param_prop[names(param_prop)=="p_rep_severe"]))
  } else {
    p_rep_severe=const_list$p_rep_severe
  }
  if(is.numeric(const_list$p_rep_death)==FALSE){
    p_rep_death=as.numeric(exp(param_prop[names(param_prop)=="p_rep_death"]))
  } else {
    p_rep_death=const_list$p_rep_death
  }

  #Get FOI and R0 values
  FOI_values=R0_values=rep(0,n_regions)
  if(const_list$type %in% c("FOI+R0 enviro","FOI enviro")){
    for(i in 1:n_regions){
      model_params=param_calc_enviro(param=param_prop,
                                     enviro_data=enviro_data[enviro_data$adm1==regions[i],])
      FOI_values[i]=as.numeric(model_params[[1]])
      if(const_list$type=="FOI+R0 enviro"){R0_values[i]=as.numeric(model_params[[2]])} else {
        R0_values[i]=const_list$R0_fixed_values[i]}
    }
  }
  if(const_list$type %in% c("FOI+R0","FOI")){
    FOI_values=exp(param_prop[c(1:n_regions)])
    if(const_list$type=="FOI+R0"){R0_values=exp(param_prop[c((n_regions+1):(2*n_regions))])
    } else {R0_values=const_list$R0_fixed_values}
  }

  #Generate modelled data over all regions
  dataset <- Generate_Dataset(input_data,FOI_values,R0_values,
                              obs_sero_data,obs_case_data,obs_outbreak_data,
                              vaccine_efficacy,p_rep_severe,p_rep_death,
                              const_list$mode_start,const_list$n_reps,const_list$dt)

  return(dataset)
}
#-------------------------------------------------------------------------------
#' @title data_match_multi
#'
#' @description Function which runs the model to create simulated data corresponding to supplied observed data for
#'   multiple parameter sets
#'
#' @details This function runs the data_match_single() function for multiple parameter sets. It takes in observed data
#'   and produces corresponding simulated datasets for comparison using multiple sets of supplied model parameters (in
#'   the param_sets input variable), population and vaccination settings (in the input_data input variable) and other
#'   settings including parameter type, time increment and environmental covariate values (in the const_list input
#'   variable). The simulated dataset will be able to be compared directly to the observed dataset - individual data
#'   values will be produced for the same years and regions (or combinations of regions), seroprevalence data will be
#'   produced for the same age groups, etc.
#'
#' @param param_sets Data frame of log values of proposed parameters, one set per row
#' @param input_data List of population and vaccination data for multiple regions, with tables to cross-reference
#'   with observed data, added using input_data_process2
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param obs_outbreak_data Outbreak Y/N data for comparison, by region and year, in format 0 = no outbreaks,
#'   1 = 1 or more outbreak(s)
#' @param const_list = List of constant parameters/flags/etc. (type,n_reps,mode_start,dt,enviro_data,R0_fixed_values,
#'   vaccine_efficacy,p_rep_severe,p_rep_death)
#'
#' @export
#'
data_match_multi <- function(param_sets=list(),input_data=list(),obs_sero_data=NULL,obs_case_data=NULL,
                             obs_outbreak_data=NULL,const_list=list()){

  assert_that(is.data.frame(param_sets))

  n_param_sets=nrow(param_sets)
  model_data_all=list()
  cat("\nSet:\n")
  for(i in 1:n_param_sets){
    cat("\t",i)
    param_prop=as.numeric(param_sets[i,])
    model_data_all[[i]] <- data_match_single(param_prop,input_data,obs_sero_data,obs_case_data,obs_outbreak_data,
                                             const_list)
  }

  return(model_data_all)
}
#-------------------------------------------------------------------------------
#' @title sero_match_graphs
#'
#' @description Function to create a series of graphs comparing modelled and observed serological data from results
#'   generated from data_match_multi() function
#'
#' @details Takes in simulated datasets produced by data_match_multi() function along with observed seroprevalence
#'   dataset used to generate the simulated datasets, and produces graph(s) (one graph per survey) showing the observed
#'   and simulated data on the same set(s) of axes, with observed data displayed as points and lines and simulated data
#'   as coloured bands.
#'
#' @param model_data Simulated datasets produced by data_match_multi() function
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param type Form in which to plot model data: "all": bands showing 95\% and 50\% of all values; "mean": bands showing
#'   95\% and 50\% confidence intervals for mean of all values
#' @param text_size1 Size of text to display on graphs
#'
#' @export
#'
sero_match_graphs <- function(model_data=list(),obs_sero_data=list(),type="mean",text_size1=1.0){

  assert_that(is.list(model_data))
  assert_that(is.data.frame(obs_sero_data))
  assert_that(type %in% c("mean","all"))
  assert_that(is.numeric(text_size1))

  n_param_sets=length(model_data)
  obs_sero_values=obs_sero_data$positives/obs_sero_data$samples
  n_sero_values=length(obs_sero_values)
  obs_sero_values[is.nan(obs_sero_values)]=0.0
  obs_sero_values_low=obs_sero_values_high=rep(NA,n_sero_values)
  for(i in 1:n_sero_values){
    if(obs_sero_data$samples[i]>0){
      CI=prop.test(x=obs_sero_data$positives[i],n=obs_sero_data$samples[i])
      obs_sero_values_low[i]=CI$conf.int[1]
      obs_sero_values_high[i]=CI$conf.int[2]
    } else {
      obs_sero_values_low[i]=0
      obs_sero_values_high[i]=0
    }
  }

  model_sero_values=array(NA,dim=c(length(obs_sero_values),n_param_sets))
  model_CI50_low=model_CI50_high=model_CI95_low=model_CI95_high=rep(0,length(obs_sero_values))
  for(i in 1:n_param_sets){
    model_sero_values[,i]=model_data[[i]]$model_sero_values
  }
  lines_all=graph_lines=c(1:length(obs_sero_values))
  if(type=="mean"){
    for(i in lines_all){
      CI_095=CI(model_sero_values[i,],ci=0.95)
      model_CI95_low[i]=CI_095[3][[1]]
      model_CI95_high[i]=CI_095[1][[1]]
      CI_050=CI(model_sero_values[i,],ci=0.50)
      model_CI50_low[i]=CI_050[3][[1]]
      model_CI50_high[i]=CI_050[1][[1]]
    }
  } else {
    n_095_low=ceiling(n_param_sets*0.025)
    n_095_high=floor(n_param_sets*0.975)
    n_050_low=ceiling(n_param_sets*0.25)
    n_050_high=floor(n_param_sets*0.75)
    for(i in lines_all){
      model_values_sorted=sort(model_sero_values[i,])
      model_CI95_low[i]=model_values_sorted[n_095_low]
      model_CI95_high[i]=model_values_sorted[n_095_high]
      model_CI50_low[i]=model_values_sorted[n_050_low]
      model_CI50_high[i]=model_values_sorted[n_050_high]
    }
  }

  if(is.null(obs_sero_data$country_zone)==FALSE){
    data_regions=names(table(obs_sero_data$country_zone))
    graph_titles=c()
    n_graphs=0
    for(region in data_regions){
      lines=lines_all[obs_sero_data$country_zone==region]
      subset=subset(obs_sero_data,obs_sero_data$country_zone==region)
      years=names(table(subset$year))
      for(year in years){
        lines2=lines[subset$year==year]
        subset2=subset(subset,year==year)
        n_graphs=n_graphs+1
        graph_lines[lines2]=n_graphs
        graph_titles=append(graph_titles,paste(region,year,sep=" "))
      }
    }
  } else {
    data_regions=names(table(obs_sero_data$adm1))
    graph_titles=c()
    n_graphs=0
    for(region in data_regions){
      lines=lines_all[obs_sero_data$adm1==region]
      subset=subset(obs_sero_data,obs_sero_data$adm1==region)
      years=names(table(subset$year))
      for(year in years){
        lines2=lines[subset$year==year]
        subset2=subset(subset,year==year)
        n_graphs=n_graphs+1
        graph_lines[lines2]=n_graphs
        graph_titles=append(graph_titles,paste(region,year,sep=" "))
      }
    }
  }

  age_values=sero_obs=sero_obs_low=sero_obs_high=sero_model_low95=NULL
  sero_model_low50=sero_model_high95=sero_model_high50=NULL
  sero_graphs=list()
  for(i in 1:n_graphs){
    lines=lines_all[graph_lines==i]

    df=data.frame(age_values=obs_sero_data$age_min[lines],sero_obs=obs_sero_values[lines],
                  samples=obs_sero_data$samples[lines],
                  sero_obs_low=obs_sero_values_low[lines],sero_obs_high=obs_sero_values_high[lines],
                  sero_model_low95=model_CI95_low[lines],sero_model_high95=model_CI95_high[lines],
                  sero_model_low50=model_CI50_low[lines],sero_model_high50=model_CI50_high[lines])
    df$samples[df$samples==0]=1

    samples=NULL
    sero_graphs[[i]] <- ggplot(data=df) + theme_bw()+labs(title=graph_titles[i])
    sero_graphs[[i]] <- sero_graphs[[i]]+geom_ribbon(data=df,aes(x=age_values,ymin=sero_model_low95,
                                                                 ymax=sero_model_high95),fill="blue",alpha=0.5)
    sero_graphs[[i]] <- sero_graphs[[i]]+geom_ribbon(data=df,aes(x=age_values,ymin=sero_model_low50,
                                                                 ymax=sero_model_high50),fill="green",alpha=0.5)
    sero_graphs[[i]] <- sero_graphs[[i]]+geom_point(data=df,aes(x=age_values,y=sero_obs,size=log(samples)),
                                                    show.legend=FALSE)
    sero_graphs[[i]] <- sero_graphs[[i]]+geom_errorbar(data=df,aes(x=age_values,ymin=sero_obs_low,ymax=sero_obs_high),
                                                       width=1.0)
    sero_graphs[[i]] <- sero_graphs[[i]]+scale_x_continuous(name="Age (min)",breaks=df$age_values,labels=df$age_values)
    sero_graphs[[i]] <- sero_graphs[[i]]+scale_y_continuous(name="Seroprevalence")
    sero_graphs[[i]] <- sero_graphs[[i]]+theme(axis.text.x = element_text(size = text_size1),
                                               axis.text.y = element_text(size = text_size1),
                                               title=element_text(size=text_size1))
  }

  return(sero_graphs)
}

#-------------------------------------------------------------------------------
#' @title case_match_graphs
#'
#' @description Function to create a series of graphs comparing modelled and observed case data from results
#'   generated from data_match_multi() function
#'
#' @details Takes in simulated datasets produced by data_match_multi() function along with observed annual case
#'   dataset used to generate the simulated datasets, and produces graph(s) (one graph per region or group of regions
#'   for which annual case data was supplied) showing the observed and simulated data on the same set(s) of axes, with
#'   observed data displayed as points and lines and simulated data as coloured bands.
#'
#' @param model_data Simulated datasets produced by data_match_multi() function
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param input_data List of population and vaccination data for multiple regions, with tables to cross-reference
#' with observed data, added using input_data_process
#' @param type Form in which to plot model data: "all": bands showing 95\% and 50\% of all values; "mean": bands showing
#'   95\% and 50\% confidence intervals for mean of all values
#' @param text_size1 Size of text to display on graphs
#'
#' @export
#'
case_match_graphs <- function(model_data=list(),obs_case_data=list(),input_data=list(),type="mean",text_size1=1.0){

  assert_that(is.list(model_data))
  assert_that(is.data.frame(obs_case_data))
  assert_that(type %in% c("mean","all"))
  assert_that(is.numeric(text_size1))

  n_param_sets=length(model_data)
  obs_case_values=obs_case_data$cases
  obs_death_values=obs_case_data$deaths
  n_case_values=length(obs_case_values)
  n_year_values=match(obs_case_data$year,input_data$years_labels)
  n_region_values=match(obs_case_data$adm1,input_data$region_labels)
  obs_case_values_low=obs_case_values_high=obs_death_values_low=obs_death_values_high=rep(NA,n_case_values)
  for(i in 1:n_case_values){
    regions=strsplit(obs_case_data$adm1[i],",")[[1]]
    n_region_values=input_data$region_labels %in% regions
    population=round(sum(input_data$pop_data[n_region_values,n_year_values[i],]),digits=0)
    CI=prop.test(x=obs_case_data$cases[i],n=population)
    obs_case_values_low[i]=round(CI$conf.int[1]*population,digits=0)
    obs_case_values_high[i]=round(CI$conf.int[2]*population,digits=0)
    CI=prop.test(x=obs_case_data$deaths[i],n=population)
    obs_death_values_low[i]=round(CI$conf.int[1]*population,digits=0)
    obs_death_values_high[i]=round(CI$conf.int[2]*population,digits=0)
  }

  model_case_values=model_death_values=array(NA,dim=c(length(obs_case_values),n_param_sets))
  for(i in 1:n_param_sets){
    model_case_values[,i]=model_data[[i]]$model_case_values
    model_death_values[,i]=model_data[[i]]$model_death_values
  }
  model_cases_CI95_low=model_cases_CI95_high=model_cases_CI50_low=model_cases_CI50_high=rep(0,n_case_values)
  model_deaths_CI95_low=model_deaths_CI95_high=model_deaths_CI50_low=model_deaths_CI50_high=rep(0,n_case_values)
  if(type=="mean"){
    for(i in 1:n_case_values){
      CI_095=CI(model_case_values[i,],ci=0.95)
      model_cases_CI95_low[i]=CI_095[3][[1]]
      model_cases_CI95_high[i]=CI_095[1][[1]]
      CI_050=CI(model_case_values[i,],ci=0.50)
      model_cases_CI50_low[i]=CI_050[3][[1]]
      model_cases_CI50_high[i]=CI_050[1][[1]]
      CI_095=CI(model_death_values[i,],ci=0.95)
      model_deaths_CI95_low[i]=CI_095[3][[1]]
      model_deaths_CI95_high[i]=CI_095[1][[1]]
      CI_050=CI(model_death_values[i,],ci=0.50)
      model_deaths_CI50_low[i]=CI_050[3][[1]]
      model_deaths_CI50_high[i]=CI_050[1][[1]]
    }
  } else {
    n_095_low=ceiling(n_param_sets*0.025)
    n_095_high=floor(n_param_sets*0.975)
    n_050_low=ceiling(n_param_sets*0.25)
    n_050_high=floor(n_param_sets*0.75)
    for(i in 1:n_case_values){
      model_case_values_sorted=sort(model_case_values[i,])
      model_cases_CI95_low[i]=model_case_values_sorted[n_095_low]
      model_cases_CI95_high[i]=model_case_values_sorted[n_095_high]
      model_cases_CI50_low[i]=model_case_values_sorted[n_050_low]
      model_cases_CI50_high[i]=model_case_values_sorted[n_050_high]
      model_death_values_sorted=sort(model_death_values[i,])
      model_deaths_CI95_low[i]=model_death_values_sorted[n_095_low]
      model_deaths_CI95_high[i]=model_death_values_sorted[n_095_high]
      model_deaths_CI50_low[i]=model_death_values_sorted[n_050_low]
      model_deaths_CI50_high[i]=model_death_values_sorted[n_050_high]
    }
  }

  data_regions=names(table(obs_case_data$adm1))
  n_graphs=length(data_regions)

  cases_graphs=deaths_graphs=list()
  years=case_obs=case_obs_low=case_obs_high=case_model_low95=case_model_low50=case_model_high95=case_model_high50=NULL
  death_obs=death_obs_low=death_obs_high=death_model_low95=death_model_low50=death_model_high95=death_model_high50=NULL
  for(i in 1:n_graphs){
    region=data_regions[i]
    lines=obs_case_data$adm1==region

    df=data.frame(years=obs_case_data$year[lines],case_obs=obs_case_values[lines],death_obs=obs_death_values[lines],
                  case_obs_low=obs_case_values_low[lines],case_obs_high=obs_case_values_high[lines],
                  case_model_low95=model_cases_CI95_low[lines],case_model_high95=model_cases_CI95_high[lines],
                  case_model_low50=model_cases_CI50_low[lines],case_model_high50=model_cases_CI50_high[lines],
                  death_obs_low=obs_death_values_low[lines],death_obs_high=obs_death_values_high[lines],
                  death_model_low95=model_deaths_CI95_low[lines],death_model_high95=model_deaths_CI95_high[lines],
                  death_model_low50=model_deaths_CI50_low[lines],death_model_high50=model_deaths_CI50_high[lines])

    cases_graphs[[i]] <- ggplot(data=df) + theme_bw()+labs(title=substr(region,1,5))
    cases_graphs[[i]] <- cases_graphs[[i]]+geom_ribbon(data=df,aes(x=years,ymin=case_model_low95,
                                                                   ymax=case_model_high95),fill="blue",alpha=0.5)
    cases_graphs[[i]] <- cases_graphs[[i]]+geom_ribbon(data=df,aes(x=years,ymin=case_model_low50,
                                                                   ymax=case_model_high50),fill="green",alpha=0.5)
    cases_graphs[[i]] <- cases_graphs[[i]]+geom_point(data=df,aes(x=years,y=case_obs))
    cases_graphs[[i]] <- cases_graphs[[i]]+geom_errorbar(data=df,aes(x=years,ymin=case_obs_low,ymax=case_obs_high),
                                                         width=0.5)
    cases_graphs[[i]] <- cases_graphs[[i]]+scale_x_continuous(name="",breaks=c(min(df$years):max(df$years)),
                                                              labels=c(min(df$years):max(df$years)))
    cases_graphs[[i]] <- cases_graphs[[i]]+scale_y_continuous(name="Cases")
    cases_graphs[[i]] <- cases_graphs[[i]]+theme(axis.text.x = element_text(size = text_size1),
                                                 axis.text.y = element_text(size = text_size1),
                                                 title=element_text(size=text_size1))


    deaths_graphs[[i]] <- ggplot(data=df) + theme_bw()+labs(title=substr(region,1,5))
    deaths_graphs[[i]] <- deaths_graphs[[i]]+geom_ribbon(data=df,aes(x=years,ymin=death_model_low95,
                                                                     ymax=death_model_high95),fill="blue",alpha=0.5)
    deaths_graphs[[i]] <- deaths_graphs[[i]]+geom_ribbon(data=df,aes(x=years,ymin=death_model_low50,
                                                                     ymax=death_model_high50),fill="green",alpha=0.5)
    deaths_graphs[[i]] <- deaths_graphs[[i]]+geom_point(data=df,aes(x=years,y=death_obs))
    deaths_graphs[[i]] <- deaths_graphs[[i]]+geom_errorbar(data=df,aes(x=years,ymin=death_obs_low,ymax=death_obs_high),
                                                           width=0.5)
    deaths_graphs[[i]] <- deaths_graphs[[i]]+scale_x_continuous(name="",breaks=c(min(df$years):max(df$years)),
                                                                labels=c(min(df$years):max(df$years)))
    deaths_graphs[[i]] <- deaths_graphs[[i]]+scale_y_continuous(name="Deaths")
    deaths_graphs[[i]] <- deaths_graphs[[i]]+theme(axis.text.x = element_text(size = text_size1),
                                                   axis.text.y = element_text(size = text_size1),
                                                   title=element_text(size=text_size1))
  }

  return(list(cases_graphs=cases_graphs,deaths_graphs=deaths_graphs))
}
