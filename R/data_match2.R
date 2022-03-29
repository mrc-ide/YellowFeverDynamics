#Alternative data-match graph plotting functions, using ggplot instead of matplot
#-------------------------------------------------------------------------------
#' @title sero_match_graphs2
#'
#' @description Function to create a series of graphs comparing modelled and observed serological data from results
#'   generated from data_match_multi
#'
#' @details TBA
#'
#' @param model_data TBA
#' @param obs_sero_data TBA
#'
#' @export
#'
sero_match_graphs2 <- function(model_data=list(),obs_sero_data=list()){
  #TODO - Add assert_that functions

  n_param_sets=length(model_data)
  obs_sero_values=obs_sero_data$positives/obs_sero_data$samples
  n_sero_values=length(obs_sero_values)
  obs_sero_values[is.nan(obs_sero_values)]=0.0
  obs_sero_values_low=obs_sero_values_high=rep(NA,n_sero_values)
  for(i in 1:n_sero_values){
    CI=prop.test(x=obs_sero_data$positives[i],n=obs_sero_data$samples[i])
    obs_sero_values_low[i]=CI$conf.int[1]
    obs_sero_values_high[i]=CI$conf.int[2]
  }

  model_sero_values=array(NA,dim=c(length(obs_sero_values),n_param_sets))
  model_CI95=array(NA,dim=c(3,length(obs_sero_values)))
  model_CI50=model_CI95
  for(i in 1:n_param_sets){
    model_sero_values[,i]=model_data[[i]]$model_sero_data$sero
  }
  lines_all=graph_lines=c(1:length(obs_sero_values))
  for(i in lines_all){
    model_CI95[,i]=CI(model_sero_values[i,],ci=0.95)
    model_CI50[,i]=CI(model_sero_values[i,],ci=0.50)
  }

  if(is.null(obs_sero_data$country_zone)==FALSE){
    data_regions=names(table(obs_sero_data$country_zone))
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
      }
    }
  } else {
    data_regions=names(table(obs_sero_data$adm1))
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
      }
    }
  }

  n_graphs=length(data_regions)
  sero_graphs=list()
  for(i in 1:n_graphs){
    lines=lines_all[graph_lines==i]
    if(is.null(obs_sero_data$country_zone)==FALSE){
      region=obs_sero_data$country_zone[lines[1]]
    } else{
      region=obs_sero_data$adm1[lines[1]]
    }
    age_values=obs_sero_data$age_min[lines]

    df=data.frame(age_values=obs_sero_data$age_min[lines],sero_obs=obs_sero_values[lines],
                  sero_obs_low=obs_sero_values_low[lines],sero_obs_high=obs_sero_values_high[lines])
    df2=data.frame(age_values=c(df$age_values,rev(df$age_values)),
                  sero_model_CI95=c(model_CI95[1,lines],rev(model_CI95[3,lines])),
                  sero_model_CI50=c(model_CI50[1,lines],rev(model_CI50[3,lines])))

    sero_graphs[[i]] <- ggplot(data=df) + theme_bw()+labs(title=region)
    sero_graphs[[i]] <- sero_graphs[[i]]+geom_polygon(data=df2,aes(x=age_values,y=sero_model_CI95),fill="blue")
    sero_graphs[[i]] <- sero_graphs[[i]]+geom_polygon(data=df2,aes(x=age_values,y=sero_model_CI50),fill="green")
    sero_graphs[[i]] <- sero_graphs[[i]]+geom_line(data=df,aes(x=age_values,y=sero_obs))
    sero_graphs[[i]] <- sero_graphs[[i]]+geom_errorbar(data=df,aes(x=age_values,ymin=sero_obs_low,ymax=sero_obs_high),
                                                       width=1.0)
    sero_graphs[[i]] <- sero_graphs[[i]]+scale_x_continuous(name="",breaks=df$age_values,labels=df$age_values)
    sero_graphs[[i]] <- sero_graphs[[i]]+scale_y_continuous(name="")
    sero_graphs[[i]] <- sero_graphs[[i]]+theme(text = element_text(size = 8))
  }

  return(sero_graphs)
}

#-------------------------------------------------------------------------------
#' @title case_match_graphs2
#'
#' @description Function to create a series of graphs comparing modelled and observed case data from results
#'   generated from data_match_multi
#'
#' @details TBA
#'
#' @param model_data TBA
#' @param obs_case_data TBA
#' @param input_data TBA
#'
#' @export
#'
case_match_graphs2 <- function(model_data=list(),obs_case_data=list(),input_data=list()){
  #TODO - Add assert_that functions

  n_param_sets=length(model_data)
  obs_case_values=obs_case_data$cases
  obs_death_values=obs_case_data$deaths
  n_case_values=length(obs_case_values)
  n_year_values=match(obs_case_data$year,input_data$years_labels)
  n_region_values=match(obs_case_data$adm1,input_data$region_labels)
  obs_case_values_low=obs_case_values_high=obs_death_values_low=obs_death_values_high=rep(NA,n_case_values)
  for(i in 1:n_case_values){
    population=round(sum(input_data$pop_data[n_region_values[i],n_year_values[i],]),digits=0)
    CI=prop.test(x=obs_case_data$cases[i],n=population)
    obs_case_values_low[i]=round(CI$conf.int[1]*population,digits=0)
    obs_case_values_high[i]=round(CI$conf.int[2]*population,digits=0)
    CI=prop.test(x=obs_case_data$deaths[i],n=population)
    obs_death_values_low[i]=round(CI$conf.int[1]*population,digits=0)
    obs_death_values_high[i]=round(CI$conf.int[2]*population,digits=0)
  }

  model_case_values=model_death_values=array(NA,dim=c(length(obs_case_values),n_param_sets))
  model_CI_cases95=model_CI_deaths95=model_CI_cases50=model_CI_deaths50=array(NA,dim=c(3,length(obs_case_values)))
  for(i in 1:n_param_sets){
    model_case_values[,i]=model_data[[i]]$model_case_data$cases
    model_death_values[,i]=model_data[[i]]$model_case_data$deaths
  }
  for(i in 1:n_case_values){
    model_CI_cases95[,i]=CI(model_case_values[i,],ci=0.95)
    model_CI_deaths95[,i]=CI(model_death_values[i,],ci=0.95)
    model_CI_cases50[,i]=CI(model_case_values[i,],ci=0.50)
    model_CI_deaths50[,i]=CI(model_death_values[i,],ci=0.50)
  }

  data_regions=names(table(obs_case_data$adm1))
  n_graphs=length(data_regions)

  cases_graphs=deaths_graphs=list()
  for(i in 1:n_graphs){
    region=data_regions[i]
    lines=obs_case_data$adm1==region

    df=data.frame(years=obs_case_data$year[lines],case_obs=obs_case_values[lines],death_obs=obs_death_values[lines],
                  case_obs_low=obs_case_values_low[lines],case_obs_high=obs_case_values_high[lines],
                  death_obs_low=obs_death_values_low[lines],death_obs_high=obs_death_values_high[lines])
    df2=data.frame(years=c(df$years,rev(df$years)),
                   case_model_CI95=c(model_CI_cases95[1,lines],rev(model_CI_cases95[3,lines])),
                   case_model_CI50=c(model_CI_cases50[1,lines],rev(model_CI_cases50[3,lines])),
                   death_model_CI95=c(model_CI_deaths95[1,lines],rev(model_CI_deaths95[3,lines])),
                   death_model_CI50=c(model_CI_deaths50[1,lines],rev(model_CI_deaths50[3,lines])))

    cases_graphs[[i]] <- ggplot(data=df) + theme_bw()+labs(title=region)
    cases_graphs[[i]] <- cases_graphs[[i]]+geom_polygon(data=df2,aes(x=years,y=case_model_CI95),fill="blue")
    cases_graphs[[i]] <- cases_graphs[[i]]+geom_polygon(data=df2,aes(x=years,y=case_model_CI50),fill="green")
    cases_graphs[[i]] <- cases_graphs[[i]]+geom_line(data=df,aes(x=years,y=case_obs))
    cases_graphs[[i]] <- cases_graphs[[i]]+geom_errorbar(data=df,aes(x=years,ymin=case_obs_low,ymax=case_obs_high),
                                                       width=0.5)
    cases_graphs[[i]] <- cases_graphs[[i]]+scale_x_continuous(name="",breaks=df$years,labels=df$years)
    cases_graphs[[i]] <- cases_graphs[[i]]+scale_y_continuous(name="")
    cases_graphs[[i]] <- cases_graphs[[i]]+theme(text = element_text(size = 8))


    deaths_graphs[[i]] <- ggplot(data=df) + theme_bw()+labs(title=region)
    deaths_graphs[[i]] <- deaths_graphs[[i]]+geom_polygon(data=df2,aes(x=years,y=death_model_CI95),fill="blue")
    deaths_graphs[[i]] <- deaths_graphs[[i]]+geom_polygon(data=df2,aes(x=years,y=death_model_CI50),fill="green")
    deaths_graphs[[i]] <- deaths_graphs[[i]]+geom_line(data=df,aes(x=years,y=death_obs))
    deaths_graphs[[i]] <- deaths_graphs[[i]]+geom_errorbar(data=df,aes(x=years,ymin=death_obs_low,ymax=death_obs_high),
                                                         width=0.5)
    deaths_graphs[[i]] <- deaths_graphs[[i]]+scale_x_continuous(name="",breaks=df$years,labels=df$years)
    deaths_graphs[[i]] <- deaths_graphs[[i]]+scale_y_continuous(name="")
    deaths_graphs[[i]] <- deaths_graphs[[i]]+theme(text = element_text(size = 8))
  }

  return(list(cases_graphs=cases_graphs,deaths_graphs=deaths_graphs))
}
