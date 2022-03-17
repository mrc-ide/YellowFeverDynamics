#-------------------------------------------------------------------------------
#' @title data_match_single
#'
#' @description Function which runs the model to create simulated data corresponding to supplied observed data for a
#'   single set of parameters with one or more repetitions
#'
#' @details TBA
#'
#' @param param_prop Log values of proposed parameters
#' @param input_data List of population and vaccination data for multiple regions, with tables to cross-reference
#' with observed data, added using input_data_process2
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

  enviro_data=const_list$enviro_data
  regions=input_data$region_labels
  n_regions=length(regions)
  p_severe=0.12
  p_death_severe=0.47
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

  #Get reporting probabilities and check they are within specified bounds
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
      model_params=param_calc_enviro(param=param_prop,enviro_data=enviro_data[enviro_data$adm1==regions[i],])
      FOI_values[i]=model_params$FOI
      if(const_list$type=="FOI+R0 enviro"){R0_values[i]=model_params$R0} else {R0_values[i]=const_list$R0_fixed_values[i]}
    }
  }
  if(const_list$type %in% c("FOI+R0","FOI")){
    FOI_values=exp(param_prop[c(1:n_regions)])
    if(const_list$type=="FOI+R0"){R0_values=exp(param_prop[c((n_regions+1):(2*n_regions))])
    } else {R0_values=const_list$R0_fixed_values}
  }

  #Set up data structures to take modelled data corresponding to observed data
  if(is.null(obs_sero_data)==FALSE){
    regions_sero=names(table(obs_sero_data$gadm36))
    model_sero_data=list()
    blank1=rep(0,nrow(obs_sero_data))
    blank2=data.frame(gadm36=obs_sero_data$gadm36,year=obs_sero_data$year,age_min=obs_sero_data$age_min,
                      age_max=obs_sero_data$age_max,samples=blank1,positives=blank1,sero=blank1)
    for(rep in 1:const_list$n_reps){
      model_sero_data[[rep]]=blank2
    }
  } else {
    model_sero_data=NULL
    regions_sero=NULL
    }
  if(is.null(obs_case_data)==FALSE){
    regions_case=names(table(obs_case_data$region))
    model_case_data=list()
    blank1=rep(0,nrow(obs_case_data))
    blank2=data.frame(region=obs_case_data$region,year=obs_case_data$year,cases=blank1,deaths=blank1)
    for(rep in 1:const_list$n_reps){
      model_case_data[[rep]]=blank2
    }
  } else {
    regions_case=NULL
    model_case_data=NULL
    }
  if(is.null(obs_outbreak_data)==FALSE){
    regions_outbreak=names(table(obs_outbreak_data$region))
    model_outbreak_data=list()
    blank1=data.frame(region=obs_outbreak_data$region,year=obs_outbreak_data$year,outbreak_yn=rep(0,nrow(obs_outbreak_data)))
    for(rep in 1:const_list$n_reps){
      model_outbreak_data[[rep]]=blank1
    }
  } else {
    regions_outbreak=NULL
    model_outbreak_data=NULL
    }

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
          severe_infs=rbinom(1,floor(infs),p_severe)
          deaths=rbinom(1,severe_infs,p_death_severe)
          annual_data$rep_deaths[i,n_year]=rbinom(1,deaths,p_rep_death)
          annual_data$rep_cases[i,n_year]=annual_data$rep_deaths[i,n_year]+rbinom(1,severe_infs-deaths,
                                                                                  p_rep_severe)
        }
      }

      if(flag_case==1){
        for(rep in 1:const_list$n_reps){
          model_case_data[[rep]]$cases[case_line_list]=model_case_data[[rep]]$cases[case_line_list]+annual_data$rep_cases[rep,]
          model_case_data[[rep]]$deaths[case_line_list]=model_case_data[[rep]]$deaths[case_line_list]+annual_data$rep_deaths[rep,]
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
          model_outbreak_data$outbreak_yn[input_data$outbreak_line_list[[n_region]]]=outbreak_risk
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
        model_sero_data[[i]]$samples[sero_line_list]=model_sero_data[[i]]$samples[sero_line_list]+sero_results$samples
        model_sero_data[[i]]$positives[sero_line_list]=model_sero_data[[i]]$positives[sero_line_list]+sero_results$positives
      }
    }
    model_output<-NULL
  }

  for(i in 1:const_list$n_reps){
    if(is.null(model_sero_data)==FALSE){
      model_sero_data[[i]]$sero=model_sero_data[[i]]$positives/model_sero_data[[i]]$samples
      model_sero_data[[i]]$samples=obs_sero_data$samples
      for(j in 1:nrow(obs_sero_data)){
        model_sero_data[[i]]$positives[j]=rbinom(1,obs_sero_data$samples[j],model_sero_data[[i]]$sero[j])
      }
    }
    if(is.null(model_outbreak_data)==FALSE){
      for(j in 1:nrow(obs_outbreak_data)){
        model_outbreak_data[[i]]$outbreak_yn[j]=rbinom(1,1,model_outbreak_data[[i]]$outbreak_yn[j])
      }
    }
  }

  return(list(regions_case=regions_case,model_case_data=model_case_data,
              regions_outbreak=regions_outbreak,model_outbreak_data=model_outbreak_data,
              regions_sero=regions_sero,model_sero_data=model_sero_data))
}
#-------------------------------------------------------------------------------
#' @title data_match_multi
#'
#' @description Function which runs the model to create simulated data corresponding to supplied observed data for
#'   multiple parameter sets
#'
#' @details TBA
#'
#' @param param_sets Data frame of log values of proposed parameters, one set per row
#' @param input_data List of population and vaccination data for multiple regions, with tables to cross-reference
#' with observed data, added using input_data_process2
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
  #TODO - Add assert_that functions

  n_param_sets=nrow(param_sets)
  frac=1.0/const_list$n_reps
  model_data_all=list()
  cat("\nSet:\n")
  for(i in 1:n_param_sets){
    cat("\t",i)
    param_prop=as.numeric(param_sets[i,])
    model_data <- data_match_single(param_prop,input_data,obs_sero_data,obs_case_data,obs_outbreak_data,const_list)
    model_data_all[[i]]=list(regions_case=model_data$regions_case,model_case_data=model_data$model_case_data[[1]],
                             regions_outbreak=model_data$regions_outbreak,model_outbreak_data=model_data$model_outbreak_data[[1]],
                             regions_sero=model_data$regions_sero,model_sero_data=model_data$model_sero_data[[1]])

    if(is.null(obs_sero_data)==FALSE){
      model_data_all[[i]]$model_sero_data$positives=model_data_all[[i]]$model_sero_data$sero=0
      for(rep in 1:const_list$n_reps){
        model_data_all[[i]]$model_sero_data$positives=model_data_all[[i]]$model_sero_data$positives+model_data$model_sero_data[[rep]]$positives
      }
      model_data_all[[i]]$model_sero_data$positives=model_data_all[[i]]$model_sero_data$positives*frac
      if(model_data_all[[i]]$model_sero_data$samples==0){
        model_data_all[[i]]$model_sero_data$sero=0
      } else {
        model_data_all[[i]]$model_sero_data$sero=model_data_all[[i]]$model_sero_data$positives/model_data_all[[i]]$model_sero_data$samples
      }
    }
    if(is.null(obs_case_data)==FALSE){
      model_data_all[[i]]$model_case_data$cases=model_data_all[[i]]$model_case_data$deaths=0
      for(rep in 1:const_list$n_reps){
        model_data_all[[i]]$model_case_data$cases=model_data_all[[i]]$model_case_data$cases+model_data$model_case_data[[rep]]$cases
        model_data_all[[i]]$model_case_data$deaths=model_data_all[[i]]$model_case_data$deaths+model_data$model_case_data[[rep]]$deaths
      }
      model_data_all[[i]]$model_case_data$cases=model_data_all[[i]]$model_case_data$cases*frac
      model_data_all[[i]]$model_case_data$deaths=model_data_all[[i]]$model_case_data$deaths*frac
    }
    if(is.null(obs_outbreak_data)==FALSE){
      #TODO
    }
  }

  return(model_data_all)
}
#-------------------------------------------------------------------------------
#' @title sero_match_graphs
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
sero_match_graphs <- function(model_data=list(),obs_sero_data=list()){
  #TODO - Add assert_that functions

  n_param_sets=length(model_data)
  obs_sero_values=obs_sero_data$positives/obs_sero_data$samples
  obs_sero_values[is.nan(obs_sero_values)]=0.0
  model_sero_values=array(NA,dim=c(length(obs_sero_values),n_param_sets))
  model_CI=array(NA,dim=c(3,length(obs_sero_values)))
  for(i in 1:n_param_sets){
    model_sero_values[,i]=model_data[[i]]$model_sero_data$sero
  }
  lines_all=graph_lines=c(1:length(obs_sero_values))
  for(i in lines_all){
    model_CI[,i]=CI(model_sero_values[i,])
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
    data_regions=names(table(obs_sero_data$gadm36))
    n_graphs=0
    for(region in data_regions){
      lines=lines_all[obs_sero_data$gadm36==region]
      subset=subset(obs_sero_data,obs_sero_data$gadm36==region)
      years=names(table(subset$year))
      for(year in years){
        lines2=lines[subset$year==year]
        subset2=subset(subset,year==year)
        n_graphs=n_graphs+1
        graph_lines[lines2]=n_graphs
      }
    }
  }

  nrows=floor(sqrt(n_graphs))
  ncols=ceiling(n_graphs/nrows)

  par(mfrow=c(nrows,ncols),mar=c(2,2,1,0))
  for(i in 1:n_graphs){
    lines=lines_all[graph_lines==i]
    if(is.null(obs_sero_data$country_zone)==FALSE){
      region=obs_sero_data$country_zone[lines[1]]
    } else{
      region=obs_sero_data$gadm36[lines[1]]
    }
    matplot(x=obs_sero_data$age_min[lines],y=obs_sero_values[lines],type="b",pch=16,col=1,
            ylim=c(0,max(obs_sero_values[lines],model_CI[1,lines],na.rm=TRUE)),cex.lab=0.5,cex.axis=0.5)
    matplot(x=obs_sero_data$age_min[lines],y=model_CI[2,lines],type="l",pch=1,col=2,add=TRUE)
    matplot(x=obs_sero_data$age_min[lines],y=model_CI[1,lines],type="l",pch=1,col=2,lty=2,add=TRUE)
    matplot(x=obs_sero_data$age_min[lines],y=model_CI[3,lines],type="l",pch=1,col=2,lty=2,add=TRUE)
    title(main=region,cex.main=0.5)
  }
  par(mfrow=c(1,1))

}

#-------------------------------------------------------------------------------
#' @title case_match_graphs
#'
#' @description Function to create a series of graphs comparing modelled and observed case data from results
#'   generated from data_match_multi
#'
#' @details TBA
#'
#' @param model_data TBA
#' @param obs_case_data TBA
#'
#' @export
#'
case_match_graphs <- function(model_data=list(),obs_case_data=list()){
  #TODO - Add assert_that functions

  n_param_sets=length(model_data)
  obs_case_values=obs_case_data$cases
  obs_death_values=obs_case_data$deaths
  model_case_values=model_death_values=array(NA,dim=c(length(obs_case_values),n_param_sets))
  model_CI_cases=model_CI_deaths=array(NA,dim=c(3,length(obs_case_values)))
  for(i in 1:n_param_sets){
    model_case_values[,i]=model_data[[i]]$model_case_data$cases
    model_death_values[,i]=model_data[[i]]$model_case_data$deaths
  }
  for(i in 1:length(obs_case_values)){
    model_CI_cases[,i]=CI(model_case_values[i,])
    model_CI_deaths[,i]=CI(model_death_values[i,])
  }

  data_regions=names(table(obs_case_data$region))
  n_graphs=length(data_regions)
  nrows=floor(sqrt(n_graphs))
  ncols=ceiling(n_graphs/nrows)

  par(mfrow=c(nrows,ncols),mar=c(2,2,1,0))
  for(region in data_regions){
    lines=obs_case_data$region==region
    matplot(x=obs_case_data$year[lines],y=obs_case_values[lines],type="b",pch=16,col=1,
            ylim=c(0,max(obs_case_values[lines],model_CI_cases[1,lines])),cex.lab=0.5,cex.axis=0.5)
    matplot(x=obs_case_data$year[lines],y=model_CI_cases[2,lines],type="l",pch=1,col=2,add=TRUE)
    matplot(x=obs_case_data$year[lines],y=model_CI_cases[1,lines],type="l",pch=1,col=2,lty=2,add=TRUE)
    matplot(x=obs_case_data$year[lines],y=model_CI_cases[3,lines],type="l",pch=1,col=2,lty=2,add=TRUE)
    title(main=region,cex.main=0.5)
  }
  par(mfrow=c(1,1))


  par(mfrow=c(nrows,ncols),mar=c(2,2,1,0))
  for(region in data_regions){
    lines=obs_case_data$region==region
    matplot(x=obs_case_data$year[lines],y=obs_death_values[lines],type="b",pch=16,col=1,
            ylim=c(0,max(obs_death_values[lines],model_CI_deaths[1,lines])),cex.lab=0.5,cex.axis=0.5)
    matplot(x=obs_case_data$year[lines],y=model_CI_deaths[2,lines],type="l",pch=1,col=2,add=TRUE)
    matplot(x=obs_case_data$year[lines],y=model_CI_deaths[1,lines],type="l",pch=1,col=2,lty=2,add=TRUE)
    matplot(x=obs_case_data$year[lines],y=model_CI_deaths[3,lines],type="l",pch=1,col=2,lty=2,add=TRUE)
    title(main=region,cex.main=0.5)
  }
  par(mfrow=c(1,1))

}
