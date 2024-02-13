#-------------------------------------------------------------------------------
#' @title Generate_Dataset_VarFR
#'
#' @description Generate annual serological and/or case/death data with annually varying FOI/R0
#'
#' @details This function is used to generate annual serological and/or case/death data based on templates;
#' it is a variant of the Generate_Dataset function in the YEP package, identical except for using the
#' Model_Run_VarFR function instead of the YEP::Model_Run function. Inputs are the same as YEP::Generate_Dataset
#' except that FOI_values and R0_values are 2-dimensional arrays (dimensions no. regions x no. years). (Parallel
#' processing functionality has also been removed)
#'
#' @param input_data List of population and vaccination data for multiple regions in standard format [TBA]
#' @param FOI_values Array of values of the force of infection due to spillover from sylvatic reservoir by region and year
#' @param R0_values Array of values of the basic reproduction number for human-human transmission by region and year
#' @param sero_template Seroprevalence data template - data frame with region, year, minimum/maximum age, vc_factor [TBA]
#'   and number of samples
#' @param case_template Annual reported case/death data template - data frame with region and year
#' @param vaccine_efficacy Fractional vaccine efficacy
#' @param p_severe_inf Probability of an infection being severe
#' @param p_death_severe_inf Probability of a severe infection resulting in death
#' @param p_rep_severe Probability of reporting of a severe but non-fatal infection
#' @param p_rep_death Probability of reporting of a fatal infection
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity (uniform by age, R0 based only)
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#'  If mode_start=3, shift some non-vaccinated individuals into recovered to give herd immunity (stratified by age)
#' @param start_SEIRV SEIRV data from end of a previous run to use as input (list of datasets, one per region)
#' @param dt Time increment in days to use in model (should be either 1.0, 2.5 or 5.0 days)
#' @param n_reps number of stochastic repetitions
#' @param deterministic TRUE/FALSE - set model to run in deterministic mode if TRUE
#' @param output_frame Flag indicating whether to output a complete data frame of results in template format (if TRUE)
#'   or calculated values only (if FALSE)
#' '
#' @export
#'
Generate_Dataset_VarFR <- function(input_data = list(),FOI_values = c(),R0_values = c(),sero_template = NULL,case_template = NULL,
                                   vaccine_efficacy = 1.0, p_severe_inf = 0.12, p_death_severe_inf = 0.39, p_rep_severe = 1.0,
                                   p_rep_death = 1.0,mode_start = 1,start_SEIRV = NULL, dt = 1.0,n_reps = 1, deterministic = FALSE,
                                   output_frame=FALSE){

  assert_that(input_data_check(input_data),msg=paste("Input data must be in standard format",
                                                     " (see https://mrc-ide.github.io/YEP/articles/CGuideAInputs.html)"))
  assert_that(any(is.null(sero_template)==FALSE,is.null(case_template)==FALSE),msg="Need at least one template")
  if(is.null(sero_template)==FALSE){
    assert_that(all(c("region","year","age_min","age_max","samples","vc_factor") %in% names(sero_template)))
  }
  if(is.null(case_template)==FALSE){
    assert_that(all(c("region","year") %in% names(case_template)))
    assert_that(p_severe_inf>=0.0 && p_severe_inf<=1.0,msg="Severe infection rate must be between 0-1")
    assert_that(p_death_severe_inf>=0.0 && p_death_severe_inf<=1.0,
                msg="Fatality rate of severe infections must be between 0-1")
    assert_that(p_rep_severe>=0.0 && p_rep_severe<=1.0,msg="Severe infection reporting probability must be between 0-1")
    assert_that(p_rep_death>=0.0 && p_rep_death<=1.0,msg="Fatal infection reporting probability must be between 0-1")
  }

  #Prune input data based on regions
  regions=regions_breakdown(c(sero_template$region,case_template$region))
  input_data=input_data_truncate(input_data,regions)
  n_regions=length(input_data$region_labels)

  #Cross-reference templates with input regions
  if(is.null(sero_template)==FALSE){
    xref_sero=template_region_xref(sero_template,input_data$region_labels)
    sero_line_list=xref_sero$line_list
  } else {
    xref_sero=data.frame(year_data_begin=rep(Inf,n_regions),year_end=rep(-Inf,n_regions))
    sero_line_list=rep(NA,n_regions)
  }
  if(is.null(case_template)==FALSE){
    xref_case=template_region_xref(case_template,input_data$region_labels)
    case_line_list=xref_case$line_list
  } else {
    xref_case=data.frame(year_data_begin=rep(Inf,n_regions),year_end=rep(-Inf,n_regions))
    case_line_list=rep(NA,n_regions)
  }
  year_data_begin=year_end=rep(NA,length(input_data$region_labels))
  for(i in 1:length(year_data_begin)){
    year_data_begin[i]=min(xref_sero$year_data_begin[i],xref_case$year_data_begin[i])
    year_end[i]=max(xref_sero$year_end[i],xref_case$year_end[i])
  }

  inv_reps=1/n_reps
  assert_that(length(FOI_values)==n_regions,msg="Length of FOI_values must match number of regions to be modelled")
  assert_that(length(R0_values)==n_regions,msg="Length of R0_values must match number of regions to be modelled")
  if(mode_start==2){assert_that(length(start_SEIRV)==n_regions,
                                msg="Number of start_SEIRV datasets must match number of regions")}

  #Set up data structures to take modelled data corresponding to observed data
  if(is.null(sero_template)){model_sero_data=NULL} else {
    blank=rep(0,nrow(sero_template))
    model_sero_data=data.frame(samples=blank,positives=blank,sero=blank)
  }
  if(is.null(case_template)){model_case_values=model_death_values=NA} else {
    model_case_values=model_death_values=rep(0,nrow(case_template))
  }

  #Set up vector of output types to get from model
  output_types=rep(NA,n_regions)
  for(n_region in 1:n_regions){
    if(is.na(case_line_list[[n_region]][1])==FALSE){
      if(is.na(sero_line_list[[n_region]][1])==FALSE){output_types[n_region]="case+sero"} else{output_types[n_region]="case"}
    } else {output_types[n_region]="sero"}
  }

  i
  #Save relevant output data from each region
  for(n_region in 1:n_regions){

    #Run model
    #cat("\n\t\tBeginning modelling region ",input_data$region_labels[n_region])
    model_output = Model_Run_VarFR(FOI_spillover = FOI_values[n_region,],R0 = R0_values[n_region,],
                                   vacc_data = input_data$vacc_data[n_region,,],pop_data = input_data$pop_data[n_region,,],
                                   years_data = c(year_data_begin[n_region]:year_end[n_region]),
                                   start_SEIRV=start_SEIRV[[n_region]],output_type = output_types[n_region],
                                   year0 = input_data$years_labels[1],mode_start = mode_start,
                                   vaccine_efficacy = vaccine_efficacy, dt = dt, n_particles = n_reps,n_threads = n_reps,
                                   deterministic = deterministic)
    #cat("\n\t\tFinished modelling region ",n_region)
    t_pts=length(model_output$year)

    #Compile case data if needed
    if(is.na(case_line_list[[n_region]][1])==FALSE){
      case_line_list_region=case_line_list[[n_region]]
      years_case=case_template$year[case_line_list_region]
      n_lines=length(case_line_list_region)

      for(n_rep in 1:n_reps){
        rep_cases=rep_deaths=rep(0,n_lines)
        for(n_line in 1:n_lines){
          pts=c(1:t_pts)[model_output$year==years_case[n_line]]
          infs=sum(model_output$C[n_rep,pts])
          if(deterministic){
            severe_infs=floor(infs)*p_severe_inf
            deaths=severe_infs*p_death_severe_inf
            rep_deaths[n_line]=round(deaths*p_rep_death)
            rep_cases[n_line]=rep_deaths[n_line]+round((severe_infs-deaths)*p_rep_severe)

          } else {
            severe_infs=rbinom(1,floor(infs),p_severe_inf)
            deaths=rbinom(1,severe_infs,p_death_severe_inf)
            rep_deaths[n_line]=rbinom(1,deaths,p_rep_death)
            rep_cases[n_line]=rep_deaths[n_line]+rbinom(1,floor(severe_infs-deaths),p_rep_severe)
          }
        }

        model_case_values[case_line_list_region]=model_case_values[case_line_list_region]+rep_cases
        model_death_values[case_line_list_region]=model_death_values[case_line_list_region]+rep_deaths
      }
    }

    #Compile seroprevalence data if necessary
    if(is.na(sero_line_list[[n_region]][1])==FALSE){
      sero_line_list_region=sero_line_list[[n_region]]
      for(n_rep in 1:n_reps){
        sero_results=sero_calculate2(sero_template[sero_line_list_region,],model_output,n_rep)
        model_sero_data$samples[sero_line_list_region]=model_sero_data$samples[sero_line_list_region]+sero_results$samples
        model_sero_data$positives[sero_line_list_region]=model_sero_data$positives[sero_line_list_region]+sero_results$positives
      }
    }
  }

  if(is.null(sero_template)==FALSE){model_sero_data$sero=model_sero_data$positives/model_sero_data$samples}
  if(is.null(case_template)==FALSE){
    model_case_values=model_case_values*inv_reps
    model_death_values=model_death_values*inv_reps
  }

  if(output_frame) { #Output complete frames of data
    return(list(model_sero_data=data.frame(region=sero_template$region,year=sero_template$year,
                                           age_min=sero_template$age_min,age_max=sero_template$age_max,
                                           samples=sero_template$samples,positives=sero_template$samples*model_sero_data$sero),
                model_case_data=data.frame(region=case_template$region,year=case_template$year,
                                           cases=model_case_values,deaths=model_death_values)))
  } else { #Minimal output for MCMC
    return(list(model_sero_values=model_sero_data$sero,model_case_values=model_case_values,
                model_death_values=model_death_values))
  }
}
