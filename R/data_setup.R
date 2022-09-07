# R file for functions relating to the creation/formatting of input data for the yellow fever model
#-------------------------------------------------------------------------------
#' @title create_input_data
#'
#' @description TBA
#'
#' @details TBA
#'
#' @param vacc_data TBA
#' @param pop_data TBA
#' '
#' @export
#'
create_input_data <- function(vacc_data=list(),pop_data=list()){

}
#-------------------------------------------------------------------------------
#' @title input_data_process
#'
#' @description Cross-reference input data with serological, annual case/death and/or outbreak data for comparison
#'
#' @details This function, used to prepare input data for functions used to calculate the likelihood of observed data,
#'   amends a list of population and vaccination data used as input for other functions,
#'   cross-referencing it with seroprevalence, case and/or outbreak data and adding connection information.
#'
#' @param input_data List of population and vaccination data for multiple regions (created using create_input_data()
#'   function and usually loaded from an RDS file)
#' @param obs_sero_data Seroprevalence data for comparison, by region, year & age group, in format no. samples/no.
#'   positives
#' @param obs_case_data Annual reported case/death data for comparison, by region and year, in format no. cases/no.
#'   deaths
#' @param obs_outbreak_data Outbreak Y/N data for comparison, by region and year, in format 0 = no outbreaks,
#'   1 = 1 or more outbreak(s)
#'
#' @export
#'
input_data_process <- function(input_data=list(),obs_sero_data=NULL,obs_case_data=NULL,obs_outbreak_data=NULL){

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
