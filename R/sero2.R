# R file for functions relating to serological data in YellowFeverDynamics package
#-------------------------------------------------------------------------------
#' @title sero_calculate2
#'
#' @description Calculate number of "samples" and number of "positives" from modelled data for specified age range(s)
#' and year(s)
#'
#' @details Takes in information on minimum and maximum ages of desired range(s), year(s) for which to calculate
#' number of "samples" (people eligible for testing) and "positives" (people who would test positive), plus vc_factor
#' (proportion of people for whom vaccination status unknown)
#'
#' @param sero_data = Data frame containing years, minimum and maximum ages, and values of vc_factor (proportion of
#' people for whom vaccination status unknown)
#' @param model_data = Output of Basic_Model_Run or Full_Model_Run
#' '
#' @export
#'
sero_calculate2 <- function(sero_data=list(),model_data=list()){

  #TODO - Add assert_that functions
  nrows=nrow(sero_data)
  output_frame=data.frame(samples=rep(NA,nrows),positives=rep(NA,nrows))

  for(i in 1:nrows){
    ages=c((sero_data$age_min[i]+1):sero_data$age_max[i])
    year=sero_data$year[i]
    vc_factor=sero_data$vc_factor[i]
    days=which(model_data$year==year)
    S_sum=sum(model_data$S[days,ages])
    E_sum=sum(model_data$E[days,ages])
    I_sum=sum(model_data$I[days,ages])
    R_sum=sum(model_data$R[days,ages])
    samples=S_sum+E_sum+I_sum+R_sum
    if(vc_factor>1){
      V_sum=sum(model_data$V[days,ages])
      if(vc_factor==1){
        samples=samples+V_sum
        positives=R_sum+V_sum
      } else {
        positives=((1.0-vc_factor)*R_sum)+(vc_factor*(samples/T_sum)*(R_sum+V_sum))
      }
    } else {
      positives=R_sum
    }
    output_frame$samples[i]=samples
    output_frame$positives[i]=positives
  }

  return(output_frame)
}
