#-------------------------------------------------------------------------------
#' @title Delay_Model_Run_Many_Reps
#'
#' @description Run delay version of SEIRV model for single region for large number of repetitions
#'
#' @details Accepts epidemiological + population parameters and model settings; runs SEIRV model
#' for one region over a specified time period for a number of repetitions and outputs time-dependent SEIRV
#' values, infection numbers and/or total force of infection values. Variation of Model_Run_Delay() used for
#' running a large number of repetitions (>20).
#'
#' @param FOI_spillover Force of infection due to spillover from sylvatic reservoir
#' @param R0 Basic reproduction number for urban spread of infection
#' @param vacc_data Vaccination coverage in each age group by year
#' @param pop_data Population in each age group by year
#' @param years_data Incremental vector of years denoting years for which to save data
#' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' @param output_type Type of data to output:
#'   "full" = SEIRVC + FOI for all steps and ages
#'   "case" = annual total new infections (C) summed across all ages
#'   "sero" = annual SEIRV summed across all ages
#'   "case+sero" = annual SEIRVC summed across all ages
#'   "case_alt" = annual total new infections not combined by age
#'   "case_alt2" = total new infections combined by age for all steps
#' @param year0 First year in population/vaccination data
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRV input in list from previous run(s)
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param dt Time increment in days to use in model (should be 1.0, 2.5 or 5.0 days)
#' @param n_reps Number of repetitions (used to set number of particles and threads)
#' @param division Number of particles to run in one go (up to 20)
#' '
#' @export
#'
Delay_Model_Run_Many_Reps <- function(FOI_spillover = 0.0,R0 = 1.0,vacc_data = list(),pop_data = list(),years_data = c(1940:1941),
                                      start_SEIRV = list(), output_type = "full", year0 = 1940, mode_start = 0,
                                      vaccine_efficacy = 1.0, dt = 1.0, n_reps=1, division=10) {
  
  assert_that(division<=20)
  n_particles0=min(division,n_reps)
  n_threads=min(division,n_particles0)
  n_divs=ceiling(n_reps/division)
  if(n_divs==1){
    n_particles_list=n_particles0
  } else {
    n_particles_list=c(rep(n_particles0,n_divs-1),n_reps-(division*(n_divs-1)))
  }
  
  n_nv=3 #Number of non-vector outputs
  N_age=length(pop_data[1,]) #Number of age groups
  nd1 <- (t_incubation+t_latent)/dt
  nd2 <- t_infectious/dt
  nd=nd1+nd2
  n_data_pts=((6+nd1+nd2)*N_age)+n_nv #Number of data values per time point in output
  step_begin=((years_data[1]-year0)*(365/dt)) #Step at which data starts being saved for final output
  step_end=((max(years_data)+1-year0)*(365/dt))-1 #Step at which to end
  t_pts_out=step_end-step_begin+1 #Number of time points in final output data
  
  if(output_type=="full"){
    dimensions=c(N_age,n_reps,t_pts_out)
    output_data=list(day=rep(NA,t_pts_out),year=rep(NA,t_pts_out),FOI_total=array(NA,c(n_reps,t_pts_out)),
                     S=array(NA,dim=dimensions),E=array(NA,dim=dimensions),I=array(NA,dim=dimensions),
                     R=array(NA,dim=dimensions),V=array(NA,dim=dimensions),C=array(NA,dim=dimensions))
  } else {
    if(output_type=="case_alt2"){
      output_data=list(day=rep(NA,t_pts_out),year=rep(NA,t_pts_out))
      output_data$C=array(0,dim=c(n_reps,t_pts_out))
    } else {
      n_years=length(years_data)
      output_data=list(year=years_data)
      if(output_type=="case+sero" || output_type=="sero"){
        output_data$V=output_data$R=output_data$I=output_data$E=output_data$S=array(0,dim=c(N_age,n_particles,n_years))
      }
      if(output_type=="case+sero" || output_type=="case"){
        output_data$C=array(0,dim=c(n_reps,n_years))
      }
      if(output_type=="case_alt"){
        output_data$C=array(0,dim=c(N_age,n_reps,n_years))
      }
    }
  }
  
  for(div in 1:n_divs){
    n_particles=n_particles_list[div]
    if(div==1){n_p0=0}else{n_p0=sum(n_particles_list[c(1:(div-1))])}
    
    x <- SEIRVModelDelay$new(pars=parameter_setup(FOI_spillover,R0,vacc_data,pop_data,year0,years_data,mode_start,
                                                  vaccine_efficacy,start_SEIRV,dt),
                             time = 0, n_particles = n_particles, n_threads = n_threads, deterministic = FALSE)
    
    x_res <- array(NA, dim = c(n_data_pts, n_particles, t_pts_out))
    for(step in step_begin:step_end){
      x_res[,,step-step_begin+1] <- x$run(step)
    }
    if(step_begin==0){x_res[2,,1]=rep(year0,n_particles)}
    
    if(output_type=="full"){
      n_p_values=c(1:n_particles)+n_p0
      dimensions=c(N_age,n_particles,t_pts_out)
      output_data$day=x_res[1,1,]
      output_data$year=x_res[2,1,]
      output_data$FOI_total[n_p_values,]=array(x_res[3,,]/dt,dim=c(n_particles,t_pts_out))
      output_data$S[,n_p_values,]=array(x_res[c((1+n_nv):(N_age+n_nv)),,],dim=dimensions)
      output_data$E[,n_p_values,]=array(x_res[c((N_age+1+n_nv):((2*N_age)+n_nv)),,],dim=dimensions)
      output_data$I[,n_p_values,]=array(x_res[c((((2+nd1)*N_age)+1+n_nv):(((3+nd1)*N_age)+n_nv)),,],dim=dimensions)
      output_data$R[,n_p_values,]=array(x_res[c((((3+nd)*N_age)+1+n_nv):(((4+nd)*N_age)+n_nv)),,],dim=dimensions)
      output_data$V[,n_p_values,]=array(x_res[c((((4+nd)*N_age)+1+n_nv):(((5+nd)*N_age)+n_nv)),,],dim=dimensions)
      output_data$C[,n_p_values,]=array(x_res[c((((5+nd)*N_age)+1+n_nv):(((6+nd)*N_age)+n_nv)),,],dim=dimensions)
    } else {
      if(output_type=="case+sero" || output_type=="sero"){
        for(n_year in 1:n_years){
          pts=c(1:t_pts_out)[x_res[2,1,]==years_data[n_year]]
          for(n_p in 1:n_particles){
            n_p2=n_p+n_p0
            output_data$S[,n_p2,n_year]=rowMeans(x_res[c((1+n_nv):(N_age+n_nv)),n_p,pts])
            output_data$E[,n_p2,n_year]=rowMeans(x_res[c((N_age+1+n_nv):((2*N_age)+n_nv)),n_p,pts])
            output_data$I[,n_p2,n_year]=rowMeans(x_res[c((((2+nd1)*N_age)+1+n_nv):(((3+nd1)*N_age)+n_nv)),n_p,pts])
            output_data$R[,n_p2,n_year]=rowMeans(x_res[c((((3+nd)*N_age)+1+n_nv):(((4+nd)*N_age)+n_nv)),n_p,pts])
            output_data$V[,n_p2,n_year]=rowMeans(x_res[c((((4+nd)*N_age)+1+n_nv):(((5+nd)*N_age)+n_nv)),n_p,pts])
          }
        }
      }
      if(output_type=="case+sero" || output_type=="case"){
        for(n_year in 1:n_years){
          pts=c(1:t_pts_out)[x_res[2,1,]==years_data[n_year]]
          for(n_p in 1:n_particles){
            n_p2=n_p+n_p0
            output_data$C[n_p2,n_year]=sum(x_res[c((((5+nd)*N_age)+1+n_nv):(((6+nd)*N_age)+n_nv)),n_p,pts])
          }
        }
      }
      if(output_type=="case_alt"){
        for(n_year in 1:n_years){
          pts=c(1:t_pts_out)[x_res[2,1,]==years_data[n_year]]
          for(n_p in 1:n_particles){
            n_p2=n_p+n_p0
            output_data$C[,n_p2,n_year]=rowSums(x_res[c((((5+nd)*N_age)+1+n_nv):(((6+nd)*N_age)+n_nv)),n_p,pts])
          }
        }
      }
      if(output_type=="case_alt2"){
        if(n_p0==0){
          output_data$day=x_res[1,1,]
          output_data$year=x_res[2,1,]
        }
        for(pt in 1:t_pts_out){
          for(n_p in 1:n_particles){
            n_p2=n_p+n_p0
            output_data$C[n_p2,pt]=sum(x_res[c((((5+nd)*N_age)+1+n_nv):(((6+nd)*N_age)+n_nv)),n_p,pt])
          }
        }
      }
    }
    x_res<-NULL
    gc()
  }
  
  return(output_data)
}