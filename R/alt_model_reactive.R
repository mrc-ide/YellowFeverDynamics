#-------------------------------------------------------------------------------
#' @title Reactive_Model_Run
#'
#' @description Run reactive version of SEIRV model
#'
#' @details Accepts epidemiological + population parameters and model settings; runs full version of SEIRV model
#' for one region over a specified time period for a number of particles/threads and outputs time-dependent SEIRV
#' values, infection numbers and total force of infection values. Alternate version incorporating case reporting
#' and reactive surveillance/control measures based on case numbers
#'
#' @param FOI_spillover Force of infection due to spillover from sylvatic reservoir
#' @param R0 Basic reproduction number for urban spread of infection
#' @param vacc_data1 Vaccination coverage in each age group by year (non-emergency)
#' @param vacc_data2 Vaccination coverage in each age group by year (emergency)
#' @param pop_data Population in each age group by year
#' @param year0 First year in population/vaccination data
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRVC input in list from previous run(s)
#' @param n_particles number of particles to use
#' @param n_threads number of threads to use
#' @param year_end year to run up to
#' @param year_data_begin year to begin saving data
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' @param dt Time increment in days to use in model (should be either 1.0 or 5.0 days)
#' @param p_rep Probabilities of an infection being reported as a case under different conditions (TBA)
#' @param outbreak_threshold1 Threshold total no. reported cases to trigger outbreak flag 1
#' @param cluster_threshold1 Threshold current infectious fraction to trigger cluster flag 1
#'
#' @export
#'
Reactive_Model_Run <- function(FOI_spillover=0.0,R0=1.0,vacc_data1=list(),vacc_data2=list(),pop_data=list(),year0=1940,
                           mode_start=0,n_particles=1,n_threads=1,year_end=2000,year_data_begin=1999,
                           vaccine_efficacy=1.0,start_SEIRV=list(),dt=1.0,p_rep=c(1.0e-6,1.0e-6),outbreak_threshold1=1,
                           cluster_threshold1=1.0) {

  assert_that(n_particles %in% c(1:20),msg="Number of particles must be an integer between 1 and 20")
  assert_that(n_threads<=n_particles,msg="Number of threads must be equal to or less than number of particles")
  assert_that(n_threads>0,msg="Number of threads must be between 1 and number of particles")

  x <- ReactiveModelOD$new(pars=parameter_setup_react(FOI_spillover,R0,vacc_data1,vacc_data2,pop_data,year0,mode_start,
                                            year_end,year_data_begin,vaccine_efficacy,start_SEIRV,dt,p_rep,
                                            outbreak_threshold1,cluster_threshold1),
                       time = 1,n_particles = n_particles,n_threads = n_threads)

  n_nv=12 #Number of non-vector outputs at beginning of output
  N_age=length(pop_data[1,]) #Number of age groups
  t_pts_all=c(1:((year_end-year0)*(365/dt))) #All output time points
  n_data_pts=(7*N_age)+n_nv #Number of data values per time point in output
  n_steps=length(t_pts_all) #Total number of output time points
  step0=(year_data_begin-year0)*(365/dt) #Step at which data starts being saved for final output
  t_pts_out=n_steps-step0 #Number of time points in final output data
  x_res <- array(NA, dim = c(n_data_pts, n_particles, t_pts_out))
  for(t in step0:n_steps){
    x_res[,,t-step0] <- x$run(t)
  }

  return(list(day=array(x_res[2,,],dim=c(n_particles,t_pts_out)),
          year=array(x_res[3,,],dim=c(n_particles,t_pts_out)),
          FOI_total=array(x_res[4,,],dim=c(n_particles,t_pts_out)),
          C_rep_total=array(x_res[5,,],dim=c(n_particles,t_pts_out)),
          flag1a=array(x_res[6,,],dim=c(n_particles,t_pts_out)),flag1b=array(x_res[7,,],dim=c(n_particles,t_pts_out)),
          flag2a=array(x_res[8,,],dim=c(n_particles,t_pts_out)),flag2b=array(x_res[9,,],dim=c(n_particles,t_pts_out)),
          flag3=array(x_res[10,,],dim=c(n_particles,t_pts_out)),
          report_rate=array(x_res[11,,],dim=c(n_particles,t_pts_out)),
          VR_check=array(x_res[12,,],dim=c(n_particles,t_pts_out)),
              S=array(x_res[c((1+n_nv):(N_age+n_nv)),,],dim=c(N_age,n_particles,t_pts_out)),
              E=array(x_res[c((N_age+1+n_nv):((2*N_age)+n_nv)),,],dim=c(N_age,n_particles,t_pts_out)),
              I=array(x_res[c(((2*N_age)+1+n_nv):((3*N_age)+n_nv)),,],dim=c(N_age,n_particles,t_pts_out)),
              R=array(x_res[c(((3*N_age)+1+n_nv):((4*N_age)+n_nv)),,],dim=c(N_age,n_particles,t_pts_out)),
              V=array(x_res[c(((4*N_age)+1+n_nv):((5*N_age)+n_nv)),,],dim=c(N_age,n_particles,t_pts_out)),
              C=array(x_res[c(((5*N_age)+1+n_nv):((6*N_age)+n_nv)),,],dim=c(N_age,n_particles,t_pts_out)),
              C_rep=array(x_res[c(((6*N_age)+1+n_nv):((7*N_age)+n_nv)),,],dim=c(N_age,n_particles,t_pts_out))))
}
#-------------------------------------------------------------------------------
#' @title Parameter setup (reactive)
#'
#' @description Set up parameters to input into model
#'
#' @details Takes in multiple inputs, outputs list for use by odin.dust SEIRV model versions.
#'
#' @param FOI_spillover Force of infection due to spillover from sylvatic reservoir
#' @param R0 Basic reproduction number for urban spread of infection
#' @param vacc_data1 Vaccination coverage in each age group by year (non-emergency)
#' @param vacc_data2 Vaccination coverage in each age group by year (emergency)
#' @param pop_data Population in each age group by year
#' @param year0 First year in population/vaccination data
#' @param mode_start Flag indicating how to set initial population immunity level in addition to vaccination
#'  If mode_start=0, only vaccinated individuals
#'  If mode_start=1, shift some non-vaccinated individuals into recovered to give herd immunity
#'  If mode_start=2, use SEIRVC input in list from previous run(s)
#' @param year_end year to run up to
#' @param year_data_begin year to begin saving data
#' @param vaccine_efficacy Proportional vaccine efficacy
#' @param start_SEIRV SEIRV data from end of a previous run to use as input
#' @param dt Time increment in days to use in model (should be either 1.0 or 5.0 days)
#' @param p_rep Probabilities of an infection being reported as a case under different conditions (TBA)
#' @param outbreak_threshold1 Threshold total no. reported cases to trigger outbreak flag 1
#' @param cluster_threshold1 Threshold current infectious fraction to trigger cluster flag 1
#' '
#' @export
#'
parameter_setup_react <- function(FOI_spillover=0.0,R0=1.0,vacc_data1=list(),vacc_data2=list(),pop_data=list(),
                                  year0=1940,mode_start=0,year_end=2000,year_data_begin=1999,vaccine_efficacy=1.0,
                                  start_SEIRV=list(),dt=1.0,p_rep=c(1.0e-6,1.0e-6),outbreak_threshold1=1,
                                  cluster_threshold1=1.0){

  assert_that(length(pop_data[,1])>1) #TODO - msg
  assert_that(length(pop_data[1,])>1) #TODO - msg
  n_years=length(pop_data[,1])-1
  N_age=length(pop_data[1,])
  assert_that(length(vacc_data1[,1])==n_years+1,msg="Population and vaccination data 1 must be for same time periods")
  assert_that(length(vacc_data1[1,])==N_age,msg="Number of age groups in population and vaccination data 1 must match")
  assert_that(length(vacc_data2[,1])==n_years+1,msg="Population and vaccination data 2 must be for same time periods")
  assert_that(length(vacc_data2[1,])==N_age,msg="Number of age groups in population and vaccination data 2 must match")
  assert_that(mode_start %in% c(0,1,2),msg="mode_start must have value 0, 1 or 2")
  assert_that(vaccine_efficacy<=1.0 && vaccine_efficacy>=0.0,msg="Vaccine efficacy must be between 0 and 1")
  if(mode_start==2){assert_that(is.null(start_SEIRV$S)==FALSE,msg="When mode_start=2, start_SEIRV data is required")}
  assert_that(year_data_begin>=year0,msg="year_data_begin must be greater than or equal to year0")
  assert_that(year_data_begin<year_end,msg="year_data_begin must be less than year_end")
  assert_that(year_end-year0<=n_years,msg="Period year0->year_end must lie within population data")
  vacc_initial=vacc_data1[1,]
  assert_that(dt %in% c(1,2.5,5),msg="dt must have value 1, 2.5 or 5 days (must have integer number of points/year")
  assert_that(length(p_rep)==2,msg="2 reporting probability values required")
  inv_365=1.0/365.0

  P0=Cas0=Sus0=Exp0=Inf0=Rec0=Vac0=rep(0,N_age)
  dP1_all=dP2_all=array(rep(0,N_age*n_years),dim=c(N_age,n_years))
  vacc_rates=array(rep(0,N_age*n_years*2),dim=c(N_age,n_years,2))
  for(i in 1:N_age){
    P0[i]=max(1.0,pop_data[1,i]) #Set all population values to a minimum of 1 to avoid NaN values appearing
  }
  for(n_year in 1:n_years){
    for(i in 1:N_age){
      dP1_all[i,n_year]=max(1.0,pop_data[n_year+1,i])*inv_365
      dP2_all[i,n_year]=max(1.0,pop_data[n_year,i])*inv_365
      if(i==1){
        vacc_rates[i,n_year,1]=vacc_data1[n_year+1,i]*inv_365
        vacc_rates[i,n_year,2]=vacc_data2[n_year+1,i]*inv_365
      } else {
        vacc_rates[i,n_year,1]=max(0.0,vacc_data1[n_year+1,i]-vacc_data1[n_year,i-1])*inv_365
        vacc_rates[i,n_year,2]=max(0.0,vacc_data2[n_year+1,i]-vacc_data2[n_year,i-1])*inv_365
      }
    }
  }

  if(mode_start==0){
    Sus0=P0*(1.0-vacc_initial)
  }
  if(mode_start==1)
  {
    if(R0>1.0){
      herd_immunity=1.0-(1.0/R0)
    } else {
      herd_immunity=0.0
    }
    for(i in 1:N_age){
      if(vacc_initial[i]<herd_immunity){
        Rec0[i]=P0[i]*(herd_immunity-vacc_initial[i])
        Sus0[i]=P0[i]*(1.0-herd_immunity)
      } else {
        Sus0[i]=P0[i]*(1.0-vacc_initial[i])
      }
    }
  }
  if(mode_start==2){
    Sus0=start_SEIRV$S
    Exp0=start_SEIRV$E
    Inf0=start_SEIRV$I
    Rec0=start_SEIRV$R
    Vac0=start_SEIRV$V
    Cas0=rep(0,N_age)
  } else {
    assert_that(length(vacc_initial)==N_age)
    Vac0=P0*vacc_initial
  }

  return(list(FOI_spillover=FOI_spillover,R0=R0,vacc_rate_annual=vacc_rates,
              Cas0=Cas0,Exp0=Exp0,Inf0=Inf0,N_age=N_age,Rec0=Rec0,Sus0=Sus0,Vac0=Vac0,dP1_all=dP1_all,dP2_all=dP2_all,
              n_years=n_years,year0=year0,vaccine_efficacy=vaccine_efficacy,dt=dt,p_rep=p_rep,
              outbreak_threshold1=outbreak_threshold1,cluster_threshold1=cluster_threshold1))
}
