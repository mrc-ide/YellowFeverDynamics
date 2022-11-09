Full_Model_Run_Delay2 <- function(FOI_spillover=0.0,R0=1.0,vacc_data=list(),pop_data=list(),year0=1940,mode_start=0,
                           n_particles=1,n_threads=1,year_end=2000,year_data_begin=1999,vaccine_efficacy=1.0,
                           start_SEIRV=list(),dt=1.0) {

  assert_that(n_particles>0)
  assert_that(n_particles<=20)
  assert_that(n_threads<=n_particles)
  assert_that(n_threads>0)

  N_age=length(pop_data[1,]) #Number of age groups
  pars=parameter_setup(FOI_spillover,R0,vacc_data,pop_data,year0,mode_start,year_end,
                       year_data_begin,vaccine_efficacy,start_SEIRV,dt)
  nd1=as.integer((t_incubation+t_latent)/dt)
  nd2=as.integer(t_infectious/dt)
  nd=nd1+nd2
  pars[[length(pars)+1]]=nd1*N_age
  names(pars)[[length(pars)]]="n_delay_steps1"
  pars[[length(pars)+1]]=nd2*N_age
  names(pars)[[length(pars)]]="n_delay_steps2"

  x <- FullModelODDelay$new(pars,time = 1,n_particles = n_particles,n_threads = n_threads)

  n_nv=4 #Number of non-vector outputs
  n_data_pts=((6+nd)*N_age)+n_nv #Number of data values per time point in output
  n_steps=(year_end-year0)*(365/dt) #Total number of output time points
  step0=(year_data_begin-year0)*(365/dt) #Step at which data starts being saved for final output
  t_pts_out=n_steps-step0 #Number of time points in final output data
  x_res <- array(NA, dim = c(n_data_pts, n_particles, t_pts_out))
  for(t in step0:n_steps){
    x_res[,,t-step0] <- x$run(t)
  }

  # return(NULL)
  return(list(day=array(x_res[2,,],dim=c(n_particles,t_pts_out)),
              year=array(x_res[3,,],dim=c(n_particles,t_pts_out)),
              FOI_total=array(x_res[4,,],dim=c(n_particles,t_pts_out)),
              S=array(x_res[c((1+n_nv):(N_age+n_nv)),,],dim=c(N_age,n_particles,t_pts_out)),
              E=array(x_res[c((N_age+1+n_nv):((2*N_age)+n_nv)),,],dim=c(N_age,n_particles,t_pts_out)),
              E_delay=array(x_res[c(((2*N_age)+1+n_nv):(((2+nd1)*N_age)+n_nv)),,],dim=c(N_age,nd1,n_particles,t_pts_out)),
              I=array(x_res[c((((2+nd1)*N_age)+1+n_nv):(((3+nd1)*N_age)+n_nv)),,],dim=c(N_age,n_particles,t_pts_out)),
              I_delay=array(x_res[c((((2+nd1)*N_age)+1+n_nv):(((2+nd)*N_age)+n_nv)),,],dim=c(N_age,nd2,n_particles,t_pts_out)),
              R=array(x_res[c((((3+nd)*N_age)+1+n_nv):(((4+nd)*N_age)+n_nv)),,],dim=c(N_age,n_particles,t_pts_out)),
              V=array(x_res[c((((4+nd)*N_age)+1+n_nv):(((5+nd)*N_age)+n_nv)),,],dim=c(N_age,n_particles,t_pts_out)),
              C=array(x_res[c((((5+nd)*N_age)+1+n_nv):(((6+nd)*N_age)+n_nv)),,],dim=c(N_age,n_particles,t_pts_out))))
}

case_data_generate_delay2 <- function(FOI_spillover=0.0,R0=1.0,vacc_data=list(),pop_data=list(),year0=1940,
                               mode_start=0,n_reps=1,year_end=2000,year_data_begin=1999,
                               vaccine_efficacy=vaccine_efficacy,start_SEIRV=list(),dt=1.0) {

  assert_that(n_reps>0)

  division=10
  n_particles0=min(division,n_reps)
  n_threads=min(10,n_particles0)
  n_divs=ceiling(n_reps/division)
  if(n_divs==1){
    n_particles_list=n_particles0
  } else {
    n_particles_list=c(rep(n_particles0,n_divs-1),n_reps-(division*(n_divs-1)))
  }

  N_age=length(pop_data[1,]) #Number of age groups
  pars=parameter_setup(FOI_spillover,R0,vacc_data,pop_data,year0,mode_start,year_end,
                       year_data_begin,vaccine_efficacy,start_SEIRV,dt)
  nd1=as.integer((t_incubation+t_latent)/dt)
  nd2=as.integer(t_infectious/dt)
  nd=nd1+nd2
  pars[[length(pars)+1]]=nd1*N_age
  names(pars)[[length(pars)]]="n_delay_steps1"
  pars[[length(pars)+1]]=nd2*N_age
  names(pars)[[length(pars)]]="n_delay_steps2"

  n_nv=4 #Number of non-vector outputs
  n_data_pts=((6+nd)*N_age)+n_nv #Number of data values per time point in output
  n_steps=(year_end-year0)*(365/dt) #Total number of output time points
  step0=(year_data_begin-year0)*(365/dt) #Step at which data starts being saved for final output
  steps=n_steps-step0
  results_data=list(year=sort(rep(c(year_data_begin:(year_end-1)),(365/dt))),C=array(data=rep(0,n_reps*steps),
                                                                                     dim=c(n_reps,steps)))
  pts_select=c((((5+nd)*N_age)+1+n_nv):(((6+nd)*N_age)+n_nv))

  for(div in 1:n_divs){
    n_particles=n_particles_list[div]
    reps=c(1:n_particles)+((div-1)*division)
    x <- FullModelODDelay$new(pars=pars,time = 1,n_particles = n_particles,n_threads = n_threads)

    x_res <- array(NA, dim = c(n_data_pts, n_particles))
    for(t in step0:n_steps){
      x_res <- x$run(t)
      if(n_particles==1){
        results_data$C[reps,t-step0]=sum(x_res[pts_select])
      } else {
        results_data$C[reps,t-step0]=colSums(x_res[pts_select,],dims=1)
      }
    }
  }

  return(results_data)
}
