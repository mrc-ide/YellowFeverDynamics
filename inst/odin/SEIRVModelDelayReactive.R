# SEIRV model with reactive functionality using delay approach
# TODO - adjust for time varying FOI/R0




time_inc <- parameter() #Time increment in days
initial(day) <- time_inc #Initial value of time in days
update(day) <- day + time_inc

#Parameters---------------------------------------------------------------------
#t_incubation <- parameter() #Length in days of yellow fever incubation period in mosquito vectors
#t_latent <- parameter() #Length in days of latent period in humans exposed to yellow dever
t_infectious <- parameter() #Length of infectious period in humans with yellow fever
FOI_spillover <- parameter() #Spillover force of infection (per day) at each time point
R0 <- parameter() #Basic reproduction number for human-human transmission at each time point
N_age <- parameter() #Number of age categories
vacc_rate_daily <- parameter() #Daily rate of vaccination by age and year (non-emergency and emergency)
vaccine_efficacy <- parameter() #Proportion of vaccinations which successfully protect the recipient
response_delay <- parameter() #Delay time in days between a flag being triggered and emergency conditions coming into effect
p_rep <- parameter() #Proportion of infections reported (2 values depending on outbreak flag conditions)
case_threshold <- parameter() #Threshold total no. reported cases to trigger outbreak flag 1
cluster_threshold <- parameter()  #Threshold current infectious fraction to trigger cluster flag 1
vacc_cov_cam <- parameter() #TBA
t_cam <- parameter() #TBA

#Initial conditions-------------------------------------------------------------
year0 <- parameter()  #Starting year
S_0 <- parameter() #Susceptible population by age group at start
E_0 <- parameter() #Exposed population by age group at start
E_delay0 <- parameter() #TBA
I_0 <- parameter() #Infectious population by age group at start
I_delay0 <- parameter() #TBA
R_0 <- parameter() #Recovered population by age group at start
V_0 <- parameter() #Vaccinated population by age group at start
dP1_all <- parameter() #Daily increase in number of people by age group (people arriving in group due to age etc.)
dP2_all <- parameter() #Daily decrease in number of people by age group (people leaving group due to age etc.)
n_years <- parameter() #Number of years for which model to be run
n_t_pts <- parameter() #Total number of time points
Pmin <- 1.0e-99 #Minimum population setting to avoid negative numbers
FOI_max <- 1.0 #Upper threshold for total force of infection to avoid more infections than people in a group
np_E_delay <- parameter()
np_I_delay <- parameter()
di1 <- np_E_delay-N_age
di2 <- np_I_delay-N_age
one <- 1
rate3 <- time_inc/response_delay
rate4 <- time_inc/t_cam
t_pt <- day/time_inc #Number of time points passed
beta <- (R0[t_pt]*time_inc)/t_infectious #Daily exposure rate
FOI_sum <-  min(FOI_max,beta*(sum(I)/P_tot) + (FOI_spillover[t_pt]*time_inc)) #Total force of infection

year_i <- floor(day/365)+1 #Number of years since start, as integer

dP1[1:N_age] <- dP1_all[i, year_i]*time_inc #Increase in population by age group over 1 time increment
dP2[1:N_age] <- dP2_all[i, year_i]*time_inc #Decrease in population by age group over 1 time increment

E_new[1:N_age] <- Binomial(as.integer(S[i]), FOI_sum) #New exposed individuals by age group

I_new[1:N_age] <- E_delay[as.integer(i+di1)]     #New infectious individuals by age group

R_new[1:N_age] <- I_delay[as.integer(i+di2)]     #New recovered individuals by age group

P_nV[1:N_age] <- S[i] + R[i] #Total vaccine-targetable population by age group
inv_P_nV[1:N_age] <- 1.0/P_nV[i]
P[1:N_age] <- P_nV[i] + V[i] #Total population by age group (excluding E+I)
P_tot <- sum(P) #Total overall population (excluding E+I)
inv_P[1:N_age] <- 1.0/P[i]
vacc_rate_cam[1:N_age] <- if(flag3==0) 0 else (vacc_cov_cam[i]*(1.0 - (V[i]*inv_P[i])))/t_cam
vacc_rate[1:N_age] <- (vacc_rate_daily[i,year_i] + (if(flag3==0) 0 else vacc_rate_cam[i]*ceiling(1-flag4)))*vaccine_efficacy*time_inc*P[i]
case_flag <- if(C_rep_total >= case_threshold) 1 else 0
p_rep_cur <- if(flag3==1) p_rep[2] else p_rep[1]
C_rep_new <- Binomial(as.integer(sum(I_new)),p_rep_cur) #Daily new reported cases across all ages
F_I_total <- sum(I)/P_tot #Total no. currently infectious people as fraction of population - check for cluster flag
cluster_flag <- if(F_I_total>=cluster_threshold) 1 else 0

#Updates to output values at each time increment--------------------------------
update(year) <- year_i + year0 - 1
update(FOI_total) <- FOI_sum
update(C_rep_total) <- C_rep_total + C_rep_new #Running total reported cases across all ages
update(flag1) <- min(one,flag1 + (case_flag*rate3)) #case flag with delay (converted to integer on use)
update(flag2) <- min(one,flag2 + (cluster_flag*rate3)) #cluster flag with delay (converted to integer on use)
update(flag3) <- if(flag1==1) 1 else (if(flag2==1) 1 else 0) #0 = No emergency (flags 1+2 not tripped); 1 = emergency (flag 1 and/or 2 tripped)
update(flag4) <- if(flag3==0) 0 else min(one,flag4+rate4)
update(report_rate) <- p_rep_cur
update(S[1]) <- max(Pmin,S[1] - E_new[1] - vacc_rate[1]*S[1]*inv_P_nV[1] + dP1[1] - (dP2[1]*S[1]*inv_P[1]))
update(S[2:N_age]) <- max(Pmin,S[i] - E_new[i] - vacc_rate[i]*S[i]*inv_P_nV[i] + (dP1[i]*S[i-1]*inv_P[i-1]) - (dP2[i]*S[i]*inv_P[i]))
update(E[1:N_age]) <- max(Pmin,E[i] + E_new[i] - I_new[i])
update(E_delay[(N_age+1):np_E_delay]) <- E_delay[i-N_age]
update(E_delay[1:N_age]) <- E_new[i]
update(I[1:N_age]) <- max(Pmin,I[i] + I_new[i] - R_new[i])
update(I_delay[(N_age+1):np_I_delay]) <- I_delay[i-N_age]
update(I_delay[1:N_age]) <- I_new[i]
update(R[1]) <- max(Pmin,R[1] + R_new[1] - vacc_rate[1]*R[1]*inv_P_nV[1] - (dP2[1]*R[1]*inv_P[1]))
update(R[2:N_age]) <- max(Pmin,R[i] + R_new[i] - vacc_rate[i]*R[i]*inv_P_nV[i] + (dP1[i]*R[i-1]*inv_P[i-1]) - (dP2[i]*R[i]*inv_P[i]))
update(V[1]) <- max(Pmin,V[1] + vacc_rate[1] - (dP2[1]*V[1]*inv_P[1]))
update(V[2:N_age]) <- max(Pmin,V[i] + vacc_rate[i] + (dP1[i]*V[i-1]*inv_P[i-1]) - (dP2[i]*V[i]*inv_P[i]))
update(C[1:N_age]) <- I_new[i]
#update(C_rep[1:N_age]) <- C_rep_new[i]

#Initial values-----------------------------------------------------------------
initial(year) <- year0
initial(FOI_total) <- FOI_spillover[1]
initial(C_rep_total) <- 0
initial(flag1) <- 0
initial(flag2) <- 0
initial(flag3) <- 0
initial(flag4) <- 0
initial(report_rate) <- p_rep[1]
initial(S[1:N_age]) <- S_0[i]
initial(E[1:N_age]) <- E_0[i]
initial(E_delay[1:np_E_delay]) <- E_delay0[i]
initial(I[1:N_age]) <- I_0[i]
initial(I_delay[1:np_I_delay]) <- I_delay0[i]
initial(R[1:N_age]) <- R_0[i]
initial(V[1:N_age]) <- V_0[i]
initial(C[1:N_age]) <- 0
#initial(C_rep[1:N_age]) <- 0

#Dimensions---------------------------------------------------------------------
dim(S) <- N_age
dim(E) <- N_age
dim(E_delay) <- np_E_delay
dim(I) <- N_age
dim(I_delay) <- np_I_delay
dim(R) <- N_age
dim(V) <- N_age
dim(C) <- N_age
#dim(C_rep) <- N_age
dim(dP1)<-N_age
dim(dP2)<-N_age
dim(E_new) <- N_age
dim(I_new) <- N_age
dim(R_new) <- N_age




dim(P_nV) <- N_age
dim(inv_P_nV) <- N_age
dim(P) <- N_age
dim(inv_P) <- N_age
dim(vacc_rate) <- N_age
#dim(C_rep_new) <- N_age
dim(S_0) <- N_age
dim(E_0) <- N_age
dim(E_delay0) <- np_E_delay
dim(I_0) <- N_age
dim(I_delay0) <- np_I_delay
dim(R_0) <- N_age
dim(V_0) <- N_age
dim(dP1_all) <- c(N_age, n_years)
dim(dP2_all) <- c(N_age, n_years)
dim(vacc_rate_daily) <- c(N_age, n_years)
dim(FOI_spillover) <- n_t_pts
dim(R0) <- n_t_pts
dim(p_rep) <- 2
dim(vacc_cov_cam) <- N_age
dim(vacc_rate_cam) <- N_age



