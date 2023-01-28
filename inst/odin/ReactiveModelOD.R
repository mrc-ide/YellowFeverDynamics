
dt <- user() #Time increment in days
initial(time) <- 0 #Initial value of time in days
update(time) <- (step + 1) * dt

#Parameters
FOI_spillover <- user() #Spillover force of infection
R0 <- user() #Basic reproduction number
N_age <- user() #Number of age categories
vacc_rate_annual[,,] <- user() #Daily rate of vaccination by age and year (non-emergency and emergency)
vaccine_efficacy <- user() #Proportion of vaccinations which successfully protect the recipient
p_rep[] <- user() #Proportion of infections reported (2 values depending on outbreak flag conditions)
outbreak_threshold1 <- user() #Threshold total no. reported cases to trigger outbreak flag 1
cluster_threshold1 <- user()  #Threshold current infectious fraction to trigger cluster flag 1

#initial conditions
year0 <- user()  #Starting year
Sus0[] <- user() #Susceptible population by age group at start
Exp0[] <- user() #Exposed population by age group at start
Inf0[] <- user() #Infectious population by age group at start
Rec0[] <- user() #Recovered population by age group at start
Vac0[] <- user() #Vaccinated population by age group at start
Cas0[] <- user() #Daily cases population by age group at start
dP1_all[,] <- user() #Daily increase in number of people by age group (people arriving in group due to age etc.)
dP2_all[,] <- user() #Daily decrease in number of people by age group (people leaving group due to age etc.)
n_years <- user() #Number of years for which model to be run

zero <- 0 #Minimum population setting to avoid negative numbers
one <- 1
Pmin <- 0 #Minimum population setting to avoid negative numbers
FOI_max <- 1.0 #Upper threshold for total force of infection to avoid more infections than people in a group
t_incubation <- 5 #Time for cases to incubate in mosquito
t_latent <- 5 #Latent period before cases become infectious
rate1 <- dt/(t_incubation+t_latent)
t_infectious <- 5 #Time cases remain infectious
rate2 <- dt/t_infectious
surv_delay <- 60 #Average time delay between symptom onset and case confirmation
rate3 <- dt/surv_delay
beta <- (R0*dt)/t_infectious #Daily exposure rate
FOI_sum <-  min(FOI_max,beta*(sum(I)/P_tot) + (FOI_spillover*dt)) #Total force of infection
year_i <- floor((step*dt)/365) + 1 #Number of years since start, as integer
dP1[1:N_age] <- dP1_all[i, as.integer(year_i)]*dt #Increase in population by age group over 1 time increment
dP2[1:N_age] <- dP2_all[i, as.integer(year_i)]*dt #Decrease in population by age group over 1 time increment

E_new[1:N_age] <- rbinom(as.integer(S[i]), FOI_sum) #New exposed individuals by age group
I_new[1:N_age] <- E[i]*rate1     #New infectious individuals by age group
R_new[1:N_age] <- I[i]*rate2     #New recovered individuals by age group
P_nV[1:N_age] <- S[i] + R[i] #Total vaccine-targetable population by age group
inv_P_nV[1:N_age] <- 1.0/P_nV[i]
P[1:N_age] <- S[i] + E[i] + I[i] + R[i] + V[i] #Total population by age group
P_tot <- sum(P) #Total overall population
inv_P[1:N_age] <- 1.0/P[i]
vacc_rate[1:N_age] <- vacc_rate_annual[i,as.integer(year_i),as.integer(flag3+1)]*vaccine_efficacy*dt*P[i] #Total no. vaccinations by age
VR_check1 <- 1
outbreak_flag1 <- min(one,max(zero,1+C_rep_total-outbreak_threshold1))
p_rep_cur <- p_rep[as.integer(flag3+1)]
C_rep_new[1:N_age] <- rbinom(as.integer(I_new[i]),p_rep_cur) #Daily new reported cases by age group
F_I_total <- sum(I)/P_tot #Total no. currently infectious people as fraction of population - check for cluster flag
cluster_flag1 <- as.integer(max(flag2a,min(one,F_I_total/cluster_threshold1)))
flag_emergency <- max(flag1b,flag2b)

#Updates to output values at each time increment
update(day) <- day + dt
update(year) <- year_i + year0 - 1
update(FOI_total) <- FOI_sum
update(C_rep_total) <- C_rep_total + sum(C_rep_new) #Running total reported cases across all ages
update(flag1a) <- outbreak_flag1 #0 = No cases so far; 1 = 1+ cases so far
update(flag1b) <- min(one,flag1b + (outbreak_flag1*rate3)) #flag1 with delay (converted to integer on use)
update(flag2a) <- cluster_flag1 #0 = No cluster; 1 = cluster (high enough no. infectious people)
update(flag2b) <- min(one,flag2b + (cluster_flag1*rate3)) #flag1 with delay (converted to integer on use)
update(flag3) <- flag_emergency #0 = No emergency (flags 1+2 not tripped); 1 = emergency (flag 1 and/or 2 tripped)
update(report_rate) <- p_rep_cur
update(VR_check) <- VR_check1
update(S[1]) <- max(Pmin,S[1] - E_new[1] - vacc_rate[1]*S[1]*inv_P_nV[1] + dP1[1] - (dP2[1]*S[1]*inv_P[1]))
update(S[2:N_age]) <- max(Pmin,S[i] - E_new[i] - vacc_rate[i]*S[i]*inv_P_nV[i] + (dP1[i]*S[i-1]*inv_P[i-1]) - (dP2[i]*S[i]*inv_P[i]))
update(E[1:N_age]) <- max(Pmin,E[i] + E_new[i] - I_new[i])
update(I[1:N_age]) <- max(Pmin,I[i] + I_new[i] - R_new[i])
update(R[1]) <- max(Pmin,R[1] + R_new[1] - vacc_rate[1]*R[1]*inv_P_nV[1] - (dP2[1]*R[1]*inv_P[1]))
update(R[2:N_age]) <- max(Pmin,R[i] + R_new[i] - vacc_rate[i]*R[i]*inv_P_nV[i] + (dP1[i]*R[i-1]*inv_P[i-1]) - (dP2[i]*R[i]*inv_P[i]))
update(V[1]) <- max(Pmin,V[1] + vacc_rate[1] - (dP2[1]*V[1]*inv_P[1]))
update(V[2:N_age]) <- max(Pmin,V[i] + vacc_rate[i] + (dP1[i]*V[i-1]*inv_P[i-1]) - (dP2[i]*V[i]*inv_P[i]))
update(C[1:N_age]) <- I_new[i]
update(C_rep[1:N_age]) <- C_rep_new[i]

#Initial values
initial(day) <- 0
initial(year) <- year0-1
initial(FOI_total) <- FOI_spillover
initial(C_rep_total) <- 0
initial(flag1a) <- 0
initial(flag1b) <- 0
initial(flag2a) <- 0
initial(flag2b) <- 0
initial(flag3) <- 0
initial(report_rate) <- p_rep[1]
initial(VR_check) <- 0
initial(S[1:N_age]) <- Sus0[i]
initial(E[1:N_age]) <- Exp0[i]
initial(I[1:N_age]) <- Inf0[i]
initial(R[1:N_age]) <- Rec0[i]
initial(V[1:N_age]) <- Vac0[i]
initial(C[1:N_age]) <- Cas0[i]
initial(C_rep[1:N_age]) <- 0

#Dimensions
dim(S) <- N_age
dim(E) <- N_age
dim(I) <- N_age
dim(R) <- N_age
dim(V) <- N_age
dim(C) <- N_age
dim(C_rep) <- N_age

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
dim(C_rep_new) <- N_age

dim(Sus0) <- N_age
dim(Exp0) <- N_age
dim(Inf0) <- N_age
dim(Rec0) <- N_age
dim(Vac0) <- N_age
dim(Cas0) <- N_age
dim(dP1_all) <- c(N_age, n_years)
dim(dP2_all) <- c(N_age, n_years)
dim(vacc_rate_annual) <- c(N_age, n_years, 2)
dim(p_rep) <- 2
