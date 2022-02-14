
dt <- user()
initial(time) <- 0
update(time) <- (step + 1) * dt

#Parameters
FOI_spillover <- user() #Spillover force of infection
R0 <- user() #Basic reproduction number
N_age <- user() #Number of age categories
vacc_rate_annual[,] <- user() #Daily rate of vaccination by age and year
vaccine_efficacy <- user()

#initial conditions
year0 <- user()
Sus0[] <- user() #Susceptible
Exp0[] <- user() #Exposed
Inf0[] <- user() #Infectious
Rec0[] <- user() #Recovered
Vac0[] <- user() #Vaccinated
Cas0[] <- user() #Daily cases
dP1_all[,] <- user() #Daily increase in number of people by age bracket due to births/ageing in from previous bracket
dP2_all[,] <- user() #Daily decrease in number of people by age bracket due to deaths/ageing into next bracket
n_years <- user()

Pmin <- 0
FOI_max <- 1.0
t_incubation <- 5 #Time for cases to incubate in mosquito (TBA)
t_latent <- 5 #Latent period before cases become infectious (TBA)
rate1=dt/(t_incubation+t_latent)
t_infectious <- 5 #Time cases remain infectious
rate2 <- dt/t_infectious
beta <- (R0*dt)/t_infectious #Daily exposure rate
FOI_sum <-  min(FOI_max,beta*(sum(I)/P_tot) + (FOI_spillover*dt)) #Total force of infection
year_i=floor((step*dt)/365) + 1
dP1[1:N_age] <- dP1_all[i, as.integer(year_i)]*dt
dP2[1:N_age] <- dP2_all[i, as.integer(year_i)]*dt

E_new[1:N_age] <- rbinom(as.integer(S[i]), FOI_sum) #New exposed individuals
I_new[1:N_age] <- E[i]*rate1     #New infectious individuals 
R_new[1:N_age] <- I[i]*rate2     #New recovered individuals
P[1:N_age] <- S[i] + E[i] + I[i] + R[i] + V[i] #Total population
P_tot <- sum(P)
inv_P[1:N_age] <- 1.0/P[i]
F_S[1:N_age] <- S[i]*inv_P[i] #Susceptible fraction
F_R[1:N_age] <- R[i]*inv_P[i] #Recovered fraction
F_V[1:N_age] <- V[i]*inv_P[i] #Vaccinated fraction
vacc_rate[1:N_age]=vacc_rate_annual[i,as.integer(year_i)]*S[i]*vaccine_efficacy*dt

#equations
update(day) <- day + dt
update(year) <- year_i + year0 - 1
update(FOI_total) <- FOI_sum
update(S[1]) <- max(Pmin,S[1] - E_new[1] - vacc_rate[1] + dP1[1] - (dP2[1]*F_S[1]))
update(S[2:N_age]) <- max(Pmin,S[i] - E_new[i] - vacc_rate[i] + (dP1[i]*F_S[i-1]) - (dP2[i]*F_S[i]))
update(E[1:N_age]) <- max(Pmin,E[i] + E_new[i] - I_new[i])
update(I[1:N_age]) <- max(Pmin,I[i] + I_new[i] - R_new[i])
update(R[1]) <- max(Pmin,R[1] + R_new[1] - (dP2[1]*F_R[1]))
update(R[2:N_age]) <- max(Pmin,R[i] + R_new[i] + (dP1[i]*F_R[i-1]) - (dP2[i]*F_R[i]))
update(V[1]) <- max(Pmin,V[1] + vacc_rate[1] - (dP2[1]*F_V[1]))
update(V[2:N_age]) <- max(Pmin,V[i] + vacc_rate[i] + (dP1[i]*F_V[i-1]) - (dP2[i]*F_V[i]))
update(C[1:N_age]) <- I_new[i]

#initial
initial(day) <- 0
initial(year) <- year0-1
initial(FOI_total) <- FOI_spillover
initial(S[1:N_age]) <- Sus0[i]
initial(E[1:N_age]) <- Exp0[i]
initial(I[1:N_age]) <- Inf0[i]
initial(R[1:N_age]) <- Rec0[i]
initial(V[1:N_age]) <- Vac0[i]
initial(C[1:N_age]) <- Cas0[i]

#dimensions
dim(S) <- N_age
dim(E) <- N_age
dim(I) <- N_age
dim(R) <- N_age
dim(V) <- N_age
dim(C) <- N_age

dim(E_new) <- N_age
dim(I_new) <- N_age
dim(R_new) <- N_age
dim(F_S) <- N_age
dim(F_R) <- N_age
dim(F_V) <- N_age
dim(dP1)<-N_age
dim(dP2)<-N_age
dim(P) <- N_age
dim(inv_P)<-N_age
dim(vacc_rate)<-N_age

dim(Sus0) <- N_age
dim(Exp0) <- N_age
dim(Inf0) <- N_age
dim(Rec0) <- N_age
dim(Vac0) <- N_age
dim(Cas0) <- N_age
dim(dP1_all) <- c(N_age, n_years)
dim(dP2_all) <- c(N_age, n_years)
dim(vacc_rate_annual) <- c(N_age, n_years)
