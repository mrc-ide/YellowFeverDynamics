# SEIRV (Susceptible, Exposed, Infectious, Recovered, Vaccinated) yellow fever model,
# incorporating the force of infection of spillover from sylvatic/non-human primate reservoirs (which can also
# represent case importation) and the reproduction number for human-human transmission. Returns SEIRV data at each
# time point (separated by increment dt) and also numbers of new infections and total force of infection at each
# time point. This version splits new infections into sylvatic and urban.

dt <- user() #Time increment in days
initial(time) <- 0 #Initial value of time in days
update(time) <- (step + 1) * dt

#Parameters
t_incubation <- user() #TBA
t_latent <- user() #TBA
t_infectious <- user() #TBA
FOI_spillover <- user() #Spillover force of infection
R0 <- user() #Basic reproduction number
N_age <- user() #Number of age categories
vacc_rate_annual[,] <- user() #Daily rate of vaccination by age and year
vaccine_efficacy <- user() #Proportion of vaccinations which successfully protect the recipient

#Initial conditions
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

Pmin <- 1.0e-99 #Minimum population setting to avoid negative numbers
FOI_max <- 1.0 #Upper threshold for total force of infection to avoid more infections than people in a group
rate1 <- dt/(t_incubation+t_latent)
rate2 <- dt/t_infectious
beta <- (R0*dt)/t_infectious #Daily exposure rate
FOI_sylvatic <- min(FOI_max,FOI_spillover*dt)
FOI_urban <- min(FOI_max,beta*((sum(I_sylvatic)+sum(I_urban))/P_tot))
year_i <- floor((step*dt)/365) + 1 #Number of years since start, as integer
dP1[1:N_age] <- dP1_all[i, as.integer(year_i)]*dt #Increase in population by age group over 1 time increment
dP2[1:N_age] <- dP2_all[i, as.integer(year_i)]*dt #Decrease in population by age group over 1 time increment

E_new_sylvatic[1:N_age] <- rbinom(as.integer(S[i]), FOI_sylvatic) #New exposed individuals by age group (sylvatic)
E_new_urban[1:N_age] <- rbinom(as.integer(S[i]), FOI_urban) #New exposed individuals by age group (urban)
I_new_sylvatic[1:N_age] <- E_sylvatic[i]*rate1     #New infectious individuals by age group (sylvatic)
I_new_urban[1:N_age] <- E_urban[i]*rate1     #New infectious individuals by age group (urban)
R_new_sylvatic[1:N_age] <- I_sylvatic[i]*rate2     #New recovered individuals by age group (sylvatic)
R_new_urban[1:N_age] <- I_urban[i]*rate2     #New recovered individuals by age group (urban)
P_nV[1:N_age] <- S[i] + R[i] #Total vaccine-targetable population by age group
inv_P_nV[1:N_age] <- 1.0/P_nV[i]
P[1:N_age] <- S[i] + E_sylvatic[i] + E_urban[i] + I_sylvatic[i] + I_urban[i] + R[i] + V[i]#Total population by age group
P_tot <- sum(P) #Total overall population
inv_P[1:N_age] <- 1.0/P[i]
vacc_rate[1:N_age] <- vacc_rate_annual[i,as.integer(year_i)]*vaccine_efficacy*dt*P[i] #Total no. vaccinations by age

#Updates to output values at each time increment
update(year) <- year_i + year0 - 1
update(FOI_sylv) <- FOI_sylvatic
update(FOI_urb) <- FOI_urban
update(S[1]) <- max(Pmin,S[1] - E_new_sylvatic[1] - E_new_urban[1] - vacc_rate[1]*S[1]*inv_P_nV[1] + dP1[1] - (dP2[1]*S[1]*inv_P[1]))
update(S[2:N_age]) <- max(Pmin,S[i] - E_new_sylvatic[i] - E_new_urban[i] - vacc_rate[i]*S[i]*inv_P_nV[i] + (dP1[i]*S[i-1]*inv_P[i-1]) - (dP2[i]*S[i]*inv_P[i]))
update(E_sylvatic[1:N_age]) <- max(Pmin,E_sylvatic[i] + E_new_sylvatic[i] - I_new_sylvatic[i])
update(E_urban[1:N_age]) <- max(Pmin,E_urban[i] + E_new_urban[i] - I_new_urban[i])
update(I_sylvatic[1:N_age]) <- max(Pmin,I_sylvatic[i] + I_new_sylvatic[i] - R_new_sylvatic[i])
update(I_urban[1:N_age]) <- max(Pmin,I_urban[i] + I_new_urban[i] - R_new_urban[i])
update(R[1]) <- max(Pmin,R[1] + R_new_sylvatic[1] + R_new_urban[1] - vacc_rate[1]*R[1]*inv_P_nV[1] - (dP2[1]*R[1]*inv_P[1]))
update(R[2:N_age]) <- max(Pmin,R[i] + R_new_sylvatic[i] + R_new_urban[i] - vacc_rate[i]*R[i]*inv_P_nV[i] + (dP1[i]*R[i-1]*inv_P[i-1]) - (dP2[i]*R[i]*inv_P[i]))
update(V[1]) <- max(Pmin,V[1] + vacc_rate[1] - (dP2[1]*V[1]*inv_P[1]))
update(V[2:N_age]) <- max(Pmin,V[i] + vacc_rate[i] + (dP1[i]*V[i-1]*inv_P[i-1]) - (dP2[i]*V[i]*inv_P[i]))
update(C_sylvatic[1:N_age]) <- I_new_sylvatic[i]
update(C_urban[1:N_age]) <- I_new_urban[i]

#Initial values
initial(year) <- year0-1
initial(FOI_sylv) <- FOI_spillover
initial(FOI_urb) <- 0
initial(S[1:N_age]) <- Sus0[i]
initial(E_sylvatic[1:N_age]) <- Exp0[i]
initial(E_urban[1:N_age]) <- 0
initial(I_sylvatic[1:N_age]) <- Inf0[i]
initial(I_urban[1:N_age]) <- 0
initial(R[1:N_age]) <- Rec0[i]
initial(V[1:N_age]) <- Vac0[i]
initial(C_sylvatic[1:N_age]) <- Cas0[i]
initial(C_urban[1:N_age]) <- 0

#Dimensions
dim(S) <- N_age
dim(E_sylvatic) <- N_age
dim(E_urban) <- N_age
dim(I_sylvatic) <- N_age
dim(I_urban) <- N_age
dim(R) <- N_age
dim(V) <- N_age
dim(C_sylvatic) <- N_age
dim(C_urban) <- N_age

dim(dP1)<-N_age
dim(dP2)<-N_age
dim(E_new_sylvatic) <- N_age
dim(E_new_urban) <- N_age
dim(I_new_sylvatic) <- N_age
dim(I_new_urban) <- N_age
dim(R_new_sylvatic) <- N_age
dim(R_new_urban) <- N_age
dim(P_nV) <- N_age
dim(inv_P_nV) <- N_age
dim(P) <- N_age
dim(inv_P) <- N_age
dim(vacc_rate) <- N_age

dim(Sus0) <- N_age
dim(Exp0) <- N_age
dim(Inf0) <- N_age
dim(Rec0) <- N_age
dim(Vac0) <- N_age
dim(Cas0) <- N_age
dim(dP1_all) <- c(N_age, n_years)
dim(dP2_all) <- c(N_age, n_years)
dim(vacc_rate_annual) <- c(N_age, n_years)
