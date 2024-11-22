# SEIRV model with functionality to split infections into sylvatic and urban





time_inc <- parameter() #Time increment in days
initial(day) <- time_inc #Initial value of time in days
update(day) <- day + time_inc

#Parameters---------------------------------------------------------------------
t_incubation <- parameter() #Length in days of yellow fever incubation period in mosquito vectors
t_latent <- parameter() #Length in days of latent period in humans exposed to yellow dever
t_infectious <- parameter() #Length of infectious period in humans with yellow fever
FOI_spillover <- parameter() #Spillover force of infection (per day)
R0 <- parameter() #Basic reproduction number for human-human transmission
N_age <- parameter() #Number of age categories
vacc_rate_daily <- parameter() #Daily rate of vaccination by age and year
vaccine_efficacy <- parameter() #Proportion of vaccinations which successfully protect the recipient







#Initial conditions-------------------------------------------------------------
year0 <- parameter()  #Starting year
S_0 <- parameter() #Susceptible population by age group at start
E_0 <- parameter() #Exposed population by age group at start

I_0 <- parameter() #Infectious population by age group at start

R_0 <- parameter() #Recovered population by age group at start
V_0 <- parameter() #Vaccinated population by age group at start
dP1_all <- parameter() #Daily increase in number of people by age group (people arriving in group due to age etc.)
dP2_all <- parameter() #Daily decrease in number of people by age group (people leaving group due to age etc.)
n_years <- parameter() #Number of years for which model to be run

Pmin <- 1.0e-99 #Minimum population setting to avoid negative numbers
FOI_max <- 1.0 #Upper threshold for total force of infection to avoid more infections than people in a group
rate1 <- time_inc/(t_incubation+t_latent)
rate2 <- time_inc/t_infectious





beta <- (R0*time_inc)/t_infectious #Daily exposure rate
FOI_sylv_cur <- min(FOI_max,FOI_spillover*time_inc)
FOI_urb_cur <- min(FOI_max,beta*((sum(I_sylv)+sum(I_urb))/P_tot))
year_i <- floor(day/365)+1 #Number of years since start, as integer

dP1[1:N_age] <- dP1_all[i, year_i]*time_inc #Increase in population by age group over 1 time increment
dP2[1:N_age] <- dP2_all[i, year_i]*time_inc #Decrease in population by age group over 1 time increment

E_new_sylv[1:N_age] <- Binomial(as.integer(S[i]), FOI_sylv_cur) #New exposed individuals by age group (sylvatic)
E_new_urb[1:N_age] <- Binomial(as.integer(S[i]), FOI_urb_cur) #New exposed individuals by age group (urban)
I_new_sylv[1:N_age] <- E_sylv[i]*rate1     #New infectious individuals by age group (sylvatic)
I_new_urb[1:N_age] <- E_urb[i]*rate1     #New infectious individuals by age group (urban)
R_new_sylv[1:N_age] <- I_sylv[i]*rate2     #New recovered individuals by age group (sylvatic)
R_new_urb[1:N_age] <- I_urb[i]*rate2     #New recovered individuals by age group (urban)
P_nV[1:N_age] <- S[i] + R[i] #Total vaccine-targetable population by age group
inv_P_nV[1:N_age] <- 1.0/P_nV[i]
P[1:N_age] <- P_nV[i] + V[i] #Total population by age group (excluding E+I)
P_tot <- sum(P) #Total overall population (excluding E+I)
inv_P[1:N_age] <- 1.0/P[i]

vacc_rate[1:N_age] <- vacc_rate_daily[i,year_i]*vaccine_efficacy*time_inc*P[i] #Total no. vaccinations by age






#Updates to output values at each time increment--------------------------------
update(year) <- year_i + year0 - 1
update(FOI_sylv) <- FOI_sylv_cur
update(FOI_urb) <- FOI_urb_cur





update(S[1]) <- max(Pmin,S[1] - E_new_sylv[1] - E_new_urb[1] - vacc_rate[1]*S[1]*inv_P_nV[1] + dP1[1] - (dP2[1]*S[1]*inv_P[1]))
update(S[2:N_age]) <- max(Pmin,S[i]-E_new_sylv[i]-E_new_urb[i]-vacc_rate[i]*S[i]*inv_P_nV[i]+(dP1[i]*S[i-1]*inv_P[i-1])-(dP2[i]*S[i]*inv_P[i]))
update(E_sylv[1:N_age]) <- max(Pmin,E_sylv[i] + E_new_sylv[i] - I_new_sylv[i])
update(E_urb[1:N_age]) <- max(Pmin,E_urb[i] + E_new_urb[i] - I_new_urb[i])

update(I_sylv[1:N_age]) <- max(Pmin,I_sylv[i] + I_new_sylv[i] - R_new_sylv[i])
update(I_urb[1:N_age]) <- max(Pmin,I_urb[i] + I_new_urb[i] - R_new_urb[i])

update(R[1]) <- max(Pmin,R[1] + R_new_sylv[1] + R_new_urb[1] - vacc_rate[1]*R[1]*inv_P_nV[1] - (dP2[1]*R[1]*inv_P[1]))
update(R[2:N_age]) <- max(Pmin,R[i]+R_new_sylv[i]+R_new_urb[i]-vacc_rate[i]*R[i]*inv_P_nV[i]+(dP1[i]*R[i-1]*inv_P[i-1])-(dP2[i]*R[i]*inv_P[i]))
update(V[1]) <- max(Pmin,V[1] + vacc_rate[1] - (dP2[1]*V[1]*inv_P[1]))
update(V[2:N_age]) <- max(Pmin,V[i] + vacc_rate[i] + (dP1[i]*V[i-1]*inv_P[i-1]) - (dP2[i]*V[i]*inv_P[i]))
update(C_sylv[1:N_age]) <- I_new_sylv[i]
update(C_urb[1:N_age]) <- I_new_urb[i]

#Initial values-----------------------------------------------------------------
initial(year) <- year0-1
initial(FOI_sylv) <- FOI_spillover
initial(FOI_urb) <- 0





initial(S[1:N_age]) <- S_0[i]
initial(E_sylv[1:N_age]) <- E_0[i]
initial(E_urb[1:N_age]) <- 0
initial(I_sylv[1:N_age]) <- I_0[i]
initial(I_urb[1:N_age]) <- 0
initial(R[1:N_age]) <- R_0[i]
initial(V[1:N_age]) <- V_0[i]
initial(C_sylv[1:N_age]) <- 0
initial(C_urb[1:N_age]) <- 0

#Dimensions---------------------------------------------------------------------
dim(S) <- N_age
dim(E_sylv) <- N_age
dim(E_urb) <- N_age
dim(I_sylv) <- N_age
dim(I_urb) <- N_age
dim(R) <- N_age
dim(V) <- N_age
dim(C_sylv) <- N_age
dim(C_urb) <- N_age
dim(dP1)<-N_age
dim(dP2)<-N_age
dim(E_new_sylv) <- N_age
dim(E_new_urb) <- N_age
dim(I_new_sylv) <- N_age
dim(I_new_urb) <- N_age
dim(R_new_sylv) <- N_age
dim(R_new_urb) <- N_age

dim(P_nV) <- N_age
dim(inv_P_nV) <- N_age
dim(P) <- N_age
dim(inv_P) <- N_age
dim(vacc_rate) <- N_age

dim(S_0) <- N_age
dim(E_0) <- N_age

dim(I_0) <- N_age

dim(R_0) <- N_age
dim(V_0) <- N_age
dim(dP1_all) <- c(N_age, n_years)
dim(dP2_all) <- c(N_age, n_years)
dim(vacc_rate_daily) <- c(N_age, n_years)



