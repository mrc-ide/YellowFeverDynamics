# SEIRV model with functionality to split infections into sylvatic and urban





dt <- user() #Time increment in days
initial(time) <- 0 #Initial value of time in days
update(time) <- time + dt

#Parameters---------------------------------------------------------------------
t_incubation <- user() #TBA
t_latent <- user() #TBA
t_infectious <- user() #TBA
FOI_spillover <- user() #Spillover force of infection (per day)
R0 <- user() #Basic reproduction number
N_age <- user() #Number of age categories
vacc_rate_daily[,] <- user() #Daily rate of vaccination by age and year
vaccine_efficacy <- user() #Proportion of vaccinations which successfully protect the recipient







#Initial conditions-------------------------------------------------------------
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
FOI_sylv_cur <- min(FOI_max,FOI_spillover*dt)
FOI_urb_cur <- min(FOI_max,beta*((sum(I_sylv)+sum(I_urb))/P_tot))
year_i <- floor(((step+1)*dt)/365) + 1 #Number of years since start, as integer

dP1[1:N_age] <- dP1_all[i, as.integer(year_i)]*dt #Increase in population by age group over 1 time increment
dP2[1:N_age] <- dP2_all[i, as.integer(year_i)]*dt #Decrease in population by age group over 1 time increment

E_new_sylv[1:N_age] <- rbinom(as.integer(S[i]), FOI_sylv_cur) #New exposed individuals by age group (sylvatic)
E_new_urb[1:N_age] <- rbinom(as.integer(S[i]), FOI_urb_cur) #New exposed individuals by age group (urban)
I_new_sylv[1:N_age] <- E_sylv[i]*rate1     #New infectious individuals by age group (sylvatic)
I_new_urb[1:N_age] <- E_urb[i]*rate1     #New infectious individuals by age group (urban)
R_new_sylv[1:N_age] <- I_sylv[i]*rate2     #New recovered individuals by age group (sylvatic)
R_new_urb[1:N_age] <- I_urb[i]*rate2     #New recovered individuals by age group (urban)
P_nV[1:N_age] <- S[i] + R[i] #Total vaccine-targetable population by age group
inv_P_nV[1:N_age] <- 1.0/P_nV[i]
P[1:N_age] <- P_nV[i] + V[i] #Total population by age group (excluding E+I)
P_tot <- sum(P) #Total overall population (excluding E+I)
inv_P[1:N_age] <- 1.0/P[i]
vacc_rate[1:N_age] <- vacc_rate_daily[i,as.integer(year_i)]*vaccine_efficacy*dt*P[i] #Total no. vaccinations by age







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







initial(S[1:N_age]) <- Sus0[i]
initial(E_sylv[1:N_age]) <- Exp0[i]
initial(E_urb[1:N_age]) <- 0
initial(I_sylv[1:N_age]) <- Inf0[i]
initial(I_urb[1:N_age]) <- 0
initial(R[1:N_age]) <- Rec0[i]
initial(V[1:N_age]) <- Vac0[i]
initial(C_sylv[1:N_age]) <- Cas0[i]
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

dim(Sus0) <- N_age
dim(Exp0) <- N_age
dim(Inf0) <- N_age
dim(Rec0) <- N_age
dim(Vac0) <- N_age
dim(Cas0) <- N_age
dim(dP1_all) <- c(N_age, n_years)
dim(dP2_all) <- c(N_age, n_years)
dim(vacc_rate_daily) <- c(N_age, n_years)


