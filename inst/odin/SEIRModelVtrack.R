# New version of YF SEIRV model where vaccinated individuals are split off as a separate
# track, to facilitate incorporation of waning immunity and related factors
# Currently just splits unvaccinated and vaccinated into 2 tracks and applies
# vaccine_efficacy as an immunity factor to 'susceptible' vaccinated population
# To add - see end

time_inc <- parameter() #Time increment in days
initial(day) <- time_inc #Initial value of time in days
update(day) <- day + time_inc

#Parameters---------------------------------------------------------------------
t_incubation <- parameter() #Length in days of yellow fever incubation period in mosquito vectors
t_latent <- parameter() #Length in days of latent period in humans exposed to yellow fever
t_infectious <- parameter() #Length of infectious period in humans with yellow fever
FOI_spillover <- parameter() #Spillover force of infection (per day) at each time point
R0 <- parameter() #Basic reproduction number for human-human transmission at each time point
N_age <- parameter() #Number of age categories
vacc_rate_daily <- parameter() #Daily rate of vaccination by age and year
vaccine_efficacy <- parameter() #Proportion of vaccinations which successfully protect the recipient







#Initial conditions-------------------------------------------------------------
year0 <- parameter()  #Starting year
S_0 <- parameter() #Susceptible population by age group at start
E_0 <- parameter() #Exposed population by age group at start
I_0 <- parameter() #Infectious population by age group at start
R_0 <- parameter() #Recovered population by age group at start



dP1_all <- parameter() #Daily increase in number of people by age group (people arriving in group due to age etc.)
dP2_all <- parameter() #Daily decrease in number of people by age group (people leaving group due to age etc.)
n_years <- parameter() #Number of years for which model to be run
n_t_pts <- parameter() #Total number of time points
Pmin <- 1.0e-99 #Minimum population setting to avoid negative numbers
FOI_max <- 1.0 #Upper threshold for total force of infection to avoid more infections than people in a group
rate1 <- time_inc/(t_incubation+t_latent) # Rate of transfer from E to I
rate2 <- time_inc/t_infectious # Rate of transfer from I to R




#Values calculated--------------------------------------------------------------
t_pt <- day/time_inc #Number of time points passed
beta <- (R0[t_pt]*time_inc)/t_infectious #Daily exposure rate
FOI_sum <-  min(FOI_max,beta*(sum(I)/P_tot) + (FOI_spillover[t_pt]*time_inc)) #Total force of infection

year_i <- floor(day/365)+1 #Number of years since start, as integer

dP1[1:N_age] <- dP1_all[i, year_i]*time_inc #Increase in population by age group over 1 time increment
dP2[1:N_age] <- dP2_all[i, year_i]*time_inc #Decrease in population by age group over 1 time increment

P_nV[1:N_age] <- S[i,1] + R[i,1] #Total vaccine-targetable population by age group (NB double-vaccination already accounted for)
inv_P_nV[1:N_age] <- 1.0/P_nV[i]
P[1:N_age] <- P_nV[i] + S[i,2] + R[i,2] #Total population by age group (excluding E+I)
P_tot <- sum(P) #Total overall population (excluding E+I)
inv_P[1:N_age] <- 1.0/P[i]
vacc_rate[1:N_age] <- vacc_rate_daily[i,year_i]*time_inc*P[i] #Total no. vaccinations by age

E_new[1:N_age,1] <- Binomial(as.integer(S[i,1]), FOI_sum) #New exposed unvaccinated individuals by age group
E_new[1:N_age,2] <- Binomial(as.integer(S[i,2]), FOI_sum*(1.0-vaccine_efficacy)) #New exposed vaccinated individuals by age group
I_new[1:N_age,1:2] <- E[i,j]*rate1     #New infectious individuals by age group
R_new[1:N_age,1:2] <- I[i,j]*rate2     #New recovered individuals by age group

S_new_V[1:N_age] <- vacc_rate[i]*S[i,1]*inv_P_nV[i]
R_new_V[1:N_age] <- vacc_rate[i]*R[i,1]*inv_P_nV[i]





#Updates to output values at each time increment--------------------------------
update(year) <- year_i + year0 - 1
update(FOI_total) <- FOI_sum






update(S[1,1]) <- max(Pmin,S[1,1] - E_new[1,1] - S_new_V[1] + dP1[1] - (dP2[1]*S[1,1]*inv_P[1]))
update(S[2:N_age,1]) <- max(Pmin,S[i,1] - E_new[i,1] - S_new_V[i] + (dP1[i]*S[i-1,1]*inv_P[i-1]) - (dP2[i]*S[i,1]*inv_P[i]))
update(S[1,2]) <- max(Pmin,S[1,2] - E_new[1,2] + S_new_V[1] - (dP2[1]*S[1,2]*inv_P[1]))
update(S[2:N_age,2]) <- max(Pmin,S[i,2] - E_new[i,2] + S_new_V[i] + (dP1[i]*S[i-1,2]*inv_P[i-1]) - (dP2[i]*S[i,2]*inv_P[i]))
update(E[1:N_age,1:2]) <- max(Pmin,E[i,j] + E_new[i,j] - I_new[i,j])
update(I[1:N_age,1:2]) <- max(Pmin,I[i,j] + I_new[i,j] - R_new[i,j])
update(R[1,1]) <- max(Pmin,R[1,1] + R_new[1,1] - R_new_V[1] - (dP2[1]*R[1,1]*inv_P[1]))
update(R[2:N_age,1]) <- max(Pmin,R[i,1] + R_new[i,1] - R_new_V[i] + (dP1[i]*R[i-1,1]*inv_P[i-1]) - (dP2[i]*R[i,1]*inv_P[i]))
update(R[1,2]) <- max(Pmin,R[1,2] + R_new[1,2] + R_new_V[1] - (dP2[1]*R[1,2]*inv_P[1]))
update(R[2:N_age,2]) <- max(Pmin,R[i,2] + R_new[i,2] - R_new_V[i] + (dP1[i]*R[i-1,2]*inv_P[i-1]) - (dP2[i]*R[i,2]*inv_P[i]))

update(C[1:N_age]) <- I_new[i,1] + I_new[i,2]



#Initial values-----------------------------------------------------------------
initial(year) <- year0
initial(FOI_total) <- FOI_spillover[1]






initial(S[1:N_age,1:2]) <- S_0[i,j]
initial(E[1:N_age,1:2]) <- E_0[i,j]

initial(I[1:N_age,1:2]) <- I_0[i,j]

initial(R[1:N_age,1:2]) <- R_0[i,j]

initial(C[1:N_age]) <- 0


#Dimensions---------------------------------------------------------------------
dim(S) <- c(N_age,2)
dim(E) <- c(N_age,2)
dim(I) <- c(N_age,2)
dim(R) <- c(N_age,2)

dim(C) <- N_age

dim(dP1)<-N_age
dim(dP2)<-N_age

dim(E_new) <- c(N_age,2)
dim(I_new) <- c(N_age,2)
dim(R_new) <- c(N_age,2)

dim(S_new_V) <- N_age
dim(R_new_V) <- N_age

dim(P_nV) <- N_age
dim(inv_P_nV) <- N_age
dim(P) <- N_age
dim(inv_P) <- N_age
dim(vacc_rate) <- N_age

dim(S_0) <- c(N_age,2)
dim(E_0) <- c(N_age,2)
dim(I_0) <- c(N_age,2)
dim(R_0) <- c(N_age,2)

dim(dP1_all) <- c(N_age, n_years)
dim(dP2_all) <- c(N_age, n_years)
dim(vacc_rate_daily) <- c(N_age, n_years)
dim(FOI_spillover) <- n_t_pts
dim(R0) <- n_t_pts

# Possible things to incorporate:
# -Split imperfect immunity and vaccine misreporting into multiple parameters
#  (it may be appropriate to apply vaccine misreporting when vacc_data is calculated)
# -Add one or more additional tracks for different vaccination groups
#  (1 dose vs 2 doses, childhood vs adult vaccination)
# -Incorporate waning of immunity over time
# -Incorporate age-dependent vaccination-associated immunity?
# -Allow recovered individuals to potentially be infected with their own immunity factor <=1?
