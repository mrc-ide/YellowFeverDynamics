# SEIRV model with rates of E->I and I->R transfer replaced by fixed delay time
# Updated based on YEP version 0.3 - time-varying epi parameters, multiple regions


#Parameters---------------------------------------------------------------------
time_inc <- parameter() #Time increment in days
n_regions <- parameter() #number of regions


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

t_pt <- day/time_inc #Number of time points passed
beta[1:n_regions] <- (R0[i, t_pt]*time_inc)/t_infectious #Daily exposure rate
FOI_sum[1:n_regions] <-  min(FOI_max, beta[i]*(sum(I[i,])/P_tot[i]) + (FOI_spillover[i, t_pt]*time_inc)) #Total force of infection
year_i <- floor(day/365)+1 #Number of years since start, as integer


dP1[1:n_regions,1:N_age] <- dP1_all[i,j,year_i]*time_inc #Increase in population by age group over 1 time increment
dP2[1:n_regions,1:N_age] <- dP2_all[i,j,year_i]*time_inc #Decrease in population by age group over 1 time increment
E_new[1:n_regions,1:N_age] <- Binomial(as.integer(S[i,j]), FOI_sum[i]) #New exposed individuals by age group
I_new[1:n_regions,1:N_age] <- E_delay[i,as.integer(j+di1)]     #New infectious individuals by age group
R_new[1:n_regions,1:N_age] <- I_delay[i,as.integer(j+di2)]     #New recovered individuals by age group





P_nV[1:n_regions,1:N_age] <- S[i,j] + R[i,j] #Total vaccine-targetable population by age group
inv_P_nV[1:n_regions,1:N_age] <- 1.0/P_nV[i,j]
P[1:n_regions,1:N_age] <- P_nV[i,j] + V[i,j] #Total population by age group (excluding E+I)
P_tot[1:n_regions] <- sum(P[i, ]) #Total overall population (excluding E+I)
inv_P[1:n_regions,1:N_age] <- 1.0/P[i,j]

vacc_rate[1:n_regions,1:N_age] <- vacc_rate_daily[i,j,year_i]*vaccine_efficacy*time_inc*P[i,j] #Total no. vaccinations by age

#Updates to output values at each time increment--------------------------------
update(day) <- day + time_inc
update(year) <- year_i + year0 - 1
update(FOI_total[1:n_regions]) <- FOI_sum[i]

update(S[1:n_regions,1]) <- max(Pmin, S[i,1] - E_new[i,1] - vacc_rate[i,1]*S[i,1]*inv_P_nV[i,1] + dP1[i,1] - (dP2[i,1]*S[i,1]*inv_P[i,1]))
update(S[1:n_regions,2:N_age]) <- max(Pmin, S[i,j] - E_new[i,j] - vacc_rate[i,j]*S[i,j]*inv_P_nV[i,j] + (dP1[i,j]*S[i,j-1]*inv_P[i,j-1]) - (dP2[i,j]*S[i,j]*inv_P[i,j]))
update(E[1:n_regions,1:N_age]) <- max(Pmin, E[i,j] + E_new[i,j] - I_new[i,j])
update(E_delay[1:n_regions,(N_age+1):np_E_delay]) <- E_delay[i,j-N_age]
update(E_delay[1:n_regions,1:N_age]) <- E_new[i,j]
update(I[1:n_regions,1:N_age]) <- max(Pmin, I[i,j] + I_new[i,j] - R_new[i,j])
update(I_delay[1:n_regions,(N_age+1):np_I_delay]) <- I_delay[i,j-N_age]
update(I_delay[1:n_regions,1:N_age]) <- I_new[i,j]
update(R[1:n_regions,1]) <- max(Pmin, R[i,1] + R_new[i,1] - vacc_rate[i,1]*R[i,1]*inv_P_nV[i,1] - (dP2[i,1]*R[i,1]*inv_P[i,1]))
update(R[1:n_regions,2:N_age]) <- max(Pmin, R[i,j] + R_new[i,j] - vacc_rate[i,j]*R[i,j]*inv_P_nV[i,j] + (dP1[i,j]*R[i,j-1]*inv_P[i,j-1]) - (dP2[i,j]*R[i,j]*inv_P[i,j]))
update(V[1:n_regions,1]) <- max(Pmin, V[i,1] + vacc_rate[i,1] - (dP2[i,1]*V[i,1]*inv_P[i,1]))
update(V[1:n_regions,2:N_age]) <- max(Pmin, V[i,j] + vacc_rate[i,j] + (dP1[i,j]*V[i,j-1]*inv_P[i,j-1]) - (dP2[i,j]*V[i,j]*inv_P[i,j]))
update(C[1:n_regions,1:N_age]) <- I_new[i,j]




















#Initial values-----------------------------------------------------------------
initial(day) <- time_inc
initial(year) <- year0
initial(FOI_total[1:n_regions]) <- FOI_spillover[i,1]
initial(S[1:n_regions,1:N_age]) <- S_0[i,j]
initial(E[1:n_regions,1:N_age]) <- E_0[i,j]
initial(E_delay[1:n_regions,1:np_E_delay]) <- E_delay0[i,j]
initial(I[1:n_regions,1:N_age]) <- I_0[i,j]
initial(I_delay[1:n_regions,1:np_I_delay]) <- I_delay0[i,j]
initial(R[1:n_regions,1:N_age]) <- R_0[i,j]
initial(V[1:n_regions,1:N_age]) <- V_0[i,j]
initial(C[1:n_regions,1:N_age]) <- 0











#Dimensions---------------------------------------------------------------------
#Updated values
dim(FOI_total) <- n_regions
dim(S) <- c(n_regions, N_age)
dim(E) <- c(n_regions, N_age)
dim(E_delay) <- c(n_regions,np_E_delay)
dim(I) <- c(n_regions, N_age)
dim(I_delay) <- c(n_regions,np_I_delay)
dim(R) <- c(n_regions, N_age)
dim(V) <- c(n_regions, N_age)
dim(C) <- c(n_regions, N_age)

dim(beta) <- n_regions
dim(FOI_sum) <- n_regions
dim(dP1) <- c(n_regions, N_age)
dim(dP2) <- c(n_regions, N_age)
dim(E_new) <- c(n_regions, N_age)
dim(I_new) <- c(n_regions, N_age)
dim(R_new) <- c(n_regions, N_age)

dim(P_nV) <- c(n_regions, N_age)
dim(inv_P_nV) <- c(n_regions, N_age)
dim(P) <- c(n_regions, N_age)
dim(P_tot) <- n_regions
dim(inv_P) <- c(n_regions, N_age)
dim(vacc_rate) <- c(n_regions, N_age)

dim(S_0) <- c(n_regions, N_age)
dim(E_0) <- c(n_regions, N_age)
dim(E_delay0) <- c(n_regions,np_E_delay)
dim(I_0) <- c(n_regions, N_age)
dim(I_delay0) <- c(n_regions,np_I_delay)
dim(R_0) <- c(n_regions, N_age)
dim(V_0) <- c(n_regions, N_age)
dim(dP1_all) <- c(n_regions, N_age, n_years)
dim(dP2_all) <- c(n_regions, N_age, n_years)
dim(vacc_rate_daily) <- c(n_regions, N_age, n_years)
dim(FOI_spillover) <- c(n_regions, n_t_pts)
dim(R0) <- c(n_regions, n_t_pts)
