######################################################## #
#I. SEIRV_age model for RSV with optimisation 
######################################################## #

rm(list=ls())
#setwd("/home/phuong/phuonght/gitKraken_project/RSV")
library(deSolve)
library(ggplot2)

#Changeable variables  ##########################
vac.frac = 0

######################################### #
# MODEL SETTINGS                     ----
######################################### #
pop_size          <- 1100000 # Belgium pop size
num_days          <- 500
num_weeks         <- 52*20
num_days_infected <- 10  #[ 8-11], 6.7 is original estimate
num_days_exposed  <- 4   # [2,6]
num_days_waning   <- 200 # [148, 164] best fit
infected_seeds    <- 10 # ? 1= 1000

######################################### #
# INITIALIZE PARAMETERS AND POPULATION ----
######################################### #
#https://steemit.com/science/@fouad/x-history-of-math-symbols

#  population states
S1 <- 1/80 # less than 1 yrs old
E1 <- infected_seeds/2/(pop_size)
I1 <- 0
R1 <- 0

V  <- 0

E2 <- (infected_seeds/2)/(pop_size) # greater than 1 yrs old #LW
I2 <- 0
R2 <- 0
S2 <- 1-S1-E1-I1-R1-E2-I2-R2-V


######################################### #
# SET FUNCTION PARAMETERS            ----
######################################### #
# set time frame
times      <- seq(0, num_weeks, by = 1)

# set initial health states
states     <- c(S1 = S1,E1 = E1, I1 = I1, R1 = R1, V = V,
                S2 = S2,E2 = E2, I2 = I2, R2 = R2)

# set parameters

params     <- c(sigma = 1/0.57, #rate of movement from latent to infectious stage
                gamma = 1/1.4, # recovery rate
                # sigma = x[1], #rate of movement from latent to infectious stage
                # gamma = x[2], # recovery rate
                # nu= 0.044,     # rate of loss of immunity
                nu=  7/ num_days_waning,     # rate of loss of immunity
                eta1 = 1/(1*52),             # aging rate
                eta2 = 1/(79*52),            # aging rate
                # eta = 1/(80*52),           # death rate in general
                mu = 1/(79*52),              # birth rate
                beta1 = 0.89,                # the degree of seasonality, range [0,1],higher value stronger seasonal drivers
                phi = 2.43,                  # phase shift?
                beta0 = 1.99,                # average transmission rate
                # R0 = 3,                    # reproduction number
                delta = 0.54,                # scaled susceptibility
                alpha = 0.76,                # scaled infectiousness
                p = vac.frac                 # vaccination proportion
                )

######################################### #
# CREATE SIRV FUNCTION               ----
######################################### #
# Note: updating the health states over time and keeping track of the changes, is handled by the
# ode-function of the 'deSolve' package.
seirv_func <- function(t, states, params) {

  with(as.list(c(states, params)),
       {
         # define variables
         # P <- (S1+E1+I1+R1+S2+E2+I2+R2)
         beta <- beta0*(1+beta1*cos(2*pi*t/52+phi))


         # calculate state changes
         dS1 <- (1-p)* mu - beta*S1*(I1+alpha*I2) - eta1*S1             + nu*R1
         dE1 <-             beta*S1*(I1+alpha*I2) - eta1*E1 - sigma*E1
         dI1 <-                                   - eta1*I1 + sigma*E1  - gamma*I1
         dR1 <-                                   - eta1*R1             + gamma*I1  - nu*R1

         dV <- p*mu  - nu*V

         dS2 <- nu*V + eta1*S1 - delta*beta*S2*(I1+alpha*I2) - eta2*S2 + nu*R2
         dE2 <-        eta1*E1 + delta*beta*S2*(I1+alpha*I2) - eta2*E2 - sigma*E2
         dI2 <-        eta1*I1                               - eta2*I2 + sigma*E2 - gamma*I2
         dR2 <-        eta1*R1                               - eta2*R2            + gamma*I2 - nu*R2

         # return (dS, dI, dR) as a vector in a list (required for the 'solve' function)
         return(list(c(dS1, dE1, dI1, dR1,dV, dS2, dE2, dI2, dR2)))
       }
  )
}
######################################### #
# SOLVE ODE                          ----
######################################### #
# use the 'ode' function of deSolve package with our SIR function, health states and parameters
out <- ode(func = seirv_func, y = states, times = times, parms = params)
plot(out)

# select final 8 epidemics
times_output <- seq(num_weeks-(52*8)-2,num_weeks,1)# burning time ## shifting peaks
out <- out[out[,1] %in% times_output,] # skip the initial time points
out[,1] <- out[,1] - min(out[,1])# rescale time points

# convert the 'out' matrix into a data-frame (to enable the use of '$' to access a column by name)
out <- as.data.frame(out)
names(out)

# calculate age-specific incidence
out$I1_incidence <- out$I1*pop_size
out$I2_incidence <- out$I2*pop_size

######################################### #
# PLOT RESULTS                       ----
######################################### #

# convert the 'out' matrix into a data-frame (to enable the use of '$' to access a column by name)
out <- as.data.frame(out)
names(out)

# plot results
par(mfrow =c(1,1))
plot(out$time,out$I2*pop_size,type='l',lwd=3, col = 4,main='model')
lines(out$time,out$I1*pop_size,type='l',lwd=3, col = 3)
legend('topleft',
       c('AG1','AG2'),
       col = c(3,4),
       lwd=2,
       cex=0.7)


######################################### #
# REFERENCE DATA                     ----
######################################### #

ref_data           <- read.csv("./RSV data/RSV_cases_time_epistat.csv")
ref_data$week      <- seq(1,dim(ref_data)[1])
ref_data_age       <- read.csv("./RSV data/RSV_cases_age_epistat.csv")
pro_case_less1yrs  <- ref_data_age$cases[ref_data_age$age==0] / sum(ref_data_age$cases)    #LW
ref_data$cases_ag1 <- round(ref_data$cases*pro_case_less1yrs)
ref_data$cases_ag2 <- ref_data$cases - ref_data$cases_ag1

# convert date into "date format"
ref_data$date      <- as.Date(ref_data$date)


# plot incidence by age
plot(ref_data$date,
     ref_data$cases_ag1,
     col = 3,lwd=2,
     type='b',ylab='incidence',main='Belgium')
lines(ref_data$date,
     ref_data$cases_ag2,
     type='b',
     col=4,lwd=2)
legend('topleft',
       c('AG1','AG2'),
       col = c(3,4),
       lwd=2,
       pch=1,
       cex=0.7)



######################################### #
# MODEL FITTING                        ----
######################################### #

# DEFINE HELP FUNCTION TO CALCULATE SUM OF SQUARES
get_sum_of_squares <- function(model_values,ref_values) {
  return(sum(sqrt((model_values-ref_values)^2)))
}

# SET GOAL: to use the optimisation package
library(Rcpp)
library(optimization)

# DEFINE FUNCTION TO RUN ODE WITH PARAMETER VECTOR X
x <- c(0.65, 0.65,0.65, 10) # delta[0.5,0.8], alpha[0.5,0.8], beta1[0,1]

get_model_output <- function(x){

  #TODO: update model parameters
  # params_fitting
  params_fit = params
  params_fit["delta"] = x[1]
  params_fit["alpha"] = x[2]
  params_fit["beta1"] = x[3]
  peak_shift          = x[4]
  
  # get output
  out <- data.frame(ode(func = seirv_func, 
                        y = states, 
                        times = seq(0,num_weeks,1), 
                        parms = params_fit))

  # apply peak shift
  out$time <- out$time + peak_shift
  
  # cutoff at 'num weeks'
  out <- out[out$time<=num_weeks,]

  # align model output and reference data
  model_selection <- seq(max(0,(nrow(out) - nrow(ref_data) +1)),nrow(out),1)
  out      <- out[model_selection,]
  out$time <- out$time - min(out$time)  # reset time

  
  # add age-specific incidence
  out$I1_incidence <- out$I1*pop_size
  out$I2_incidence <- out$I2*pop_size
  
  # return output
  return(out)
}


# HELP FUNCTION TO CALCULATE SUM OF SQUARES GIVEN PARAMETERS 'X'
get_parameter_score <- function(x) {

  # get model output given the parameters in "x"
  model_out <- get_model_output(x)

  # get model score
  model_incidence     <- c(model_out$I1_incidence, 
                           model_out$I2_incidence)
  reference_incidence <- c(ref_data$cases_ag1,ref_data$cases_ag2)
  model_score <- get_sum_of_squares(model_incidence,reference_incidence)

  # return model score
  return(model_score)
}

# # try some combinations
get_parameter_score(c(0.65,0.65,0.65,0))
get_parameter_score(c(0.65,0.65,0.65,15))

# HELP FUNCTION TO VISUALIZE THE MODEL FITTING
plot_model_fit <- function(x){

  # get model output given the parameters in "x"
  model_out <- get_model_output(x)

  # get model score
  model_score <- get_parameter_score(x)

  # plot reference data
  par(mfrow=c(2,1))
  y_lim <- range(ref_data$cases_ag1,ref_data$cases_ag2)*1.2
  plot(ref_data$week,
       ref_data$cases_ag1,
       col = 1,lwd=1,pch=20,
       ylim = y_lim,
       type='p',ylab='incidence',main='Belgium - AG1')
  lines(model_out$time,model_out$I1_incidence,col=3,lwd=2)
  
  plot(ref_data$week,
       ref_data$cases_ag2,
       col = 1,lwd=1,pch=20,
       ylim = y_lim,
       type='p',ylab='incidence',main='Belgium - AG2')
  lines(model_out$time,model_out$I2_incidence,col=4,lwd=2)
  
}

# try some combinations
plot_model_fit(c(0.65, 0.65, 0.65, 42))
plot_model_fit(c(0.65, 0.67, 0.62, 5))

#TODO: update model parameters
# try other function, by specifing lower and upper values
opt_param_sa <- optim_sa(fun = get_parameter_score, 
                         start = c(0.65, 0.65 , 0.65, 0),
                         trace = FALSE,
                         lower = c(0.5, 0.5, 0,   0),#
                         upper = c(0.8, 0.8, 1, 200),
                         control = list(dyn_rf = FALSE,
                                        rf = 1.2,
                                        t0 = 10, 
                                        nlimit = 500, 
                                        r = 0.6, 
                                        t_min = 0.1
                         ))$par
plot_model_fit(opt_param_sa)
opt_param_sa

# [1]    0.54    0.76    0.89 -153.37
plot_model_fit(c(0.54,0.76,0.89,153.37)) # delta[0.5,0.8], alpha[0.5,0.8], beta1[0,1]

