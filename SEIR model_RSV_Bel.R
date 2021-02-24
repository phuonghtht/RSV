# ########################## #
#SEIR model for RSV-------
########################## #
rm(list=ls())
setwd("/home/phuong/phuonght/gitKraken_project/RSV")
library(deSolve)
library(ggplot2)

######################################### #
#MODEL SETTINGS                     ----
######################################### #
pop_size          <- 11176549 # Belgium pop size (average 2011-2018)
ag1_size          <- pop_size*1/80
num_days          <- 500
num_weeks         <- 52*30
num_days_infected <- 10  #[ 8-11], 6.7 is original estimate
num_days_exposed  <- 4   # [2,6]
num_days_waning   <- 200 # [148, 164] best fit
infected_seeds    <- 10# ? 1= 1000

######################################### #
# INITIALIZE PARAMETERS AND POPULATION ----
######################################### #
#https://steemit.com/science/@fouad/x-history-of-math-symbols

#  population states

# sum both groups = 1
S1 <- 1/80 # less than 1 yrs old
E1 <- (infected_seeds/2)/(pop_size)
I1 <- 0
R1 <- 0

E2 <- (infected_seeds/2)/(pop_size)
I2 <- 0 # greater than 1 yrs old
R2 <- 0
S2 <- 1-S1-E1-I1-R1-E2-I2-R2

#
# #sum each group =1
# E1 <- infected_seeds/2/(pop_size)
# I1 <- 0
# R1 <- 0
# S1 <- 1-E1-I1-R1
#
# E2 <- infected_seeds/2/(pop_size)
# I2 <- 0
# R2 <- 0
# S2 <- 1-E2-I2-R2


######################################### #
# SET FUNCTION PARAMETERS            ----
######################################### #
# set time frame
times      <- seq(0, num_weeks, by = 1)

# set initial health states
states     <- c(S1 = S1,E1 = E1, I1 = I1, R1 = R1,
                S2 = S2,E2 = E2, I2 = I2, R2 = R2)

# set parameters
params     <- c(sigma = 1/0.57,              #rate of movement from latent to infectious stage
                gamma = 1/1.4,               # recovery rate
                nu=  7/num_days_waning,      # rate of loss of immunity
                eta1 = 1/(1*52),             # aging rate AG1
                eta2 = 1/(79*52),            # aging rate AG2
                # eta = 1/(80*52),           # death rate in general
                mu = 1/(80*52),              # birth rate
                beta1 = 0.65,                # the degree of seasonality, range [0,1],higher value stronger seasonal drivers
                phi = 2.43,                  # phase shift?
                beta0 = 1.99,                #average transmission rate
                # R0 = 3,                   # reproduction number
                delta = 0.65,                # scaled susceptibility
                # delta = 0.7,  # change delta make change the pattern
                alpha = 0.65              # scaled infectiousness
)

######################################### #
# CREATE SIRV FUNCTION               ----
######################################### #
# define a function for the ODE-solver with 3 function parameters (required)
#  * times  = current time point (not yet used here, in the absence of seasonality)
#  * states = current health states, to calculate updates
#  * params = epidemiogical parameters, used in the differential equations
#
# Note: updating the health states over time and keeping track of the changes, is handled by the
# ode-function of the 'deSolve' package.
sirv_func <- function(t, states, params) {

  with(as.list(c(states, params)),
       {

         # define variables

         # https://notebook.community/piklprado/ode_examples/Qualitative%20analysis%20and%20Bifurcation%20diagram%20Tutorial%20-%20R
         beta <- beta0*(1+beta1*cos(2*pi*t/52+phi))


         # calculate state changes
         dS1 <- mu - beta*S1*(I1+alpha*I2) - eta1*S1             + nu*R1
         dE1 <-      beta*S1*(I1+alpha*I2) - sigma*E1 - eta1*E1
         dI1 <-                              sigma*E1 - gamma*I1 - eta1*I1
         dR1 <-                                         gamma*I1 - nu*R1 - eta1*R1
         dS2 <- eta1*S1 - delta*beta*S2*(I1+alpha*I2) - eta2*S2             + nu*R2
         dE2 <- eta1*E1 + delta*beta*S2*(I1+alpha*I2) - sigma*E2 - eta2*E2
         dI2 <- eta1*I1                               + sigma*E2 - gamma*I2 - eta2*I2
         dR2 <- eta1*R1                                          + gamma*I2 - nu*R2  - eta2*R2


         # return (dS, dI, dR) as a vector in a list (required for the 'solve' function)
         return(list(c(dS1, dE1, dI1, dR1,dS2, dE2, dI2, dR2)))
       }
  )
}

# # plot the forcing seasonal function
# t = seq(1,200,1) # weekly
# beta0 = 1.99
# beta1 = 0.65 #[0,1] # the higiher beta1 the higher variation of beta = stronger seasonality driver
# phi= 1.5# peak position ~ week 40 , but Belgium peak at week 48 or 49 (1.5,2.5)
# y  <- beta0*(1+beta1*cos(2*pi*t/52+phi))
# plot(t,y,type="l", xlab="time", ylab="wave")

######################################### #
# SOLVE ODE                          ----
######################################### #
# use the 'ode' function of deSolve package with our SIR function, health states and parameters
out <- ode(func = sirv_func, y = states, times = times, parms = params)
plot(out)
# summary(out)
times_output <- seq(num_weeks-(52*8)-2,num_weeks,1)# burning time ## shifting peaks
times_output <- times_output[1:378] # take exactly 8 years (378 wks) for later scaling of data
out <- out[out[,1] %in% times_output,] # skip the initial time points
out[,1] <- out[,1] - min(out[,1])# rescale time points

par(mfrow = c(1,1))
matplot(out[,1], out[,c(4,8)], type = "l", xlab = "Time (weeks)", ylab = "population fraction")
legend("topright", col = 1:2, lty = 1:2, legend = c("I1","I2"))


######################################### #
# PLOT RESULTS                       ----
######################################### #

# convert the 'out' matrix into a data-frame (to enable the use of '$' to access a column by name)
out <- as.data.frame(out)
names(out)


# Total population
pop <- out[,"S1"] + out[,"E1"] + out[,"I1"] + out[,"R1"]+
  out[,"S2"] + out[,"E2"] + out[,"I2"] + out[,"R2"]
time <- out[,"time"]
plot(time,pop,type='l',lwd=3,ylim = c(0,2))# checking the consistence of pop

#BELGIUM DATA-----------------
case_dt <- read.csv("./RSV data/RSV_cases_time_epistat.csv")
case_dt$week      <- seq(1,dim(case_dt)[1])
age_dt <- read.csv("./RSV data/RSV_cases_age_epistat.csv")
# total cases: 63301,4yrs: 822,3yr:1794,2yr:4559 1yr:13751,0yr:42375
pro_case_less1yrs  <- age_dt$cases[age_dt$age==0] / sum(age_dt$cases)    #LW
case_dt$cases_ag1 <- round(case_dt$cases*pro_case_less1yrs)
case_dt$cases_ag2 <- case_dt$cases - case_dt$cases_ag1

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

# plot real data
par(mfrow = c(1,1))
plot(case_dt$week,
     case_dt$cases_ag1,type = "b",pch=19,
     ylim = c(0, max(case_dt$cases_ag1)*1.5))

# plot initial model
lines(out$time,out$I1*pop_size,col=2,lwd=2,type="b")

# score for initial model
get_sum_of_squares(out$I1[1:length(case_dt$cases)]*pop_size,case_dt$cases_ag1)

# DEFINE FUNCTION TO RUN ODE WITH PARAMETER VECTOR X
# x <- c(1.99, 0.65,-154, 1/200) #beta0[?] beta1(0,1],phi[1,2],shifting time, nu[1/160,1/230]#160-230 days

get_model_output <- function(x){

  # params_fitting
  params_fit = params
  params_fit["beta0"] = x[1]
  params_fit["beta1"] = x[2]
  time_shift =x[3]
  params_fit["nu"] = x[4]


  # get output
  out <- data.frame(ode(func = sirv_func, y = states, times = seq(0,num_weeks,1), parms = params_fit))

  # shift in time (fill with 0)
  out <- approx(x   = out$time +time_shift,
                y   = out$I1,
                xout = seq(0,num_weeks,1),
                rule = 2)
  names(out) <- c('time','I1')

  # rescale time points
  out$time <- out$time - min(out$time)

  # return output
  return(out)
}

# HELP FUNCTION TO CALCULATE SUM OF SQUARES GIVEN PARAMETERS 'X'
get_parameter_score <- function(x) {

  # get model output given the parameters in "x"
  model_out <- get_model_output(x)

  # get model score
  model_score <- get_sum_of_squares(model_out$I1[length(case_dt$cases)]*pop_size,case_dt$cases_ag1)

  # return model score
  return(model_score)
}


# HELP FUNCTION TO VISUALIZE THE MODEL FITTING
plot_model_fit <- function(x){

  # get model output given the parameters in "x"
  model_out <- get_model_output(x)

  # get model score
  model_score <- get_sum_of_squares(model_out$I1[1:length(case_dt$cases)]*pop_size,case_dt$cases_ag1)

  # plot reference data
  plot(case_dt$week,
       case_dt$cases_ag1,
       ylim = c(0,max(case_dt$cases_ag1)*2),type = "b",pch=19,
       main = paste('score:',round(model_score,digits = 2)))
  # plot initial model
  lines(model_out$time,model_out$I1*pop_size,col=2,lwd=1,type = "b")
}

# try other function, by specifing lower and upper values
# num_days_infected <- 10 #[ 8-11]
# num_days_exposed  <- 4  # [2,6]

opt_param_sa <- optim_sa(fun = get_parameter_score, start = c(1.9,0.65,-465,7/350),
                         trace = FALSE,
                         lower = c(1,0.1,-600,7/350),#
                         upper = c(3,1,-300,7/160),
                         control = list(dyn_rf = FALSE,
                                        rf = 1.2,
                                        t0 = 10, nlimit = 500, r = 0.6, t_min = 0.1
                         ))$par
opt_param_sa
plot_model_fit(opt_param_sa)

# [1]    2.16    0.75 -465.79    0.0

# some good fit?
plot_model_fit(c(1.78,0.63,-569,0.03))#33833.14
plot_model_fit(c(1.79,0.68,-674,0.035))#32986.37
plot_model_fit(c(1.8,0.7,-569,0.03))#31835.53
plot_model_fit(c(1.89,0.85,-465,0.015))#27885
plot_model_fit(c(1.81,0.74,-571,0.03)) #24836.87

