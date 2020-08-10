########################## #
#SEIR model for RSV-------
########################## #

rm(list = ls())
setwd("/home/phuong/phuonght/gitKraken_project/RSV")
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyverse)

######################################### #
# MODEL SETTINGS                     ----
######################################### #
pop_size          <- 11000000 # https://statbel.fgov.be/en, 11492641, growth rate 0.54%
num_days          <- 500
num_weeks         <- 52*20
num_days_infected <- 10  #[ 8-11], 6.7 is original estimate
num_days_exposed  <- 4   # [2,6]
num_days_waning   <- 160  # [148, 164] best fit 
infected_seeds    <- 10# 

######################################### #
# INITIALIZE PARAMETERS AND POPULATION ----
######################################### #
#https://steemit.com/science/@fouad/x-history-of-math-symbols

#  population states

S1 <- 1/80 # less than 1 yrs old
E1 <- infected_seeds/(pop_size)
I1 <- 0
R1 <- 0

E2 <- 0 #  [0,2)
I2 <- 0
R2 <- 0
S2 <- 0

E3 <- 0 # > 2
I3 <- 0
R3 <- 0
S3 <- 1-S1-E1-I1-R1-S2-E2-I2-R2-E3-I3-R3


######################################### #
# SET FUNCTION PARAMETERS            ----
######################################### #
# set time frame
times      <- seq(0, num_weeks, by = 1)

# set initial health states
states     <- c(S1 = S1,E1 = E1, I1 = I1, R1 = R1,
                S2 = S2,E2 = E2, I2 = I2, R2 = R2,
                S3 = S3,E3 = E3, I3 = I3, R3 = R3)

# set parameters
params     <- c(sigma = 7/num_days_exposed, #rate of movement from latent to infectious stage
                gamma = 7/num_days_infected, # recovery rate
                # sigma = x[1], #rate of movement from latent to infectious stage
                # gamma = x[2], # recovery rate
                # nu= 0.044,     # rate of loss of immunity
                nu=  7/ num_days_waning,    # rate of loss of immunity
                eta1 = 1/(1*52),            # aging rate 1 yr old 
                eta2 = 1/(1*52),            # aging rate 2 yrs old
                eta3 = 1/(78*52),           # aging rate > 2 yrs old
                # eta = 1/(80*52),             # death rate in general
                mu = 1/(80*52),             # birth rate 1.66 2020
                beta1 = 0.65,               # the degree of seasonality, range [0,1],higher value stronger seasonal drivers
                phi = 2.43,                 # phase shift?
                beta0 = 2.05,               # average transmission rate# this
                
                # # for 3 classes of age
                # beta1 = 0.522,               # the degree of seasonality, range [0,1],higher value stronger seasonal drivers
                # phi = 2.43,                 # phase shift?
                # beta0 = 530,               #average transmission
                
                # R0 = 3,                   # reproduction number
                delta2 = 0.228,             # scaled susceptibility
                delta3 = 0.6,
                alpha2 = 1,                 # scaled infectiousness
                alpha3 = 0.6                # scaled infectiousness
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
         # P <- (S1+E1+I1+R1+S2+E2+I2+R2)
         beta <- beta0*(1+beta1*cos(2*pi*t/52+phi))
         # beta <- beta0*(1+beta1*sin(2*pi*t/52))#
         
         
         # calculate state changes
         dS1 <- mu - beta*S1*(I1+alpha2*I2+alpha3*I3) - eta1*S1 + nu*R1
         dE1 <- beta*S1*(I1+alpha2*I2+alpha3*I3) - eta1*E1 - sigma*E1 
         dI1 <- sigma*E1 - eta1*I1 - gamma*I1
         dR1 <- gamma*I1 - eta1*R1 - nu*R1
         
         dS2 <- eta1*S1 - delta2*beta*S2*(I1 +alpha2*I2 + alpha3*I3) - eta2*S2 + nu*R2
         dE2 <- eta1*E1 + delta2*beta*S2*(I1 +alpha2*I2 + alpha3*I3) - eta2*E2 - sigma*E2
         dI2 <- eta1*I1 + sigma*E2 - eta2*I2 - gamma*I2
         dR2 <- eta1*R1 + gamma*I2 - eta2*R2 - nu*R2
         
         dS3 <- eta2*S2 - delta3*beta*S3*(I1 + alpha2*I2 + alpha3*I3) - eta3*S3 + nu*R3
         dE3 <- eta2*E2 + delta3*beta*S3*(I1 + alpha2*I2 + alpha3*I3) - eta3*E3 - sigma*E3
         dI3 <- eta2*I2 + sigma*E3 - eta3*I3 - gamma*I3
         dR3 <- eta2*R2 + gamma*I3 - eta3*R3 - nu*R3
         
         # return (dS, dI, dR) as a vector in a list (required for the 'solve' function)
         return(list(c(dS1, dE1, dI1, dR1,dS2, dE2, dI2, dR2, dS3, dE3, dI3, dR3)))
       }
  )
}

######################################### #
# SOLVE ODE                          ----
######################################### #
# use the 'ode' function of deSolve package with our SIR function, health states and parameters
out <- ode(func = sirv_func, y = states, times = times, parms = params)
# plot(out)
# out <- as.data.frame(out)
# par(mfrow=c(1,1))
# plot(out$time, out$I1)

# summary(out)
times_output <- seq(num_weeks-(52*8)-2,num_weeks,1)# burning time ## shifting peaks
out <- out[out[,1] %in% times_output,] # skip the initial time points
out[,1] <- out[,1] - min(out[,1])# rescale time points

# par(mfrow = c(1,1))
# matplot(out[,1], out[,2:9], type = "l", xlab = "time", ylab = "population fraction")
# legend("topright", col = 1:8, lty = 1:8, legend = c("S", "E", "I","R","S2", "E2", "I2","R2"))


######################################### #
# PLOT RESULTS                       ----
######################################### #

# convert the 'out' matrix into a data-frame (to enable the use of '$' to access a column by name)
out <- as.data.frame(out)
# out$S <- out[,"S1"] + out[,"S2"] + out[,"S3"]
# out$E <- out[,"E1"] + out[,"E2"] + out[,"E3"]
# out$I <- out[,"I1"] + out[,"I2"] + out[,"I3"]
# out$R <- out[,"R1"] + out[,"R2"] + out[,"R3"]

# Total population
# pop <- out[,"S"] + out[,"E"] + out[,"I"] + out[,"R"]

# weekly incidence
inc <- out[,"I1"]
inc2 <- out[,"I2"]
inc3 <- out[,"I3"]

# inc <- inc*pop_size
# inc2<- inc2*pop_size

time <- out[,"time"]
# plot(time,pop,type='l',lwd=3,ylim = c(0,2))# checking the consistence of pop

# Plo"infectious" proportion I
ggplot(out%>%gather(key = "I", value = "value", starts_with("I")),
       aes(x=time, y = value, colour = I))+
  geom_point()

ggplot(out, aes(time))+
  geom_line(aes(y = I1), color = "darkred")+
  geom_line(aes(y = I2), color = "steelblue")+
  ylab("")
  

# par(mfrow = c(1,1))
# matplot(out[,1], out[,2:5], type = "l", xlab = "time", ylab = "population fraction")
# legend("topright", col = 1:4, lty = 1:4, legend = c(colnames(out)[2:5]))

# # plot susceptible class
# par(mfrow = c(1,1))
# plot(out$S,
#      type = 'l',
#      xlab = 'Time',
#      ylab = 'Population fraction',
#      ylim = c(0,1.3),
#      col  = 4)
# lines(out$R, col=3)
# lines(out$I, col=2)
# lines(out$E, col=1)
#
# # add legend
# legend(x      = 'top',
#        legend = c('Susceptible','Infectious','Recovered',"Exposed"),
#        col    = c(4,2,3,1),
#        lwd    = 1,
#        ncol   = 3,
#        cex    = 0.9)
# 

#BELGIUM DATA-----------------
bel_data <- read.csv("./RSV data/RSV_cases_time_epistat.csv")
# total cases: 63301 , less than 2yrs: 56126, 1yr:42375,2yrs:13751
pro_case_less2yrs <- 56126/63301
pro_case_2yrs <- 13751/63301
pro_case_1yr <- 42375/63301

bel_data$week <- seq(1,dim(bel_data)[1])
# bel_data$cases.less2yrs <- round(bel_data$cases*pro_case_less2yrs)
bel_data$cases.1yr <- round(bel_data$cases*pro_case_1yr)
bel_data$cases.2yrs <- round(bel_data$cases*pro_case_2yrs)

## score for initial model
 # ggplot(bel_data, aes(week, cases.less2yrs))+
#   geom_point()


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

# # 1. create (dummy) function which takes 2 parameters and returns a score
# hi <- function(x){get_sum_of_squares(0:10,seq(x[1],x[2],length=11))} 
# 
# # 2. use optimization with Nelder-Mead
# round(optim_nm(fun = hi, k = 2)$par)# k= 2 parameters


# plot real data
par(mfrow = c(1,1))
plot(bel_data$week,
     bel_data$cases.1yr,
     ylim = c(0, max(bel_data$cases.1yr)*1.5))
# plot(bel_data$week,
#      bel_data$cases.2yrs,
#      ylim = c(0, max(bel_data$cases.2yrs)*1.5))

# plot initial model
lines(out$time,out$I1*pop_size,col=2,lwd=2)
# plot(out$time,out$I1*pop_size,col=2,lwd=2)

# score for initial model
get_sum_of_squares(out$I1[1:length(bel_data$cases)]*pop_size,bel_data$cases.1yr)

# DEFINE FUNCTION TO RUN ODE WITH PARAMETER VECTOR X
x <- c(2.05,0.65,0) #beta0, beta1,

get_model_output <- function(x){
  
   # params_fitting
  params_fit = params
  params_fit["beta0"] = x[1]
  params_fit["beta1"] = x[2]
  # params_fit["beta1"] = x[3]
  
  
   # get output
  out <- data.frame(ode(func = sirv_func, y = states, times = seq(0,num_weeks,1), parms = params_fit))
  
  # shift in time (fill with 0)
  out <- approx(x   = out$time + x[3],
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
  model_score <- get_sum_of_squares(model_out$I1[1:length(bel_data$cases)]*pop_size,bel_data$cases.1yr)
  
  # return model score
  return(model_score)
}

# try some combinations
get_parameter_score(c(2.05,0.65,0))
get_parameter_score(c(2.01,0.65,-52*11.96))


# HELP FUNCTION TO VISUALIZE THE MODEL FITTING
plot_model_fit <- function(x){
  
  # get model output given the parameters in "x"
  model_out <- get_model_output(x)
  
  # get model score
  model_score <- get_sum_of_squares(model_out$I1[1:length(bel_data$cases)]*pop_size,bel_data$cases.1yr)
  
  # plot reference data
  plot(bel_data$week,
       bel_data$cases.1yr,
       ylim = c(0,max(bel_data$cases.1yr)*2),
       main = paste('score:',round(model_score,digits = 2)))
  
  # plot initial model
  lines(model_out$time,model_out$I1*pop_size,col=2,lwd=2)
  # plot(model_out$time,model_out$I1*pop_size,col=2,lwd=2)
  
}

# try some combinations
plot_model_fit(c(2.05,0.65,0))
plot_model_fit(c(2.05,0.65,-52*11.96))


# 
# # Nelder-Mead optimisation
# opt_param <- optim_nm(fun = get_parameter_score, start = c(7/6,7/11,-154))$par
# plot_model_fit(opt_param)
# 
# # initial values have an effect...
# opt_param <- optim_nm(fun = get_parameter_score, start = c(7/6,7/11,-102))$par
# plot_model_fit(opt_param)

# try other function, by specifing lower and upper values
# num_days_infected <- 10 #[ 8-11]
# num_days_exposed  <- 4  # [2,6]
opt_param_sa <- optim_sa(fun = get_parameter_score, start = c(1.6, 0.4,-52*11.96),
                         trace = F, 
                         lower = c(1.6,0.4,-52*11.96),
                         upper = c(2.2,1,52*11.96),
                         control = list(dyn_rf = T,
                                        rf = 1.2,
                                        t0 = 10, nlimit = 500, r = 0.6, t_min = 0.1
                         ))$par
opt_param_sa
get_parameter_score(opt_param_sa)
plot_model_fit(opt_param_sa)


