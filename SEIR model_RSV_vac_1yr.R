########################## #
#I.SEIR model for RSV-------
########################## #
rm(list=ls())
setwd("/home/phuong/phuonght/gitKraken_project/RSV")
library(deSolve)
library(ggplot2)

######################################### #
#MODEL SETTINGS                     ----
######################################### #
pop_size          <- 1100000 # Belgium pop size
num_days          <- 500
num_weeks         <- 52*20
num_days_infected <- 10  #[ 8-11], 6.7 is original estimate
num_days_exposed  <- 4   # [2,6]
num_days_waning   <- 200 # [148, 164] best fit 
infected_seeds    <- 10# ? 1= 1000

######################################### #
# INITIALIZE PARAMETERS AND POPULATION ----
######################################### #
#https://steemit.com/science/@fouad/x-history-of-math-symbols

#  population states
# less than 1 yrs old:13759

# sum both groups = 1
S1 <- 1/80 # less than 1 yrs old
E1 <- infected_seeds/(pop_size)
I1 <- infected_seeds/(pop_size)
R1 <- 0

E2 <- 0
I2 <- (infected_seeds-3)/(pop_size) # greater than 1 yrs old
R2 <- 0
S2 <- 1-S1-E1-I1-R1-E2-I2-R2


# #sum each group =1
# E1 <- 0 # <1 yr
# I1 <- infected_seeds/(pop_size*(1/80))
# R1 <- 0
# S1 <- 1-E1-I1-R1
# 
# E2 <- 0 # 1-2 yrs
# I2 <- infected_seeds/(pop_size*(1/80))
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
params     <- c(sigma = 1/0.57, #rate of movement from latent to infectious stage
                gamma = 1/1.4, # recovery rate
                # sigma = x[1], #rate of movement from latent to infectious stage
                # gamma = x[2], # recovery rate
                # nu= 0.044,     # rate of loss of immunity
                nu=  7/ num_days_waning,     # rate of loss of immunity
                eta1 = 1/(1*52),             # aging rate
                eta2 = 1/(79*52),            # aging rate
                # eta = 1/(80*52),             # death rate in general
                mu = 1/(80*52),              # birth rate
                # mu = 0.0135/52,            # birth rate Perth Australia
                beta1 = 0.65,                # the degree of seasonality, range [0,1],higher value stronger seasonal drivers
                phi = 2.43,                  # phase shift?
                beta0 = 1.99,                #average transmission rate
                # R0 = 3,                     # reproduction number
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
         # P <- (S1+E1+I1+R1+S2+E2+I2+R2)
         # S1 <- 1- E1-I1-R1
         # S2 <- 1- E2-I2-R2
         # https://notebook.community/piklprado/ode_examples/Qualitative%20analysis%20and%20Bifurcation%20diagram%20Tutorial%20-%20R
         beta <- beta0*(1+beta1*cos(2*pi*t/52+phi))
         # beta <- beta0*(1+beta1*sin(2*pi*t/52+phi))
         # beta <- beta0*(1+beta1*u)

         # calculate state changes
         dS1 <- mu - beta*S1*(I1+alpha*I2) - eta1*S1 + nu*R1
         dE1 <- beta*S1*(I1+alpha*I2) - eta1*E1 - sigma*E1 
         dI1 <- sigma*E1 - eta1*I1 - gamma*I1
         dR1 <- gamma*I1 - eta1*R1 - nu*R1
         
         dS2 <- eta1*S1 - delta*beta*S2*(I1+alpha*I2) - eta2*S2 + nu*R2
         dE2 <- eta1*E1 + delta*beta*S2*(I1+alpha*I2) - eta2*E2 - sigma*E2
         dI2 <- eta1*I1 + sigma*E2 - eta2*I2 - gamma*I2
         dR2 <- eta1*R1 + gamma*I2 - eta2*R2 - nu*R2
         
         
         # u = u*(1-u^2-v^2)-2*pi*v
         # v= v*(1-u^2-v^2)+2*pi*u
         
         # return (dS, dI, dR) as a vector in a list (required for the 'solve' function)
         return(list(c(dS1, dE1, dI1, dR1,dS2, dE2, dI2, dR2)))
       }
  )
}

######################################### #
# SOLVE ODE                          ----
######################################### #
# use the 'ode' function of deSolve package with our SIR function, health states and parameters
out <- ode(func = sirv_func, y = states, times = times, parms = params)
plot(out)
# summary(out)
times_output <- seq(num_weeks-(52*8)-2,num_weeks,1)# burning time ## shifting peaks
out <- out[out[,1] %in% times_output,] # skip the initial time points
out[,1] <- out[,1] - min(out[,1])# rescale time points

# par(mfrow = c(1,1))
# matplot(out[,1], out[,c(4,8)], type = "l", xlab = "Time (weeks)", ylab = "population fraction")
# legend("topright", col = 1:2, lty = 1:2, legend = c("I1","I2"))


######################################### #
# PLOT RESULTS                       ----
######################################### #

# convert the 'out' matrix into a data-frame (to enable the use of '$' to access a column by name)
out <- as.data.frame(out)
out$S <- out[,"S1"] + out[,"S2"]
out$E <- out[,"E1"] + out[,"E2"]
out$I <- out[,"I1"] + out[,"I2"]
out$R <- out[,"R1"] + out[,"R2"]

# Total population
pop <- out[,"S1"] + out[,"E1"] + out[,"I1"] + out[,"R1"]+out[,"S2"] + out[,"E2"] + out[,"I2"] + out[,"R2"]

# # weekly incidence
# inc1 <- out[,"I1"]
# inc2 <- out[,"I2"]
# inc1 <- inc*pop_size
# inc2<- inc1*pop_size
# 
# time <- out[,"time"]
# plot(time,pop,type='l',lwd=3,ylim = c(0,2))# checking the consistence of pop
# 
# # Plot I1 (<1yrs) and I2(>1yrs)
# par(mfrow=c(1,2))
# plot(time,inc1,type='l',lwd=3, col = 1)
# plot(time,inc2,type='l',lwd=3, col = 1)


#BELGIUM DATA-----------------
bel_data <- read.csv("./RSV data/RSV_cases_time_epistat.csv")
# total cases: 63301 , less than 2yrs: 56126, 1yr:13751
pro_case_less1yrs <- 13751/63301
bel_data$week <- seq(1,dim(bel_data)[1])
bel_data$cases.less1yrs <- round(bel_data$cases*pro_case_less1yrs)
## score for initial model
# ggplot(bel_data, aes(week, cases.less1yrs))+
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

# plot real data
par(mfrow = c(1,1))
plot(bel_data$week,
     bel_data$cases.less1yrs,
     ylim = c(0, max(bel_data$cases.less1yrs)*1.5))

# plot initial model
lines(out$time,out$I1*pop_size,col=2,lwd=2)

# score for initial model
get_sum_of_squares(out$I1[1:length(bel_data$cases)]*pop_size,bel_data$cases.less1yrs)

# DEFINE FUNCTION TO RUN ODE WITH PARAMETER VECTOR X
x <- c(0.65, 0.65,0.65, -154) # delta[0.5,0.8], alpha[0.5,0.8], beta1[0,1]

get_model_output <- function(x){

  # params_fitting
  params_fit = params
  params_fit["delta"] = x[1]
  params_fit["alpha"] = x[2]
  params_fit["beta1"] = x[3]

  # get output
  out <- data.frame(ode(func = sirv_func, y = states, times = seq(0,num_weeks,1), parms = params_fit))

  # shift in time (fill with 0)
  out <- approx(x   = out$time + x[4],
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
  model_score <- get_sum_of_squares(model_out$I1[length(bel_data$cases)]*pop_size,bel_data$cases.less1yrs)

  # return model score
  return(model_score)
}

# # try some combinations
# get_parameter_score(c(0.65,0.65,0.65,0))
# get_parameter_score(c(0.65,0.65,0.65,-154))

# HELP FUNCTION TO VISUALIZE THE MODEL FITTING
plot_model_fit <- function(x){

  # get model output given the parameters in "x"
  model_out <- get_model_output(x)

  # get model score
  model_score <- get_sum_of_squares(model_out$I1[1:length(bel_data$cases)]*pop_size,bel_data$cases.less1yrs)

  # plot reference data
  plot(bel_data$week,
       bel_data$cases.less1yrs,
       ylim = c(0,max(bel_data$cases.less1yrs)*2),
       main = paste('score:',round(model_score,digits = 2)))

  # plot initial model
  lines(model_out$time,model_out$I1*pop_size,col=2,lwd=2)

}

# # try some combinations
# plot_model_fit(c(0.65,0.65,0.65,0))
# plot_model_fit(c(0.6,0.6,0.65,-154))

# try other function, by specifing lower and upper values
# num_days_infected <- 10 #[ 8-11]
# num_days_exposed  <- 4  # [2,6]
opt_param_sa <- optim_sa(fun = get_parameter_score, start = c(0.6, 0.6,0.8,-200),
                         trace = FALSE,
                         lower = c(0.5, 0.5,0,-519),#
                         upper = c(0.8, 0.8,1,-100),
                         control = list(dyn_rf = FALSE,
                                        rf = 1.2,
                                        t0 = 10, nlimit = 500, r = 0.6, t_min = 0.1
                         ))$par
plot_model_fit(opt_param_sa)
opt_param_sa
# [1]    0.54    0.76    0.89 -149.08
plot_model_fit(c(0.54,0.76,0.89,-153.37)) # delta[0.5,0.8], alpha[0.5,0.8], beta1[0,1]

########################################## #
# II. VACCINATION                   -----------
########################################## #

rm(list=ls())
setwd("/home/phuong/phuonght/gitKraken_project/RSV")
library(deSolve)
library(ggplot2)

#Changeable variables  ##########################
vac.frac = 1

######################################### #
# MODEL SETTINGS                     ----
######################################### #
pop_size          <- 1100000 # Belgium pop size
num_days          <- 500
num_weeks         <- 52*20
num_days_infected <- 10  #[ 8-11], 6.7 is original estimate
num_days_exposed  <- 4   # [2,6]
num_days_waning   <- 200 # [148, 164] best fit
infected_seeds    <- 10# ? 1= 1000

######################################### #
# INITIALIZE PARAMETERS AND POPULATION ----
######################################### #
#https://steemit.com/science/@fouad/x-history-of-math-symbols

#  population states
# less than 1 yrs old:13759

S1 <- 1/80 # less than 1 yrs old
E1 <- infected_seeds/(pop_size)
I1 <- 0
R1 <- 0
V  <- 0
E2 <- (infected_seeds-3)/(pop_size) # greater than 1 yrs old
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
                mu = 1/(80*52),              # birth rate
                beta1 = 0.89,                # the degree of seasonality, range [0,1],higher value stronger seasonal drivers
                phi = 2.43,                  # phase shift?
                beta0 = 1.99,                #average transmission rate
                # R0 = 3,                    # reproduction number
                delta = 0.54,                # scaled susceptibility
                alpha = 0.76,                # scaled infectiousness
                p = vac.frac               # vaccination proportion
                )

######################################### #
# CREATE SIRV FUNCTION               ----
######################################### #
# Note: updating the health states over time and keeping track of the changes, is handled by the
# ode-function of the 'deSolve' package.
sirv_func <- function(t, states, params) {

  with(as.list(c(states, params)),
       {
         # define variables
         # P <- (S1+E1+I1+R1+S2+E2+I2+R2)
         beta <- beta0*(1+beta1*cos(2*pi*t/52+phi))


         # calculate state changes
         dS1 <- (1-p)* mu - beta*S1*(I1+alpha*I2) - eta1*S1 + nu*R1
         dE1 <- beta*S1*(I1+alpha*I2) - eta1*E1 - sigma*E1
         dI1 <- sigma*E1 - eta1*I1 - gamma*I1
         dR1 <- gamma*I1 - eta1*R1 - nu*R1

         dV <- p*mu - nu*V

         dS2 <- nu*V + eta1*S1 - delta*beta*S2*(I1+alpha*I2) - eta2*S2 + nu*R2
         dE2 <- eta1*E1 + delta*beta*S2*(I1+alpha*I2) - eta2*E2 - sigma*E2
         dI2 <- eta1*I1 + sigma*E2 - eta2*I2 - gamma*I2
         dR2 <- eta1*R1 + gamma*I2 - eta2*R2 - nu*R2

         # return (dS, dI, dR) as a vector in a list (required for the 'solve' function)
         return(list(c(dS1, dE1, dI1, dR1,dV, dS2, dE2, dI2, dR2)))
       }
  )
}
######################################### #
# SOLVE ODE                          ----
######################################### #
# use the 'ode' function of deSolve package with our SIR function, health states and parameters
out <- ode(func = sirv_func, y = states, times = times, parms = params)
plot(out)
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
out$S <- out[,"S1"] + out[,"S2"]
out$E <- out[,"E1"] + out[,"E2"]
out$I <- out[,"I1"] + out[,"I2"]
out$R <- out[,"R1"] + out[,"R2"]

# Total population
pop <- out[,"S1"] + out[,"E1"] + out[,"I1"] + out[,"R1"]+
  out[,"S2"] + out[,"E2"] + out[,"I2"] + out[,"R2"] + out[,"V"]

time <- out[,"time"]
plot(time,pop,type='l',lwd=3,ylim = c(0,2))# checking the consistence of pop


par(mfrow =c(1,1))
plot(time,out$I1*pop_size,type='l',lwd=3, col = 1)
# plot(time,out$I*pop_size,type='l',lwd=3, col = 1)

#BELGIUM DATA-----------------
bel_data <- read.csv("./RSV data/RSV_cases_time_epistat.csv")
# total cases: 63301 , less than 2yrs: 56126, 1yr:13751
pro_case_less1yrs <- 13751/63301
bel_data$week <- seq(1,dim(bel_data)[1])
bel_data$cases.less1yrs <- round(bel_data$cases*pro_case_less1yrs)
## score for initial model
# ggplot(bel_data, aes(week, cases.less1yrs))+
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

# plot real data
par(mfrow = c(1,1))
plot(bel_data$week,
     bel_data$cases.less1yrs,
     xlab = "Time (weeks)",
     ylab = "Number of RSV reported",
     main= paste("Vaccination =",vac.frac),
     ylim = c(0, max(bel_data$cases.less1yrs)*1.5))

# plot initial model
lines(out$time,out$I1*pop_size,col=2,lwd=2)

# score for initial model
get_sum_of_squares(out$I1[1:length(bel_data$cases)]*pop_size,bel_data$cases.less1yrs)

# DEFINE FUNCTION TO RUN ODE WITH PARAMETER VECTOR X
x <- c(0.65, 0.65,0.65, -154) # delta[0.5,0.8], alpha[0.5,0.8], beta1[0,1]

get_model_output <- function(x){

  # params_fitting
  params_fit = params
  params_fit["delta"] = x[1]
  params_fit["alpha"] = x[2]
  params_fit["beta1"] = x[3]

  # get output
  out <- data.frame(ode(func = sirv_func, y = states, times = seq(0,num_weeks,1), parms = params_fit))

  # shift in time (fill with 0)
  out <- approx(x   = out$time + x[4],
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
  model_score <- get_sum_of_squares(model_out$I1[length(bel_data$cases)]*pop_size,bel_data$cases.less1yrs)

  # return model score
  return(model_score)
}

# # try some combinations
# get_parameter_score(c(0.65,0.65,0.65,0))
# get_parameter_score(c(0.65,0.65,0.65,-154))

# HELP FUNCTION TO VISUALIZE THE MODEL FITTING
plot_model_fit <- function(x){

  # get model output given the parameters in "x"
  model_out <- get_model_output(x)

  # get model score
  model_score <- get_sum_of_squares(model_out$I1[1:length(bel_data$cases)]*pop_size,bel_data$cases.less1yrs)

  # # plot reference data
  # plot(bel_data$week,
  #      bel_data$cases.less1yrs,
  #      ylim = c(0,max(bel_data$cases.less1yrs)*2),
  #      main = paste('score:',round(model_score,digits = 2)))
  #
  # # plot initial model
  # lines(model_out$time,model_out$I1*pop_size,col=2,lwd=2)

}

# # try some combinations
# plot_model_fit(c(0.65,0.65,0.65,0))
# plot_model_fit(c(0.6,0.6,0.65,-154))

# try other function, by specifing lower and upper values
# num_days_infected <- 10 #[ 8-11]
# num_days_exposed  <- 4  # [2,6]
# opt_param_sa <- optim_sa(fun = get_parameter_score, start = c(0.6, 0.6,0.8,-200),
#                          trace = FALSE,
#                          lower = c(0.5, 0.5,0,-519),#
#                          upper = c(0.8, 0.8,1,-100),
#                          control = list(dyn_rf = FALSE,
#                                         rf = 1.2,
#                                         t0 = 10, nlimit = 500, r = 0.6, t_min = 0.1
#                          ))$par
# plot_model_fit(opt_param_sa)
opt_param_sa
# [1]    0.54    0.76    0.89 -153.37
plot_model_fit(c(0.54,0.76,0.89,-153.37)) # delta[0.5,0.8], alpha[0.5,0.8], beta1[0,1]

