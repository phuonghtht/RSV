############################################################################ #
# EXAMPLE: SIR MODEL WITH ODE-SOLVER AND OPTIMISATION FUNCTION
#
#   --> 1. POPULATION MODEL WITH HOMOGENEOUS RANDOM MIXING
#   --> 2. PARAMETER ESTIMATION AND MODEL FITTING
#
# Author(s):    L. Willem, University of Antwerp, Belgium
# Packages:    'deSolve' (General Solver for Ordinary Differential Equations)
############################################################################ #

## set working directory (or open RStudio with this script)
# setwd("C:\\User\\path\\to\\the\\rcode\\folder") ## WINDOWS
# setwd("/Users/path/to/the/rcode/folder") ## MAC

# clear global environment
rm(list = ls())

# if the 'deSolve' package is not installed --> install
if(!'deSolve' %in% installed.packages()[,1]){ install.packages('deSolve')}

# load the 'deSolve' package
library(deSolve)


######################################### #
# MODEL SETTINGS                     ----
######################################### #
pop_size          <- 100000
num_days          <- 100
R0                <- 2
num_days_infected <- 7
infected_seeds    <- 100


######################################### #
# INITIALIZE PARAMETERS AND POPULATION ----
######################################### #
# recovery parameter
gamma  <- 1/num_days_infected

# transmission parameter
beta   <- R0*gamma

# population states
S      <- 1 - (infected_seeds/pop_size)
I      <- infected_seeds/pop_size
R      <- 0


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
sirv_func <- function(times, states, params) {

  # rename states and parameters
  S     <- states['S']
  I     <- states['I']
  beta  <- params['beta']
  gamma <- params['gamma']

  # calculate state changes
  dS <- -beta * S * I
  dI <-  beta * S * I - gamma * I
  dR <-                 gamma * I

  # return (dS, dI, dR) as a vector in a list (required for the 'solve' function)
  return(list(c(dS, dI, dR)))
}


######################################### #
# SET FUNCTION PARAMETERS            ----
######################################### #
# set time frame
times      <- seq(0, num_days, by = 1)

# set initial health states
states     <- c(S = S, I = I, R = R)

# set parameters
params     <- c(beta = beta, gamma = gamma)



######################################### #
# SOLVE ODE                          ----
######################################### #
# use the 'ode' function of deSolve package with our SIR function, health states and parameters
out <- ode(func = sirv_func, y = states, times = times, parms = params)

# convert the 'out' matrix into a data-frame (to enable the use of '$' to access a column by name)
out <- as.data.frame(out)


######################################### #
# PLOT RESULTS                       ----
######################################### #

# # open pdf-stream (optional)
# pdf('SIR_example.pdf',6,4) 

# plot susceptible class
plot(out$S,
     type = 'l',
     xlab = 'Time',
     ylab = 'Population fraction',
     ylim = c(0,1.15),
     col  = 4)
lines(out$R, col=3)
lines(out$I, col=2)

# add legend
legend(x      = 'top', 
       legend = c('Susceptible','Infectious','Recovered'), 
       col    = c(4,2,3), 
       lwd    = 1,
       ncol   = 3,
       cex    = 0.9)

# # close pdf-stream (optional)
# dev.off()


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

# example
# 1. create (dummy) function which takes 2 parameters and returns a score
hi <- function(x){get_sum_of_squares(0:10,seq(x[1],x[2],length=11))} 

# 2. use optimization with Nelder-Mead
round(optim_nm(fun = hi, k = 2)$par)


# SO... APPLY THIS TO ODE MODEL

# create artificial reference data
# example: increase incidence, add stochastic noise and shift peak in time
set.seed(2020) # set seed, to make the example reproducible
out_ref <- data.frame(ode(func = sirv_func, y = states, times = times, parms = c(beta = 0.38, gamma = 0.14)))
reference_data <- approx(x   = out_ref$time + 10,
                         y   = out_ref$I * runif(n = nrow(out_ref),min = 0.9,max=1.1),
                         xout = seq(0,num_days,1),
                         rule = 2)
names(reference_data) <- c('time','I')

# plot reference data
plot(reference_data$time,
       reference_data$I,
     ylim = c(0,0.5))

# plot initial model
lines(out$time,out$I,col=2,lwd=2)

# score for initial model
get_sum_of_squares(out$I,reference_data$I)


# DEFINE FUNCTION TO RUN ODE WITH PARAMETER VECTOR X
x <- c(0.4,0.2,0)
get_model_output <- function(x){
  
  # get output
  out <- data.frame(ode(func = sirv_func, y = states, times = seq(0,num_days,1), parms = c(beta=x[1],gamma=x[2])))
  
  # shift in time (fill with 0)
  out <- approx(x   = out$time + x[3],
                y   = out$I,
                xout = seq(0,num_days,1),
                rule = 2)
  names(out) <- c('time','I')
  
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
  model_score <- get_sum_of_squares(model_out$I,reference_data$I)
  
  # return model score
  return(model_score)
}

# try some combinations
get_parameter_score(c(0.4,0.2,0))
get_parameter_score(c(0.38, 0.14,15))
get_parameter_score(c(0.38, 0.14,-6))

# HELP FUNCTION TO VISUALIZE THE MODEL FITTING
plot_model_fit <- function(x){
  
  # get model output given the parameters in "x"
  model_out <- get_model_output(x)
  
  # get model score
  model_score <- get_sum_of_squares(model_out$I,reference_data$I)
  
  # plot reference data
  plot(reference_data$time,
       reference_data$I,
       ylim = c(0,0.5),
       main = paste('score:',round(model_score,digits = 2)))
  
  # plot initial model
  lines(model_out$time,model_out$I,col=2,lwd=2)
  
}

# try some combinations
plot_model_fit(c(0.4,0.2,0))
plot_model_fit(c(0.38, 0.14,15))
plot_model_fit(c(0.38, 0.14,-6))


# Nelder-Mead optimisation
opt_param <- optim_nm(fun = get_parameter_score, start = c(0.1,0.1,1))$par
plot_model_fit(opt_param)

# initial values have an effect...
opt_param <- optim_nm(fun = get_parameter_score, start = c(0.1,0.1,-1))$par
plot_model_fit(opt_param)

# try other function, by specifing lower and upper values
opt_param_sa <- optim_sa(fun = get_parameter_score, start = c(0.1,0.1,-1),
                         trace = FALSE, 
                         lower = c(0.01, 0.01,-10),
                         upper = c(0.8, 0.8,20),
                         control = list(dyn_rf = FALSE,
                                        rf = 1.2,
                                        t0 = 10, nlimit = 500, r = 0.6, t_min = 0.1
                         ))$par
plot_model_fit(opt_param_sa)

# to be continued...
