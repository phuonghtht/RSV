########################## #
#SEIR model for RSV-------
########################## #

setwd("/home/phuong/phuonght/gitKraken_project/RSV")

library(deSolve)

######################################### #
# MODEL SETTINGS                     ----
######################################### #

##  latency period: 4 days [2-6]
##  infectious period : 6.7 days [ 1-21], 10 days [7-12]
##  duration of immunity: 230 days / 200 days

######################################### #
# MODEL SETTINGS                     ----
######################################### #
pop_size          <- 100000
num_days          <- 500
num_weeks         <- 52*5
R0                <- 3
num_days_infected <- 6.7
num_days_exposed  <- 4
num_days_waning   <- 230 # best fit : 23.5 weeks => 164.5days
num_weeks_waning  <- 23.5
infected_seeds    <- 1 # ? 

######################################### #
# INITIALIZE PARAMETERS AND POPULATION ----
######################################### #
#https://steemit.com/science/@fouad/x-history-of-math-symbols
# recovery parameter
gamma  <- 1/(num_days_infected/7)
# transmission parameter
beta0 <- 1.99 #average transmission rate
beta1 <- 0.65 #degree of seasonality [ 0,1], higher value stronger seasonal driver
phi <- 2.43 # shift phase
# beta <- beta0*(1+beta1*cos(2*pi*t/52)+phi)## t need to be updated 

#rate of movement from latent to infectious stage
sigma <- 1/(num_days_exposed/7)
# duration of immunity
nu <- 1/(num_days_waning/7)
# death rate 
eta <- 1/(80*52)
# birth rate
mu <- 1/(80*52)

# # population states
S <- 1 - (infected_seeds/pop_size)
E <- infected_seeds/pop_size
I <- 0
R <- 0

######################################### #
# SET FUNCTION PARAMETERS            ----
######################################### #
# set time frame
times      <- seq(0, num_weeks, by = 1)

# set initial health states
states     <- c(S = S,E = E, I = I, R = R)

# set parameters
params     <- c(gamma = gamma, nu = nu,mu = mu,
                beta0 = beta0, beta1 = beta1, phi = phi)

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

  # rename states and parameters
  S     <- states['S']
  E     <- states['E']
  I     <- states['I']
  R     <- states['R']
  beta <- beta0*(1+beta1*cos(2*pi*t/52)+phi)
  # gamma <- params['gamma']
  # nui   <- params['nui']
  # v   <- params['v']
  


  # calculate state changes
  dS <- mu -beta*S*I + nu*R - eta*S
  dE <- beta*S*I - sigma*E - eta*E
  dI <- sigma*E - gamma*I - eta*I
  dR <- gamma*I - nu*R - eta*R

  # return (dS, dI, dR) as a vector in a list (required for the 'solve' function)
  return(list(c(dS, dE, dI, dR)))
}

######################################### #
# SOLVE ODE                          ----
######################################### #
# use the 'ode' function of deSolve package with our SIR function, health states and parameters
out <- ode(func = sirv_func, y = states, times = times, parms = params)
plot(out)

# convert the 'out' matrix into a data-frame (to enable the use of '$' to access a column by name)
out <- as.data.frame(out)

# total population
pop<-out[,"S"]+out[,"E"]+out[,"I"]+out[,"R"]
# weekly incidence
# inc <- params["report"]*params["gamma"]*out[,"E"]
inc <- params["sigma"]*out[,"E"]*pop_size
# make a new panel
par(mfrow=c(1,2))
# more plots
time <- out[,"time"]
plot(time,pop,type='l',lwd=3)
plot(time,inc,type='l',lwd=3)


######################################### #
# PLOT RESULTS                       ----
######################################### #

# # open pdf-stream (optional)
# pdf('SIR_example.pdf',6,4) 

# plot susceptible class
par(mfrow = c(1,1))
plot(out$S,
     type = 'l',
     xlab = 'Time',
     ylab = 'Population fraction',
     ylim = c(0,1.3),
     col  = 4)
lines(out$R, col=3)
lines(out$I, col=2)
lines(out$E, col=1)

# add legend
legend(x      = 'top', 
       legend = c('Susceptible','Infectious','Recovered',"Exposed"), 
       col    = c(4,2,3,1), 
       lwd    = 1,
       ncol   = 3,
       cex    = 0.9)

# # close pdf-stream (optional)
# dev.off()

