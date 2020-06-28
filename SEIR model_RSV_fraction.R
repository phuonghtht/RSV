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
R0                <- 3
num_days_infected <- 6.7
num_days_exposed  <- 4
num_days_waning   <- 230
infected_seeds    <- 100 # ? 

######################################### #
# INITIALIZE PARAMETERS AND POPULATION ----
######################################### #
# recovery parameter
nui  <- 1/num_days_infected
# transmission parameter
beta   <- R0*nui
#rate of movement from latent to infectious stage
gamma <- 1/ num_days_exposed
# duration of immunity
v <- 1/ num_days_waning 

# # population states
S <- 1 - (infected_seeds/pop_size)
E <- infected_seeds/pop_size
I <- 0
R <- 0

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
  E     <- states['E']
  I     <- states['I']
  R     <- states['R']
  beta  <- params['beta']
  gamma <- params['gamma']
  nui   <- params['nui']


  # calculate state changes
  dS <- -beta*S*I + v*R
  dE <- beta*S*I - gamma*E
  dI <- gamma*E - nui*I
  dR <- nui*I - v*R

  # return (dS, dI, dR) as a vector in a list (required for the 'solve' function)
  return(list(c(dS, dE, dI, dR)))
}


######################################### #
# SET FUNCTION PARAMETERS            ----
######################################### #
# set time frame
times      <- seq(0, num_days, by = 1)

# set initial health states
states     <- c(S = S,E = E, I = I, R = R)

# set parameters
params     <- c(beta = beta, gamma = gamma, nui = nui, v = v)



######################################### #
# SOLVE ODE                          ----
######################################### #
# use the 'ode' function of deSolve package with our SIR function, health states and parameters
out <- ode(func = sirv_func, y = states, times = times, parms = params)
# plot(out)

# convert the 'out' matrix into a data-frame (to enable the use of '$' to access a column by name)
out <- as.data.frame(out)


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

