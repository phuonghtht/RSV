########################## #
#SEIR model for RSV-------
########################## #

setwd("/home/phuong/Dropbox/2.Master of Epidemiology_Antwerp/INTERNSHIP/SIR model")

library(deSolve)

######################################### #
# MODEL SETTINGS                     ----
######################################### #
pop_size          <- 100000
num_days          <- 500
R0                <- 3 
num_days_infected <- 6.7
num_days_exposed  <- 4
num_days_waning   <- 230
infected_seeds    <- 100

##  latency period: 4 days [2-6]
##  infectious period : 6.7 days [ 1-21], 10 days [7-12]
##  duration of immunity: 230 days / 200 days

######################################### #
# INITIALIZE PARAMETERS AND POPULATION ----
######################################### #
initP <- 10000 # population size
initE <- 1 # Exposed
initI <- 0 # Infectious
initR <- 0 # Immune
initS <- initP - initE - initI - initR # Susceptible
# set time frame
times      <- seq(0, num_days, by = 1)
# set initial health states
states <- c(S = initS, E = initE, I = initI, R = initR)
# set parameters
params     <- c(gamma = 1/ num_days_exposed,#rate of movement from latent to infectious stage
                nui = 1/num_days_infected,  # recovery rate
                v =  1/ num_days_waning)    # rate of loss of immunity
                

######################################### #
# CREATE FUNCTION  to solve the equation             ----
######################################### #
# define a function for the ODE-solver with 3 function parameters (required)
#  * times  = current time point (not yet used here, in the absence of seasonality)
#  * states = current health states, to calculate updates
#  * params = epidemiogical parameters, used in the differential equations
#
# Note: updating the health states over time and keeping track of the changes, is handled by the
# ode-function of the 'deSolve' package.

# set up a function to solve the equations
RSV<-function(times, states, params) 
{
  with(as.list(c(states, params)),
       {
         
         # define variables
         P <- (S+E+I+R)
         beta<- R0*nui
         lam <- beta*I/P
         
         
         # rate of change
         dS <- -lam*S + v*R
         dE <- lam*S - gamma*E
         dI <- gamma*E - nui*I
         dR <- nui*I - v*R
         
         # return the rate of change
         list(c(dS, dE, dI, dR))
       }
  ) 
  
}


# run the model
out <- ode(y = states, times = times, func = RSV, parms = params)
# a simple plot of the model output
# plot(out)


# convert the 'out' matrix into a data-frame (to enable the use of '$' to access a column by name)
out1 <- as.data.frame(out)
# plot susceptible class
par(mfrow = c(1,1))
plot(out1$S,
     type = 'l',
     xlab = paste('Time'),
     ylab = 'Population',
     ylim = c(0,12000),
     col  = 4)
lines(out1$R, col=3)
lines(out1$I, col=2)
lines(out1$E, col=1)

# add legend
legend(x      = 'top', 
       legend = c('Susceptible','Infectious','Recovered',"Exposed"), 
       col    = c(4,2,3,1), 
       lwd    = 1,
       ncol   = 3,
       cex    = 0.9)

