########################## #
#SEIR model for RSV-------
########################## #

setwd("/home/phuong/phuonght/gitKraken_project/RSV")

library(deSolve)

######################################### #
# MODEL SETTINGS                     ----
######################################### #
pop_size          <- 100000
num_days          <- 500
num_weeks         <- 52*7 # belgium data : 2012- 2018 * 2011 data is incomplete?
num_days_infected <- 10 #[ 8-11], 6.7 is original estimate
num_days_exposed  <- 4 # [2,6]
num_days_waning   <- 160 # [148, 164] best fit 
infected_seeds    <- 1 # ? 


######################################### #
# INITIALIZE PARAMETERS AND POPULATION ----
######################################### #
initP <- 100000 # population size
initE <- 1 # Exposed
initI <- 0 # Infectious
initR <- 0 # Immune
initS <- initP - initE - initI - initR # Susceptible
# set time frame
times      <- seq(1, num_weeks, by = 1)
# set initial health states
states <- c(S = initS, E = initE, I = initI, R = initR)
# set parameters
params     <- c(sigma = 7/ num_days_exposed, #rate of movement from latent to infectious stage
                gamma = 7/num_days_infected, # recovery rate
                nu=  7/ num_days_waning,     # rate of loss of immunity
                eta1 = 1/(2*52),             # aging rate 
                eta2 = 1/(78*52),            # aging rate 
                eta = 1/(80*52),            # death rate in general
                mu = 1/(80*52),              # birth rate
                beta1 = 0.65,                # the degree of seasonality, range [0,1],higher value stronger seasonal drivers
                phi = 2.43,                  # phase shift?
                beta0 = 1.99,                #average transmission rate
                R0 = 3,                      # reproduction number
                deta = 0.65,                 # scaled susceptibility
                alpha = 0.65                 # scaled infectiousness
)  

                

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
RSV<-function(t, states, params) 
{
  with(as.list(c(states, params)),
       {
         
         # define variables
         P <- (S+E+I+R)
         beta <- beta0*(1+beta1*cos(2*pi*t/52+phi))
         lam <- beta*I/P
         
         
         # rate of change
         dS <- mu*P -lam*S + nu*R -eta*S
         dE <- lam*S - sigma*E - eta*E
         dI <- sigma*E - gamma*I - eta*I
         dR <- gamma*I - nu*R - eta*R
         
         # return the rate of change
         list(c(dS, dE, dI, dR))
       }
  ) 
  
}


# run the model
out <- ode(y = states, times = times, func = RSV, parms = params)
# a simple plot of the model output
plot(out)

# some more model outputs
# total population
pop<-out[,"S"]+out[,"E"]+out[,"I"]+out[,"R"]
# weekly incidence
inc <- params["sigma"]*out[,"E"]
# make a new panel
par(mfrow=c(1,2))
# more plots
time<-out[,"time"]
plot(time,pop,type='l',lwd=3, ylim = c(0,pop_size*2))
plot(time,inc,type='l',lwd=3)

# convert the 'out' matrix into a data-frame (to enable the use of '$' to access a column by name)
out1 <- as.data.frame(out)
# plot susceptible class
par(mfrow = c(1,1))
plot(out1$S,
     type = 'l',
     xlab = paste('Time'),
     ylab = 'Population',
     ylim = c(0,120000),
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

