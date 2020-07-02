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
#https://steemit.com/science/@fouad/x-history-of-math-symbols

#  population states
S <- 1 - (infected_seeds/pop_size)
E <- infected_seeds/pop_size
I <- 0
R <- 0

######################################### #
# SET FUNCTION PARAMETERS            ----
######################################### #
# set time frame
times      <- seq(1, num_weeks, by = 1)

# set initial health states
states     <- c(S = S,E = E, I = I, R = R)

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
         P <- (S+E+I+R)
         beta <- beta0*(1+beta1*cos(2*pi*t/52+phi))
         
         
         # calculate state changes
         dS <- mu *P- beta*S*I + nu*R - eta*S
         dE <- beta*S*I - sigma*E - eta*E
         dI <- sigma*E - gamma*I - eta*I
         dR <- gamma*I - nu*R - eta*R
         
         # return (dS, dI, dR) as a vector in a list (required for the 'solve' function)
         return(list(c(dS, dE, dI, dR)))
       }
  )
}

######################################### #
# SOLVE ODE                          ----
######################################### #
# use the 'ode' function of deSolve package with our SIR function, health states and parameters
out <- ode(func = sirv_func, y = states, times = times, parms = params)
plot(out)
summary(out)

# par(mfrow = c(1,1))
# matplot(out[,1], out[,2:5], type = "l", xlab = "time", ylab = "population fraction")
# legend("topright", col = 1:4, lty = 1:4, legend = c("S", "E", "I","R"))


######################################### #
# PLOT RESULTS                       ----
######################################### #

# Total population
pop <- out[,"S"] + out[,"E"] + out[,"I"] + out[,"R"]
# weekly incidence
# inc <- params["report"]*params["gamma"]*out[,"E"]
inc <- params["sigma"]*out[,"E"]
# make a new panel
par(mfrow=c(1,2))
# more plots
time <- out[,"time"]
plot(time,pop,type='l',lwd=3,ylim = c(0,2))
plot(time,inc,type='l',lwd=3)

# convert the 'out' matrix into a data-frame (to enable the use of '$' to access a column by name)
out <- as.data.frame(out)
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


