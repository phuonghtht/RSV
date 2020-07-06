########################## #
#SEIR model for RSV-------
########################## #

setwd("/home/phuong/phuonght/gitKraken_project/RSV")
library(deSolve)

######################################### #
# MODEL SETTINGS                     ----
######################################### #
pop_size          <- 10000000
num_days          <- 500
num_weeks         <- 52*8 # belgium data : 2012- 2018 * 2011 data is incomplete?
num_days_infected <- 10   #[ 8-11], 6.7 is original estimate
num_days_exposed  <- 4    # [2,6]
num_days_waning   <- 160  # [148, 164] best fit 
infected_seeds    <- 1 # ? 1= 1000

######################################### #
# INITIALIZE PARAMETERS AND POPULATION ----
######################################### #
#https://steemit.com/science/@fouad/x-history-of-math-symbols

#  population states
# S2 <- 0.9 # change this initial stage do change the pattern of peaks
S1 <- 0.1
E1 <- infected_seeds/pop_size
I1 <- 0
R1 <- 0

E2 <- infected_seeds/(2*pop_size)#
I2 <- 0
R2 <- 0
S2 <- 1 - S1- E1- E2

######################################### #
# SET FUNCTION PARAMETERS            ----
######################################### #
# set time frame
times      <- seq(26, num_weeks, by = 1)

# set initial health states
states     <- c(S1 = S1,E1 = E1, I1 = I1, R1 = R1,
                S2 = S2,E2 = E2, I2 = I2, R2 = R2)

# set parameters
params     <- c(sigma = 7/ num_days_exposed, #rate of movement from latent to infectious stage
                gamma = 7/num_days_infected, # recovery rate
                nu=  7/ num_days_waning,     # rate of loss of immunity
                eta1 = 1/(2*52),             # aging rate 
                eta2 = 1/(78*52),            # aging rate 
                eta = 1/(80*52),             # death rate in general
                mu = 1/(80*52),              # birth rate
                beta1 = 0.65,                # the degree of seasonality, range [0,1],higher value stronger seasonal drivers
                phi = 2.43,                  # phase shift?
                beta0 = 1.99,                #average transmission rate
                R0 = 3,                      # reproduction number
                delta = 0.65,                # scaled susceptibility
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
         P <- (S1+E1+I1+R1+S2+E2+I2+R2)
         beta <- beta0*(1+beta1*cos(2*pi*t/52+phi))
         
         
         # calculate state changes
         dS1 <- mu - beta*S1*(I1+alpha*I2) - eta1*S1 + nu*R1
         dE1 <- beta*S1*(I1+alpha*I2) - eta1*E1 - sigma*E1 
         dI1 <- sigma*E1 - eta1*I1 - gamma*I1
         dR1 <- gamma*I1 - eta1*R1 - nu*R1
         dS2 <- eta1*S1 - delta*beta*S2*(I1+alpha*I2) - eta2*S2 + nu*R2
         dE2 <- eta1*E1 + delta*beta*S2*(I1+alpha*I2) - eta2*E2 - sigma*E2
         dI2 <- eta1*I1 + sigma*E2 - eta2*I2 - gamma*I2
         dR2 <- eta1*R1 + gamma*I2 - eta2*R2 - nu*R2
         
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
# plot(out)
# summary(out)

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
pop <- out[,"S1"] + out[,"E1"] + out[,"I1"] + out[,"R1"]+out[,"S2"] + out[,"E2"] + out[,"I2"] + out[,"R2"]
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

# par(mfrow = c(1,1))
# matplot(out[,1], out[,10:13], type = "l", xlab = "time", ylab = "population fraction")
# legend("topright", col = 1:8, lty = 1:8, legend = c("S", "E", "I","R"))

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

