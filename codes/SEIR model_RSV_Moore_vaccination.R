########################## #
#SEIRV V1 -> R1 -------
#Only vaccinate newborns
########################## #

rm(list=ls())
setwd("/home/phuong/phuonght/gitKraken_project/RSV")
library(deSolve)
library(ggplot2)
library(dplyr)

######################################### #
# MODEL SETTINGS                     ----
######################################### #
pop_size          <- 2200000*3/4 # Western Australia pop.# metropolitan
# pop_size          <- 1861923; #2014 perth pop, 0-79 years,
num_days          <- 500
num_weeks         <- 52*190
num_days_infected <- 10  #[ 8-11], 6.7 is original estimate
num_days_exposed  <- 4   # [2,6]
num_days_waning   <- 160  # [148, 164] best fit
infected_seeds    <- 10 # ? 1= 1000
p                 <- 0# a praction p of newborn will be vaccinated

######################################### #
# INITIALIZE PARAMETERS AND POPULATION ----
######################################### #
#https://steemit.com/science/@fouad/x-history-of-math-symbols

#  population states

# # sum both groups = 1
# S1 <- 2/80                         # less than 2 yrs old
# E1 <- infected_seeds/2/(pop_size)
# I1 <- 0
# R1 <- 0
# V1 <- 1-p
# E2 <- (infected_seeds/2)/(pop_size) # greater than 2 yrs old
# I2 <- 0
# R2 <- 0
# S2 <- 1-S1-E1-I1-R1-E2-I2-R2- V1


# each age group sum = 1
E1 <- infected_seeds/2/(pop_size) # less than 2 yrs old
I1 <- 0
R1 <- 0
S1 <- 1-E1-I1-R1

E2 <-infected_seeds/2/(pop_size) # greater than 2 yrs old
I2 <- 0
R2 <- 0
S2 <- 1-E2-I2-R2

######################################### #
# SET FUNCTION PARAMETERS            ----
######################################### #
# set time frame
times      <- seq(0, num_weeks, by = 1)

# set initial health states
states     <- c(S1 = S1,E1 = E1, I1 = I1, R1 = R1,
                S2 = S2,E2 = E2, I2 = I2, R2 = R2)

# set parameters
params     <- c(#sigma = num_days_exposed,   #rate of movement from latent to infectious stage
                sigma = 1/0.57,              #best fit
                # gamma = num_days_waning,   # recovery rate
                gamma = 1/1.4,               # best fit
                # nu=  7/ num_days_waning,   # rate of loss of immunity
                nu=  0.0438,                 # rate of loss of immunity
                eta1 = 1/(2*52),             # aging rate
                eta2 = 1/(78*52),            # aging rate
                # eta = 1/(80*52),           # death rate in general
                mu = 1/(78*52),              # birth rate/deathrate
                beta1 = 0.6495,              # the degree of seasonality, range [0,1],higher value stronger seasonal drivers
                phi = 2.43,                  # phase shift?
                beta0 = 1.9896,              #average transmission rate
                vac = 1/(150/7),            # vaccine protection = 150 days
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
         beta <- beta0*(1+beta1*cos(2*pi*t/52+phi))
         # calculate state changes
         dS1 <- mu*(1-p) - beta*S1*(I1+alpha*I2) - eta1*S1             + nu*R1
         dE1 <-            beta*S1*(I1+alpha*I2) - sigma*E1 - eta1*E1
         dI1 <-                                    sigma*E1 - gamma*I1 - eta1*I1
         dR1 <-                                               gamma*I1 - nu*R1 - eta1*R1 + p*mu
    
         dS2 <-         eta1*S1 - delta*beta*S2*(I1+alpha*I2) - eta2*S2             + nu*R2 
         dE2 <-         eta1*E1 + delta*beta*S2*(I1+alpha*I2) - sigma*E2 - eta2*E2
         dI2 <-         eta1*I1                               + sigma*E2 - gamma*I2 - eta2*I2
         dR2 <-         eta1*R1                                          + gamma*I2 - nu*R2  - eta2*R2

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
times_output <- seq(num_weeks-(52*6)-36,num_weeks,1)# burning time ## shifting peaks
times_output <- times_output[1:312] # take exactly 6 years for later scaling of data
out <- out[out[,1] %in% times_output,] # skip the initial time points
out[,1] <- out[,1] - min(out[,1])# rescale time points

par(mfrow = c(1,1))
matplot(out[,1], out[,c(4,8)], type = "l", xlab = "time", ylab = "population fraction",
        main = "")
legend("topright", col = 1:2, lty = 1:2, legend = c("I1","I2"))

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
# pop <- out[,"S1"] + out[,"E1"] + out[,"I1"] + out[,"R1"]
# pop <- out[,"S"] + out[,"E"] + out[,"I"] + out[,"R"]

# weekly incidence
inc1 <- out[,"I1"]
inc2 <- out[,"I2"]
inc.less2yrs <- inc1*pop_size*(2/80)
# inc.less2yrs <- inc1*pop_size*(2/80)/(sum(inc.less2yrs)/3394)# scaling the incidence for fiiting the data
sum(inc.less2yrs )
inc.older2yrs<- inc2*pop_size*(78/80)
time <- out[,"time"]
# plot(time,pop,type='l',lwd=3,ylim = c(0,2))# checking the consistence of pop

# Plot incidence (<2yrs)

# jpeg(paste("~/Dropbox/2.Master of Epidemiology_Antwerp/year 2/thesis/data analysis/vac p=",p,".jpg"))
par(mfrow=c(1,2))
plot(time,inc.less2yrs,type='l',lwd=3, col = 1,xlab = "Time (weeks)",
     ylab = "Incidence",
     main = paste("AG1 vaccination p=", p))
plot(time,inc.older2yrs,type='l',lwd=2, col = 2,xlab = "Time (weeks)",
     main = paste("AG2 vaccination p=", p),
     ylab = "Incidence",
     ylim = c(0,max(inc.older2yrs)))
# dev.off()

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

# #Western Australia DATA-----------------
# # this data were extracted from Moore's paper.
# report.data <- read.csv("./RSV data/Australia_RSV_2yrs_Hannah Moore.csv")
#
# # plot the reported  data
# par(mfrow = c(1,1))
# plot(report.data$week,
#      report.data$cases,
#      ylim = c(0, max(report.data$cases)*1.5))
#
# # plot initial model
# lines(out$time,inc.less2yrs,col=2,lwd=2)


# ########################## #
# #SEIRV V1 -> S2   -------
# ########################## #
# 
# rm(list=ls())
# setwd("/home/phuong/phuonght/gitKraken_project/RSV")
# library(deSolve)
# library(ggplot2)
# library(dplyr)
# 
# ######################################### #
# # MODEL SETTINGS                     ----
# ######################################### #
# pop_size          <- 2200000*3/4 # Western Australia pop.# metropolitan
# # pop_size          <- 1861923; #2014 perth pop, 0-79 years,
# num_days          <- 500
# num_weeks         <- 52*190
# num_days_infected <- 10  #[ 8-11], 6.7 is original estimate
# num_days_exposed  <- 4   # [2,6]
# num_days_waning   <- 160  # [148, 164] best fit
# infected_seeds    <- 10 # ? 1= 1000
# p                 <- 0# a praction p of newborn will be vaccinated
# 
# ######################################### #
# # INITIALIZE PARAMETERS AND POPULATION ----
# ######################################### #
# #https://steemit.com/science/@fouad/x-history-of-math-symbols
# 
# #  population states
# 
# # # sum both groups = 1
# # S1 <- 2/80                         # less than 2 yrs old
# # E1 <- infected_seeds/2/(pop_size)
# # I1 <- 0
# # R1 <- 0
# # V1 <- 1-p
# # E2 <- (infected_seeds/2)/(pop_size) # greater than 2 yrs old
# # I2 <- 0
# # R2 <- 0
# # S2 <- 1-S1-E1-I1-R1-E2-I2-R2- V1
# 
# 
# # each age group sum = 1
# E1 <- infected_seeds/2/(pop_size) # less than 2 yrs old
# I1 <- 0
# R1 <- 0
# V1 <- 0
# S1 <- 1-E1-I1-R1-V1
# 
# E2 <-infected_seeds/2/(pop_size) # greater than 2 yrs old
# I2 <- 0
# R2 <- 0
# S2 <- 1-E2-I2-R2
# 
# ######################################### #
# # SET FUNCTION PARAMETERS            ----
# ######################################### #
# # set time frame
# times      <- seq(0, num_weeks, by = 1)
# 
# # set initial health states
# states     <- c(S1 = S1,E1 = E1, I1 = I1, R1 = R1, V1=V1,
#                 S2 = S2,E2 = E2, I2 = I2, R2 = R2)
# 
# # set parameters
# params     <- c(#sigma = num_days_exposed,   #rate of movement from latent to infectious stage
#   sigma = 1/0.57,              #best fit
#   # gamma = num_days_waning,   # recovery rate
#   gamma = 1/1.4,               # best fit
#   # nu=  7/ num_days_waning,   # rate of loss of immunity
#   nu=  0.0438,                 # rate of loss of immunity
#   eta1 = 1/(2*52),             # aging rate
#   eta2 = 1/(78*52),            # aging rate
#   # eta = 1/(80*52),           # death rate in general
#   mu = 1/(78*52),              # birth rate/deathrate
#   beta1 = 0.6495,              # the degree of seasonality, range [0,1],higher value stronger seasonal drivers
#   phi = 2.43,                  # phase shift?
#   beta0 = 1.9896,              #average transmission rate
#   vac = 1/(150/7),             # vaccine protection = 150 days
#   delta = 0.65,                # scaled susceptibility
#   alpha = 0.65                 # scaled infectiousness
# )
# 
# ######################################### #
# # CREATE SIRV FUNCTION               ----
# ######################################### #
# # define a function for the ODE-solver with 3 function parameters (required)
# #  * times  = current time point (not yet used here, in the absence of seasonality)
# #  * states = current health states, to calculate updates
# #  * params = epidemiogical parameters, used in the differential equations
# #
# # Note: updating the health states over time and keeping track of the changes, is handled by the
# # ode-function of the 'deSolve' package.
# sirv_func <- function(t, states, params) {
# 
#   with(as.list(c(states, params)),
#        {
# 
#          # define variables
#          beta <- beta0*(1+beta1*cos(2*pi*t/52+phi))
# 
#          # calculate state changes
#          dS1 <- mu*(1-p) - beta*S1*(I1+alpha*I2) - eta1*S1             + nu*R1
#          dE1 <-            beta*S1*(I1+alpha*I2) - sigma*E1 - eta1*E1
#          dI1 <-                                    sigma*E1 - gamma*I1 - eta1*I1
#          dR1 <-                                               gamma*I1 - nu*R1 - eta1*R1
#          dV1 <- p*mu - nu*V1
#          dS2 <- nu*V1 + eta1*S1 - delta*beta*S2*(I1+alpha*I2) - eta2*S2             + nu*R2 + nu*V1
#          dE2 <-         eta1*E1 + delta*beta*S2*(I1+alpha*I2) - sigma*E2 - eta2*E2
#          dI2 <-         eta1*I1                               + sigma*E2 - gamma*I2 - eta2*I2
#          dR2 <-         eta1*R1                                          + gamma*I2 - nu*R2  - eta2*R2
# 
#          # return (dS, dI, dR) as a vector in a list (required for the 'solve' function)
#          return(list(c(dS1, dE1, dI1, dR1,dV1,dS2, dE2, dI2, dR2)))
#        }
#   )
# }
# 
# ######################################### #
# # SOLVE ODE                          ----
# ######################################### #
# # use the 'ode' function of deSolve package with our SIR function, health states and parameters
# out <- ode(func = sirv_func, y = states, times = times, parms = params)
# plot(out)
# # summary(out)
# times_output <- seq(num_weeks-(52*6)-36,num_weeks,1)# burning time ## shifting peaks
# times_output <- times_output[1:312] # take exactly 6 years for later scaling of data
# out <- out[out[,1] %in% times_output,] # skip the initial time points
# out[,1] <- out[,1] - min(out[,1])# rescale time points
# 
# par(mfrow = c(1,1))
# matplot(out[,1], out[,c(4,8)], type = "l", xlab = "time", ylab = "population fraction",
#         main = "")
# legend("topright", col = 1:2, lty = 1:2, legend = c("I1","I2"))
# 
# ######################################### #
# # PLOT RESULTS                       ----
# ######################################### #
# 
# # convert the 'out' matrix into a data-frame (to enable the use of '$' to access a column by name)
# out <- as.data.frame(out)
# out$S <- out[,"S1"] + out[,"S2"]
# out$E <- out[,"E1"] + out[,"E2"]
# out$I <- out[,"I1"] + out[,"I2"]
# out$R <- out[,"R1"] + out[,"R2"]
# 
# 
# # weekly incidence
# inc1 <- out[,"I1"]
# inc2 <- out[,"I2"]
# inc.less2yrs <- inc1*pop_size*(2/80)
# # inc.less2yrs <- inc1*pop_size*(2/80)/(sum(inc.less2yrs)/3394)# scaling the incidence for fiiting the data
# sum(inc.less2yrs )
# inc.older2yrs<- inc2*pop_size*(78/80)
# time <- out[,"time"]
# # plot(time,pop,type='l',lwd=3,ylim = c(0,2))# checking the consistence of pop
# 
# # Plot incidence (<2yrs)
# 
# par(mfrow=c(1,2))
# # jpeg(paste("~/Dropbox/2.Master of Epidemiology_Antwerp/year 2/thesis/data analysis/vac p=",p,".jpg"))
# plot(time,inc.less2yrs,type='l',lwd=3, col = 1,xlab = "Time (weeks)",
#      ylab = "Incidence",
#      main = paste("AG1 vaccination p=", p))
# plot(time,inc.older2yrs,type='l',lwd=2, col = 2,xlab = "Time (weeks)",
#      main = paste("AG2 vaccination p=", p),
#      ylab = "Incidence",
#      ylim = c(0,max(inc.older2yrs)))
# # dev.off()
# 
# 
# 
# # #Western Australia DATA-----------------
# # # this data were extracted from Moore's paper.
# # report.data <- read.csv("./RSV data/Australia_RSV_2yrs_Hannah Moore.csv")
# #
# # # plot the reported  data
# # par(mfrow = c(1,1))
# # plot(report.data$week,
# #      report.data$cases,
# #      ylim = c(0, max(report.data$cases)*1.5))
# #
# # # plot initial model
# # lines(out$time,inc.less2yrs,col=2,lwd=2)

