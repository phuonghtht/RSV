###################################### #
#I. SEIRV_age model for RSV -------
###################################### #
rm(list=ls())
#setwd("/home/phuong/phuonght/gitKraken_project/RSV")
library(deSolve)
library(ggplot2)


#Changeable variables  ##########################
vac.frac = 0

######################################### #
# MODEL SETTINGS                     ----
######################################### #
pop_size          <- 1100000 # Belgium pop size
num_days          <- 500
num_weeks         <- 52*20
num_days_infected <- 10  #[ 8-11], 6.7 is original estimate
num_days_exposed  <- 4   # [2,6]
num_days_waning   <- 200 # [148, 164] best fit
infected_seeds    <- 10 # ? 1= 1000


######################################### #
# INITIALIZE PARAMETERS AND POPULATION ----
######################################### #
#https://steemit.com/science/@fouad/x-history-of-math-symbols

#  population states

S1 <- 1/80 # less than 1 yrs old
E1 <- infected_seeds/2/(pop_size)
I1 <- 0
R1 <- 0

V  <- 0

E2 <- (infected_seeds/2)/(pop_size) # greater than 1 yrs old #LW
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
                mu = 1/(79*52),              # birth rate
                beta1 = 0.89,                # the degree of seasonality, range [0,1],higher value stronger seasonal drivers
                phi = 2.43,                  # phase shift?
                beta0 = 1.99,                # average transmission rate
                # R0 = 3,                    # reproduction number
                delta = 0.54,                # scaled susceptibility
                alpha = 0.76,                # scaled infectiousness
                p = vac.frac                 # vaccination proportion
                )

######################################### #
# CREATE SIRV FUNCTION               ----
######################################### #
# Note: updating the health states over time and keeping track of the changes, is handled by the
# ode-function of the 'deSolve' package.
seirv_func <- function(t, states, params) {

  with(as.list(c(states, params)),
       {
         # define variables
         # P <- (S1+E1+I1+R1+S2+E2+I2+R2)
         beta <- beta0*(1+beta1*cos(2*pi*t/52+phi))


         # calculate state changes
         dS1 <- (1-p)* mu - beta*S1*(I1+alpha*I2) - eta1*S1                         + nu*R1
         dE1 <-             beta*S1*(I1+alpha*I2) - eta1*E1 - sigma*E1
         dI1 <-                                   - eta1*I1 + sigma*E1  - gamma*I1
         dR1 <-                                   - eta1*R1             + gamma*I1  - nu*R1

         dV <- p*mu  - nu*V

         dS2 <- nu*V + eta1*S1 - delta*beta*S2*(I1+alpha*I2) - eta2*S2 + nu*R2
         dE2 <-        eta1*E1 + delta*beta*S2*(I1+alpha*I2) - eta2*E2 - sigma*E2
         dI2 <-        eta1*I1                               - eta2*I2 + sigma*E2 - gamma*I2
         dR2 <-        eta1*R1                               - eta2*R2            + gamma*I2 - nu*R2

         # return (dS, dI, dR) as a vector in a list (required for the 'solve' function)
         return(list(c(dS1, dE1, dI1, dR1,dV, dS2, dE2, dI2, dR2)))
       }
  )
}

######################################### #
# SOLVE ODE                          ----
######################################### #
# use the 'ode' function of deSolve package with our SIR function, health states and parameters
out <- ode(func = seirv_func, y = states, times = times, parms = params)
plot(out)

# select final 8 epidemics
times_output <- seq(num_weeks-(52*8)-2,num_weeks,1)# burning time ## shifting peaks
out <- out[out[,1] %in% times_output,] # skip the initial time points
out[,1] <- out[,1] - min(out[,1])# rescale time points

# convert the 'out' matrix into a data-frame (to enable the use of '$' to access a column by name)
out <- as.data.frame(out)
names(out)

# calculate age-specific incidence
out$I1_incidence <- out$I1*pop_size
out$I2_incidence <- out$I2*pop_size

######################################### #
# REFERENCE DATA                     ----
######################################### #

ref_data <- read.csv("./RSV data/RSV_cases_time_epistat.csv")
# total cases: 63301 , less than 2yrs: 56126, 1yr:13751
#  pro_case_less1yrs <- 13751/63301
ref_data$week <- seq(1,dim(ref_data)[1])
# ref_data$cases.less1yrs <- round(ref_data$cases*pro_case_less1yrs)
## score for initial model
# ggplot(ref_data, aes(week, cases.less1yrs))+
#   geom_point()

# LW
ref_data_age       <- read.csv("./RSV data/RSV_cases_age_epistat.csv")
pro_case_less1yrs  <- ref_data_age$cases[ref_data_age$age==0] / sum(ref_data_age$cases)    #LW
ref_data$cases_ag1 <- round(ref_data$cases*pro_case_less1yrs)
ref_data$cases_ag2 <- ref_data$cases - ref_data$cases_ag1

# convert date into "date format"
ref_data$date      <- as.Date(ref_data$date)


######################################### #
# PLOT RESULTS                       ----
######################################### #

# plot results
par(mfrow =c(1,1))
plot(out$time,out$I2_incidence,type='l',lwd=3, col = 4,main='model')
lines(out$time,out$I1_incidence,type='l',lwd=3, col = 3)
legend('topleft',
       c('AG1','AG2'),
       col = c(3,4),
       lwd=2,
       cex=0.7)




# plot reference and model data
par(mfrow=c(2,1))
y_lim <- range(ref_data$cases_ag1,ref_data$cases_ag2)*1.2
plot(ref_data$week,
     ref_data$cases_ag1,
     col = 1,lwd=1,pch=20,
     ylim = y_lim,
     type='p',ylab='incidence',main='Belgium - AG1')
lines(out$time,out$I1_incidence,col=3,lwd=2)

plot(ref_data$week,
     ref_data$cases_ag2,
     col = 1,lwd=1,pch=20,
     ylim = y_lim,
     type='p',ylab='incidence',main='Belgium - AG2')
lines(out$time,out$I2_incidence,col=4,lwd=2)



