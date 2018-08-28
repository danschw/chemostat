# In this code I will model a chemostat in which a spores re added at t=0

#Clear environment
rm(list=ls())
pd=par()
#Load necessary library
library(deSolve)
library(ggplot2)

## Model with multi stage sporulation
#Define system of differential equations for resource and bacteria
chemostat <- function(t,current,params){
  with(as.list(c(params,current)), {
    #list the equations for the derivatives of each state variable
 
    #N are the bacteria 
    dS <- -D * S 
    
    results <- c( dS)
    return (list(results))
  })
}

###define parameters

#dilution rate
D <- 0.1


# all parameters together
params <- c( D)

#declare initial conditions for state variables
current0 <- c( S = 1e6)

#define the time to run model
t <- seq (0, 60, by=1)
#solve the model and save the output
out1 <- lsoda(current0, t, chemostat, params)

plot (out1, col=1, log='y')

#dilution rate
D <- 0.25
#solve the model and save the output
out2 <- lsoda(current0, t, chemostat, params)

points (out2,type = 'l',col=2)

#dilution rate
D <- 0.5
#solve the model and save the output
out3 <- lsoda(current0, t, chemostat, params)

points (out3,type = 'l',col=3)
abline(h=5e5)
abline(h=1e5)
abline(h=1e4)
