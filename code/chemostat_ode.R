# In this code I will model a chemostat in which a sporulating bacteria os grown

#Clear environment
rm(list=ls())

#Load necessary library
library(deSolve)
library(ggplot2)

## Starting off without sporulation
#Define system of differential equations for resource and bacteria
chemostat <- function(t,current,params){
  with(as.list(c(params,current)), {
    #list the equations for the derivatives of each state variable
    #R is the resource
    dR <- D * (R_in - R) - mu * R * N
    
    #N are the bacteria 
    dN <- mu * R * N - D * N

    results <- c(dR, dN)
    return (list(results))
  })
}





#declare initial conditions for state variables
current0 <- c(R = R_in, N = 10)

###declare parameters
# resource in inflow
R_in <- 50

# Dilution rate
params <- c(D = 0.5, mu = 3)

#define the time to run model
t <- seq (0, 1000, by=0.1)

#solve the model and save the output
out <- lsoda(current0, t, chemostat, params)

plot (out, log='y')


