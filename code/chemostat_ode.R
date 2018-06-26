# In this code I will model a chemostat in which a sporulating bacteria os grown

#Clear environment
rm(list=ls())

#Load necessary library
library(deSolve)
library(ggplot2)

## Model with multi stage sporulation
#Define system of differential equations for resource and bacteria
chemostat <- function(t,current,params){
  with(as.list(c(params,current)), {
    #list the equations for the derivatives of each state variable
    #R is the resource
    dR <- D * (R_in - R) - mu * R * N
    
    #N are the bacteria 
    dN <- mu * R * N - D * N - prob_spor * N + prob_actv * Sp3
    
    # Sporulation takes 3 time steps. modeled as 3 population maturig through the process
    dSp1 <- prob_spor * N - Sp1
    
    dSp2 <- Sp1 - Sp2
    
    dSp3 <- Sp2 - prob_actv * Sp3
    

    results <- c(dR, dN, dSp1, dSp1, dSp3)
    return (list(results))
  })
}



###define parameters
# resource in inflow
R_in <- 1000
#dilution rate
D <- 2
# growth rate
mu <- 2
# probaility of sporulation
prob_spor <- 0.001
#probability of reactivation of spores
prob_actv <- 0.001

# all parameters together
params <- c(R_in, D, mu, prob_spor, prob_actv)

#declare initial conditions for state variables
current0 <- c(R = R_in, N = 100, Sp1=0, Sp2=0, Sp3=0)

#define the time to run model
t <- seq (0, 5000, by=0.1)

#solve the model and save the output
out <- lsoda(current0, t, chemostat, params)

plot (out)


