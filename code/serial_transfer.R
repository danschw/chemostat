# implementing the srial transfer model of:
#Chaudhry, W. N., Pleška, M., Shah, N. N., Weiss, H., McCall, I. C., Meyer, J. R., … Levin, B. R. 
#(2018). Leaky resistance and the conditions for the existence of lytic bacteriophage. 
#PLOS Biology, 16(8), e2005971. https://doi.org/10.1371/journal.pbio.2005971

#Clear environment
rm(list=ls())
pd=par()
#Load necessary library
library(deSolve)
library(tidyverse)
library(gridExtra)

## Models with multi stage sporulation
#Define system of differential equations for resource and bacteria
serial <- function(t,current,params){
  with(as.list(c(params,current)), {
    #list the equations for the derivatives of each state variable
    #R is the resource (μg per mL)
    dR =  -cmax*N*R/(R+h) 
    
    #N are the bacteria 
    dN = e*cmax*N*R/(R+h)-m*N


    results <- c(dR, dN )
    return (list(results))
  })
}


event <- function(t, current, parms) {
  current['R'] <- 1e8
  current['N'] <- 0.1 * current['N']
  return(current)
  }

###define parameters
#Maximum consumption rate of resources by bacteria
cmax = 0.8
#Half-saturation constant of resource uptake by bacteria
h = 50
#Efficiency of converting resources to bacteria
e = 1
#mortality
m=0.1


# all parameters together
params <- c (cmax=cmax, h=h, e=e, m=m)

#declare initial conditions for state variables
current0 <- c(R = 1e8, N = 1e6)

#define the time to run model
t <- seq (0, 2400, by=1)
t.event <- t[t%%120==0][-1]
# run the model
 out <- lsoda(current0, t, serial, params)

out <- lsoda(current0, t, serial, params,
             events = list(func = event, time=t.event))
plot(out)
# plot(out, log='y')
