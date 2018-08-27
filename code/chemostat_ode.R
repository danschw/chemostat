# In this code I will model a chemostat in which a sporulating bacteria os grown

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
    #R is the resource
    dR <- D * (R_in - R) - mu * R * N 
    
    #N are the bacteria 
    dN <- mu * R * N - D * N - prob_spor * N + prob_actv * S
    
    # Sporulation takes 6 time steps. modeled as 6 population maturig through the process
    dSp1 <- prob_spor * N - Sp1 - D*Sp1
    
    dSp2 <- Sp1  - Sp2 - D*Sp2
    dSp3 <- Sp2 - Sp3 - D*Sp3
    dSp4 <- Sp3 - Sp4 - D*Sp4
    dSp5 <- Sp4 - Sp5 - D*Sp5

    dS <- Sp5 - prob_actv * S - D*S
    

    results <- c(dR, dN, dSp1, dSp2,dSp3,dSp4,dSp5, dS)
    return (list(results))
  })
}



###define parameters
# resource in inflow
R_in <- 1e6
#dilution rate
D <- 0.3
# growth rate
mu <- D
# probaility of sporulation
prob_spor <- 0.01
#probability of reactivation of spores
prob_actv <- 0.01

# all parameters together
params <- c(R_in, D, mu, prob_spor, prob_actv)

#declare initial conditions for state variables
current0 <- c(R = R_in, N = 100, Sp1=1, Sp2=1,Sp3=1, Sp4=1,Sp5=1,S=1)

#define the time to run model
t <- seq (0, 200, by=1)

#solve the model and save the output
out <- lsoda(current0, t, chemostat, params)

plot (out)



library(tidyr)
out.l <- gather(as.data.frame(out), key='population', value='Conc', c("R","N","Sp1","Sp2","Sp3","Sp4","Sp5","S"))
# ggplot(out.l[out.l$population!='R',], aes(x=time, y=Conc))+
ggplot(out.l, aes(x=time, y=Conc))+
  geom_line(aes(color=population))+
  scale_y_log10()+
  theme_bw()+
  scale_colour_brewer(palette='Set1')

rel <- as.data.frame(out)
rel$all <- apply(X = rel[,3:9],MARGIN = 1,FUN = sum)
rel$rel <- rel$S/rel$all
qplot(data=rel, x=time, y=rel)+geom_line()


#####
# spore ratio
###define parameters
# resource in inflow
R_in <- 1e6

# growth rate
mu <- D
# probaility of sporulation
prob_spor <- 0.01
#probability of reactivation of spores
prob_actv <- 0.01

#declare initial conditions for state variables
current0 <- c(R = R_in, N = 100, Sp1=1, Sp2=1,Sp3=1, Sp4=1,Sp5=1,S=1)

#define the time to run model
t <- seq (0, 500, by=1)


d.rates <- data.frame(rates=c(0.01,0.03,0.06,0.09,0.15,0.3,0.6,0.9,1), ratio=NA)

for (i in d.rates$rates){
  #dilution rate
  D <- i
  # all parameters together
  params <- c(R_in, D, mu, prob_spor, prob_actv)
  #solve the model and save the output
  out <- lsoda(current0, t, chemostat, params)
  rel <- as.data.frame(out)
  rel$all <- apply(X = rel[,3:9],MARGIN = 1,FUN = sum)
  rel$rel <- rel$S/rel$all
  d.rates$ratio[d.rates$rates==i] <- rel$rel[length(rel$rel)]
}

qplot(data=d.rates, x=rates, y=ratio)+geom_line()+
  # scale_y_log10(breaks=c(0.1,0.01, 0.001,0.0001))+
  ylab("spore frequency")+
  xlab("Dilution rate (hr^-1)")+ theme_bw()+
  annotation_logticks(sides = 'l')
  
