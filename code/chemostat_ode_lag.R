# In this code I will model a chemostat in which a sporulating bacteria os grown

#Clear environment
rm(list=ls())
pd=par()
#Load necessary library
library(deSolve)
library(tidyverse)

## Model with multi stage sporulation
#Define system of differential equations for resource and bacteria
chemostat <- function(t,current,params){
  with(as.list(c(params,current)), {
    if (t <= t.spo)
      Ilag <- 0
    else {
      Ilag <- lagvalue((t-t.spo),2) 
      Ilag <- prob_spor * Ilag
     
      # print(paste(t,Ilag))
    }
    
    
    
    #list the equations for the derivatives of each state variable
    #R is the resource
    dR <- D * (R_in - R) - mu * R * N 
    
    #N are the cells in vegetative growth
    dN <- mu * R * N + prob_actv * S - prob_spor * N  - D * N 
    
    # I are cells intiated to sporulate
    dI <- prob_spor * N - Ilag - D*I
    
    # S are the spores
    dS <- Ilag - prob_actv * S - D*S
    
   
    results <- c(dR, dN, dI, dS)
    return (list(results))
  })
}



###define parameters
# resource in inflow
R_in <- 1e6
#dilution rate
D <- 1
# growth rate
mu <- D
# time of sporulation parameter
t.spo <- 10

# probaility of sporulation
prob_spor <- 0.01
#probability of reactivation of spores
prob_actv <- 0.01

# all parameters together
params <- c(R_in, D, mu, t.spo, prob_spor, prob_actv)

#declare initial conditions for state variables
current0 <- c(R = R_in, N = 100, I=0,S=0)

#define the time to run model
t <- seq (0, 200, by=.1)

#solve the model and save the output
out <- dede(current0, t, chemostat, params)

plot (out,log='y')




out.l <- gather(as.data.frame(out), key='population', value='Conc', c("R","N","I","S"))
# ggplot(out.l[out.l$population!='R',], aes(x=time, y=Conc))+
ggplot(out.l, aes(x=time, y=Conc))+
  geom_line(aes(color=population))+
  scale_y_log10()+
  theme_bw()+
  scale_colour_brewer(palette='Set1')

rel <- as.data.frame(out)
rel$all <- apply(X = rel[,3:5],MARGIN = 1,FUN = sum)
rel$rel <- rel$S/rel$all
qplot(data=rel, x=time, y=rel)+geom_line()


#####
# spore ratio
###define parameters
# resource in inflow
R_in <- 1e6
#dilution rate
D <- 0.25
# growth rate
mu <- D
# time of sporulation parameter
t.spo <- 2

# probaility of sporulation
prob_spor <- 0.01
#probability of reactivation of spores
prob_actv <- 0.01

# all parameters together
params <- c(R_in, D, mu, t.spo, prob_spor, prob_actv)

#declare initial conditions for state variables
current0 <- c(R = R_in, N = 100, I=0,S=0)

#define the time to run model
t <- seq (0, 50, by=.1)


d.rates <- data.frame(dilution=rep(c(0.01,0.03,0.06,0.09,0.15,0.3,0.6,0.9,1),10), 
                      spo.time= as.numeric(rep(1:10,each=9)),ratio=NA)

for (i in 1:nrow(d.rates)){
  #dilution rate
  D <- d.rates$dilution[i]
  mu <- D
  t.spo <- d.rates$spo.time[i]
  # all parameters together
  params <- c(R_in, D, mu, t.spo, prob_spor, prob_actv)
  #solve the model and save the output
  out <- dede(current0, t, chemostat, params)
  rel <- as.data.frame(out)
  rel$all <- apply(X = rel[,3:5],MARGIN = 1,FUN = sum)
  rel$rel <- rel$I/rel$all
  d.rates$ratio[i] <- rel$rel[length(rel$rel)]
}

ggplot(d.rates, aes(x=dilution, y=ratio))+
  geom_line(aes(color=as.factor(spo.time)),size=1)+
  scale_y_log10()+
  ylab("spore frequency")+
  xlab("Dilution rate (hr^-1)")+ theme_bw()+
  annotation_logticks(sides = 'l')
  
