# In this code I will model a chemostat in which a sporulating bacteria os grown

#Clear environment
rm(list=ls())
pd=par()
#Load necessary library
library(deSolve)
library(tidyverse)

## Models with multi stage sporulation
#Define system of differential equations for resource and bacteria
chemostat1 <- function(t,current,params){
  with(as.list(c(params,current)), {
    #list the equations for the derivatives of each state variable
    #R is the resource
    dR <- D * (R_in - R) - mu * R * (N+Nmut) 
    
    #N are the sporulating bacteria 
    dN <- mu * R * N - D * N - prob_spor * N + prob_actv * S - loss * N
    # Nmut are the spo- mutants
    dNmut <- mu * R * Nmut - D*Nmut +loss * N
    
    #S are the spores 
    dS <- prob_spor * N - prob_actv * S - D*S
    
    
    results <- c(dR, dN,dNmut, dS)
    return (list(results))
  })
}

chemostat1v <- function(t,current,params){
  with(as.list(c(params,current)), {
    #list the equations for the derivatives of each state variable
    #R is the resource
    dR <- D * (R_in - R) - mu * R * (N+Nmut) 
    
    #N are the sporulating bacteria 
    dN <- mu * R * N - D * N - prob_spor * N + prob_actv * S - loss * N - ads * N * V 
    # Nmut are the spo- mutants
    dNmut <- mu * R * Nmut - D*Nmut +loss * N- ads * Nmut * V 
    
    #S are the spores 
    dS <- prob_spor * N - prob_actv * S - D*S
    #V are the phages
    dV <- ads * N * V * (burst-1) + ads * Nmut * V * (burst-1) - D*V
    
    results <- c(dR, dN,dNmut, dS, dV)
    return (list(results))
  })
}
chemostat3 <- function(t,current,params){
  with(as.list(c(params,current)), {
    #list the equations for the derivatives of each state variable
    #R is the resource
    dR <- D * (R_in - R) - mu * R * N 
    
    #N are the bacteria 
    dN <- mu * R * N - D * N - prob_spor * N + prob_actv * S
    
    # Sporulation takes 3 time steps. modeled as 3 population maturing through the process
    dI1 <- prob_spor * N - I1 - D*I1
    dI2 <- I1  - I2 - D*I2
    dS <- I2 - prob_actv * S - D*S
    
    results <- c(dR, dN, dI1, dI2, dS)
    return (list(results))
  })
}
chemostat6 <- function(t,current,params){
  with(as.list(c(params,current)), {
    #list the equations for the derivatives of each state variable
    #R is the resource
    dR <- D * (R_in - R) - mu * R * N 
    
    #N are the bacteria 
    dN <- mu * R * N - D * N - prob_spor * N + prob_actv * S
    
    # Sporulation takes 6 time steps. modeled as 6 population maturig through the process
    dI1 <- prob_spor * N - I1 - D*I1
    dI2 <- I1  - I2 - D*I2
    dI3 <- I2 - I3 - D*I3
    dI4 <- I3 - I4 - D*I4
    dI5 <- I4 - I5 - D*I5
    dS <- I5 - prob_actv * S - D*S
    
    results <- c(dR, dN, dI1, dI2,dI3,dI4,dI5, dS)
    return (list(results))
  })
}
chemostat9 <- function(t,current,params){
  with(as.list(c(params,current)), {
    #list the equations for the derivatives of each state variable
    #R is the resource
    dR <- D * (R_in - R) - mu * R * N 
    
    #N are the bacteria 
    dN <- mu * R * N - D * N - prob_spor * N + prob_actv * S
    
    # Sporulation takes 9 time steps. modeled as 9 population maturing through the process
    dI1 <- prob_spor * N - I1 - D*I1
    dI2 <- I1 - I2 - D*I2
    dI3 <- I2 - I3 - D*I3
    dI4 <- I3 - I4 - D*I4
    dI5 <- I4 - I5 - D*I5
    dI6 <- I5 - I6 - D*I6
    dI7 <- I6 - I7 - D*I7
    dI8 <- I7 - I8 - D*I8
    dS  <- I8 - prob_actv * S - D*S
    
    results <- c(dR, dN, dI1, dI2,dI3,dI4,dI5,dI6, dI7, dI8, dS)
    return (list(results))
  })
}


###define parameters
# resource in inflow
R_in <- 1e6
#dilution rate
D <- 0.1
# growth rate
mu <- D
# probaility of sporulation
prob_spor <- 0.01
#probability of reactivation of spores
prob_actv <- 0.01
#loss of sporulation by mutation
loss <- 1e-6
#adsorption constant rate (per hour per ml)
ads <- 1e-8
# busrt size
burst <- 100

# all parameters together
params1 <- c(R_in, D, mu, prob_spor, prob_actv, loss)
params1v <- c(R_in, D, mu, prob_spor, prob_actv, loss, ads, burst)

#declare initial conditions for state variables
current1 <- c(R = R_in, N = 1e5, Nmut=0, S=0)
current1v <- c(R = R_in, N = 1e5, Nmut=0, S=0, V=1)
current3 <- c(R = R_in, N = 100, I1=1, I2=1,S=1)
current6 <- c(R = R_in, N = 100, I1=1, I2=1,I3=1, I4=1,I5=1,S=1)
current9 <- c(R = R_in, N = 100, I1=1, I2=1,I3=1, I4=1,I5=1,I6=1, I7=1,I8=1,S=1)

#define the time to run model
t <- seq (0, 5000, by=1)

#solve the model and save the output
out <- lsoda(current1, t, chemostat1, params1)
plot (out)
out <- lsoda(current1v, t, chemostat1v, params1v)
plot (out)
out <- lsoda(current3, t, chemostat3, params)
plot (out)
out <- lsoda(current6, t, chemostat6, params)
plot (out)
out <- lsoda(current9, t, chemostat9, params)
plot (out)


out.l <- gather(as.data.frame(out), key='population', value='Conc', colnames(out)[-1])
# ggplot(out.l[out.l$population!='R',], aes(x=time, y=Conc))+
ggplot(out.l, aes(x=time, y=Conc))+
  geom_line(aes(color=population))+
  scale_y_log10()+
  theme_bw()+
  scale_colour_brewer(palette='Set1')

rel <- as.data.frame(out)
rel$all <- apply(X = rel[,3:ncol(rel)],MARGIN = 1,FUN = sum)
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


d.rates <- data.frame(dilution=rep(c(0.01,0.03,0.06,0.09,0.15,0.3,0.6,0.9,1),4), 
                      spo.time= as.numeric(rep(c(1,3,6,9),each=9)),ratio=NA, N=NA, S=NA)

for (i in 1:nrow(d.rates)){
  #dilution rate
  D <- d.rates$dilution[i]
  mu <- D
  # all parameters together
  params <- c(R_in, D, mu, prob_spor, prob_actv)
  t.spo <- d.rates$spo.time[i]
  if (t.spo==1){
    current1 <- c(R = R_in, N = 100, S=1)
    out <- lsoda(current1, t, chemostat1, params)
  
    } else if (t.spo==3){
    current3 <- c(R = R_in, N = 100, I1=1, I2=1,S=1)
    out <- lsoda(current3, t, chemostat3, params)
  
    } else if (t.spo==6){
    current6 <- c(R = R_in, N = 100, I1=1, I2=1,I3=1, I4=1,I5=1,S=1)
    out <- lsoda(current6, t, chemostat6, params)
  
    } else if (t.spo==9){
    current9 <- c(R = R_in, N = 100, I1=1, I2=1,I3=1, I4=1,I5=1,I6=1, I7=1,I8=1,S=1)
    out <- lsoda(current9, t, chemostat9, params)
  
    }

  rel <- as.data.frame(out)
  rel$all <- apply(X = rel[,3:ncol(rel)],MARGIN = 1,FUN = sum)
  rel$rel <- rel$S/rel$all
  d.rates$ratio[i] <- rel$rel[nrow(rel)]
  d.rates$N[i] <- rel$N[nrow(rel)]
  d.rates$S[i] <- rel$S[nrow(rel)]
}
p1 <- 
  ggplot(d.rates, aes(x=dilution, y=ratio))+
    geom_line(aes(color=as.factor(spo.time)),size=1)+
    geom_point(aes(shape=as.factor(spo.time),color=as.factor(spo.time)),size=3)+
    scale_color_discrete(name="Sporulation time")+
    scale_shape_discrete(name="Sporulation time")+
    scale_y_log10(breaks=c(1e-1,1e-2,1e-3,1e-4))+
    ylab("spore frequency (log)")+
    xlab("Dilution rate (hr^-1)")+ theme_bw()+
    theme(legend.position = "bottom")+
    annotation_logticks(sides = 'l')+
    ggtitle("spore frequency")

p2 <- 
  ggplot(d.rates, aes(x=dilution, y=S))+
    geom_line(aes(color=as.factor(spo.time)),size=1)+
    geom_point(aes(shape=as.factor(spo.time),color=as.factor(spo.time)),size=3)+
    scale_color_discrete(name="Sporulation time")+
    scale_shape_discrete(name="Sporulation time")+
    scale_y_log10(limits=c(1e1,1e6))+
    ylab("Spore conc. (log)")+
    xlab("Dilution rate (hr^-1)")+ theme_bw()+
    theme(legend.position = "bottom")+
    annotation_logticks(sides = 'l')+
    ggtitle("Spore concentration")
p3 <- 
  ggplot(d.rates, aes(x=dilution, y=N))+
    geom_line(aes(color=as.factor(spo.time)),size=1)+
    geom_point(aes(shape=as.factor(spo.time),color=as.factor(spo.time)),size=3)+
    scale_color_discrete(name="Sporulation time")+
    scale_shape_discrete(name="Sporulation time")+
    scale_y_log10(limits=c(1e1,1e6))+ 
    ylab("Vegetative cell conc.")+
    xlab("Dilution rate (hr^-1)")+ theme_bw()+
    theme(legend.position = "bottom")+
    annotation_logticks(sides = 'l')+
    ggtitle("Vegetative cell concentration")

library(gridExtra)
grid.arrange(p1,p2,p3,nrow=1)

