# In this code I will model a chemostat in which a sporulating bacteria is grown
# Adding a phage into the chemostat
#Clear environment
rm(list=ls())
pd=par()
#Load necessary library
library(deSolve)
library(tidyverse)
library(gridExtra)

## Models with multi stage sporulation
#Define system of differential equations for resource and bacteria
chemostat0 <- function(t,current,params){
  with(as.list(c(params,current)), {
    #list the equations for the derivatives of each state variable
    #R is the resource
    dR <- D * (R_in - R) - mu * R * N 
    
    #N are the bacteria 
    dN <- mu*R * N - ads * N * V  -  D*N

    #V are the phages
    dV <- ads * N * V * (burst-1) - D*V
    results <- c(dR, dN, dV)
    return (list(results))
  })
}
chemostat1 <- function(t,current,params){
  with(as.list(c(params,current)), {
    #list the equations for the derivatives of each state variable
    #R is the resource
    dR <- D * (R_in - R) - mu * R * N 
    
    #N are the bacteria 
    dN <- mu*R * N + prob_actv * S - D * N - prob_spor * N - ads * N * V
    #S are the spores 
    dS <- prob_spor * N - prob_actv * S - D*S
    
    #V are the phages
    dV <- ads * N * V * (burst-1) - D*V
    results <- c(dR, dN, dS, dV)
    return (list(results))
  })
}
chemostat3 <- function(t,current,params){
  with(as.list(c(params,current)), {
    #list the equations for the derivatives of each state variable
    #R is the resource
    dR <- D * (R_in - R) - mu * R * N 
    
    #N are the bacteria 
    dN <- mu * R * N - D * N - prob_spor * N + prob_actv * S - ads * N * V
    
    # Sporulation takes 3 time steps. modeled as 3 population maturing through the process
    dI1 <- prob_spor * N - I1 - D*I1
    dI2 <- I1  - I2 - D*I2
    dS <- I2 - prob_actv * S - D*S
    #V are the phages
    dV <- ads * N * V * (burst-1) - D*V
    results <- c(dR, dN, dI1, dI2, dS,dV)
    return (list(results))
  })
}
chemostat6 <- function(t,current,params){
  with(as.list(c(params,current)), {
    #list the equations for the derivatives of each state variable
    #R is the resource
    dR <- D * (R_in - R) - mu * R * N 
    
    #N are the bacteria 
    dN <- mu * R * N - D * N - prob_spor * N + prob_actv * S  - ads * N * V
    
    # Sporulation takes 6 time steps. modeled as 6 population maturig through the process
    dI1 <- prob_spor * N - I1 - D*I1
    dI2 <- I1  - I2 - D*I2
    dI3 <- I2 - I3 - D*I3
    dI4 <- I3 - I4 - D*I4
    dI5 <- I4 - I5 - D*I5
    dS <- I5 - prob_actv * S - D*S
    #V are the phages
    dV <- ads * N * V * (burst-1) - D*V
    
    results <- c(dR, dN, dI1, dI2,dI3,dI4,dI5, dS, dV)
    return (list(results))
  })
}
chemostat9 <- function(t,current,params){
  with(as.list(c(params,current)), {
    #list the equations for the derivatives of each state variable
    #R is the resource
    dR <- D * (R_in - R) - mu * R * N 
    
    #N are the bacteria 
    dN <- mu * R * N - D * N - prob_spor * N + prob_actv * S  - ads * N * V
    
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
    #V are the phages
    dV <- ads * N * V * (burst-1) - D*V
    
    results <- c(dR, dN, dI1, dI2,dI3,dI4,dI5,dI6, dI7, dI8, dS, dV)
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

#adsorption constant rate (per hour per ml)
ads <- 1e-8

# busrt size
burst <- 100

# all parameters together
params <- c(R_in, D, mu, prob_spor, prob_actv, ads, burst)
params0 <- c(R_in, D, mu,  ads, burst)
#declare initial conditions for state variables
current0 <- c(R = R_in, N = 100,  V=1)
current1 <- c(R = R_in, N = 100, S=1, V=1)
current3 <- c(R = R_in, N = 100, I1=1, I2=1,S=1,V=1)
current6 <- c(R = R_in, N = 100, I1=1, I2=1,I3=1, I4=1,I5=1,S=1,V=1)
current9 <- c(R = R_in, N = 100, I1=1, I2=1,I3=1, I4=1,I5=1,I6=1, I7=1,I8=1,S=1, V=1)

#define the time to run model
t <- seq (0, 100, by=1)

l.plot <- list()
for(i in as.character(c(0,1,3,6,9))){
  if (i=='0'){
    out <- lsoda(current0, t, chemostat0, params0)
  } else if (i=='1') {
    out <- lsoda(current1, t, chemostat1, params)
  } else if (i=='3') {
    out <- lsoda(current3, t, chemostat3, params)
  } else if (i=='6') {
    out <- lsoda(current6, t, chemostat6, params)
  } else if (i=='9') {
    out <- lsoda(current9, t, chemostat9, params)
  }
  out.l <- gather(as.data.frame(out), key='population', value='Conc',colnames(out)[-1])
  out.l$population <- factor(out.l$population, levels = c('R','N','V','S'))
  l.plot[[i]] <- 
    ggplot(filter(out.l, population =='R'| population == 'N' | population == 'S'| population =='V'), aes(x=time, y=Conc))+
    geom_line(aes(color=population), size=2)+
    scale_y_log10()+
    theme_bw()+
    scale_colour_brewer(palette='Set1')+
    theme(legend.position = "bottom",
          legend.title = element_blank())+
    guides(col = guide_legend(nrow = 2))+
    ggtitle(paste(i, 'stage(s)'))
}

pdf(file = "./figures/Phage_sporeStages.pdf", width = 10,paper = 'a4r')
grid.arrange(grobs=l.plot,nrow=1)
dev.off()

