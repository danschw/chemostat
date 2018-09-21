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
# current3 <- c(R = R_in, N = 100, I1=1, I2=1,S=1)
# current6 <- c(R = R_in, N = 100, I1=1, I2=1,I3=1, I4=1,I5=1,S=1)
# current9 <- c(R = R_in, N = 100, I1=1, I2=1,I3=1, I4=1,I5=1,I6=1, I7=1,I8=1,S=1)

#define the time to run model
t <- seq (0, 5000, by=1)

#solve the model and save the output
out <- lsoda(current1, t, chemostat1, params1)
plot (out)
out <- lsoda(current1v, t, chemostat1v, params1v)
plot (out)
# out <- lsoda(current3, t, chemostat3, params)
# plot (out)
# out <- lsoda(current6, t, chemostat6, params)
# plot (out)
# out <- lsoda(current9, t, chemostat9, params)
# plot (out)


out.l <- gather(as.data.frame(out), key='population', value='Conc', colnames(out)[-1])
# ggplot(out.l[out.l$population!='R',], aes(x=time, y=Conc))+
ggplot(out.l, aes(x=time, y=Conc))+
  geom_line(aes(color=population))+
  scale_y_log10()+
  theme_bw()+
  scale_colour_brewer(palette='Set1')



################
# is sporulation lost because spores have an equal mortality (washout)?
# to test this I will add a parameter that reduces sporulation washout: 
# pDs probability of spore dilution. 0<= pDs <=1.  
# pDs=0 spores are all retained. pDs=1 spores dilution equals all other cells.

chemostat1s <- function(t,current,params){
  with(as.list(c(params,current)), {
    #list the equations for the derivatives of each state variable
    #R is the resource
    dR <- D * (R_in - R) - mu * R * (N+Nmut) 
    
    #N are the sporulating bacteria 
    dN <- mu * R * N - D * N - prob_spor * N + prob_actv * S - loss * N
    # Nmut are the spo- mutants
    dNmut <- mu * R * Nmut - D*Nmut +loss * N
    
    #S are the spores 
    dS <- prob_spor * N - prob_actv * S - pDs*D*S
    
    
    results <- c(dR, dN,dNmut, dS)
    return (list(results))
  })
}

# spore retention
pDs <- 0.5
# all parameters together
params1s <- c(R_in, D, mu, prob_spor, prob_actv, loss, pDs )
#declare initial conditions for state variables
current1s <- c(R = R_in, N = 1e5, Nmut=0, S=0)

#define the time to run model
t.model <- 5000
t <- seq (0, t.model, by=1)

#solve the model and save the output
out <- lsoda(current1s, t, chemostat1s, params1s)
# plot (out)

out.l <- gather(as.data.frame(out), key='population', value='Conc', colnames(out)[-1])
# ggplot(out.l[out.l$population!='R',], aes(x=time, y=Conc))+
ggplot(out.l, aes(x=time, y=Conc))+
  geom_line(aes(color=population))+
  # scale_y_log10()+
  theme_bw()+
  scale_colour_brewer(palette='Set1')

# bifuricate over spore retention parameter
i.pDs <- sort(c(0.01,seq(0,1,0.1)))
#define the time to run model
t.model <- 5000
t <- seq (0, t.model, by=1)
#data frame to collect data
d.res <- data.frame(matrix(ncol = 5, nrow = length(i.pDs)))
colnames(d.res) <-  c('pDs','coex.t','N.t','S.t','Nmut.t')
d.res$pDs <- i.pDs
# list to save plots
l.plot <- list()
for (i in 1:length(i.pDs)){
  # spore retention
  pDs <- i.pDs[i]
  # all parameters together
  params1s <- c(R_in, D, mu, prob_spor, prob_actv, loss, pDs )
  #declare initial conditions for state variables
  current1s <- c(R = R_in, N = 1e5, Nmut=0, S=0)
  #solve the model and save the output
  out <- lsoda(current1s, t, chemostat1s, params1s)
  # save data
    #time of equal densities (difference is the lowest)
  d.res$coex.t[i] <-
    out[which.min(abs(out[,'N']+out[,'S']-out[,'Nmut'])),'time']
    #final densities
  d.res$N.t[i] <-
    out[nrow(out),'N']
  d.res$Nmut.t[i] <-
    out[nrow(out),'Nmut']
  d.res$S.t[i] <-
    out[nrow(out),'S']
  # plots
  out.l <- gather(as.data.frame(out), key='population', value='Conc', colnames(out)[-1])
  l.plot[[i]] <- 
    ggplot(out.l, aes(x=time, y=Conc))+
    geom_line(aes(color=population))+
    # scale_y_log10()+
    theme_bw()+
    scale_colour_brewer(palette='Set1')+
    ggtitle(paste0('pDs=',pDs))
}

# dev.off()
# plot(d.res$pDs,d.res$coex.t, type='l')

# barplot(log10(as.matrix(d.res[,3:5])), beside = T, ylim = range(pretty(-10:6)), col=rainbow(12), 
#         ylab = paste0("log of final density (t=",t.model,')'),
#         main="final densities in chemostat models\nwith varying spore washout (fraction of dilution rate)")
# legend("bottomright", legend = as.character(d.res$pDs), fill =rainbow(12), ncol=3)

pdf(file = "./figures/spore_washout_densities.pdf", width = 10,paper = 'a4r')
matplot(d.res$pDs,log10(as.matrix(d.res[,3:5])),type='l',col = c(2:4),lty=1, lwd=5,
        xlab = "spore washout (fraction of dilution rate)", ylab = paste0("log of final density (t=",t.model,')'),
        main="final densities in chemostat models with varying spore washout")
legend("right", legend = c('N','S', 'Nmut'),,col = c(2:4),lty=1, lwd=5)
dev.off()

pdf(file = "./figures/spore_washout_plots.pdf", width = 10,paper = 'a4r')
grid.arrange(grobs=l.plot,nrow=4, ncol=3)
dev.off()

###############
# Will the virus favor sporulation under any of these conditions?
chemostat1sv <- function(t,current,params){
  with(as.list(c(params,current)), {
    #list the equations for the derivatives of each state variable
    #R is the resource
    dR <- D * (R_in - R) - mu * R * (N+Nmut) 
    
    #N are the sporulating bacteria 
    dN <- mu * R * N - D * N - prob_spor * N + prob_actv * S - loss * N - ads * N * V 
    # Nmut are the spo- mutants
    dNmut <- mu * R * Nmut - D*Nmut +loss * N- ads * Nmut * V 
    
    #S are the spores 
    dS <- prob_spor * N - prob_actv * S -  pDs*D*S
    #V are the phages
    dV <- ads * N * V * (burst-1) + ads * Nmut * V * (burst-1) - D*V
    
    results <- c(dR, dN,dNmut, dS, dV)
    return (list(results))
  })
}


# spore retention
pDs <- 0.01
# all parameters together
params1sv <- c(R_in, D, mu, prob_spor, prob_actv, loss, pDs, ads, burst)

#declare initial conditions for state variables
current1sv <- c(R = R_in, N = 1e5, Nmut=0, S=0, V=1)

#define the time to run model
t.model <- 10000
t <- seq (0, t.model, by=1)

#solve the model and save the output
out <- lsoda(current1sv, t, chemostat1sv, params1sv)
plot (out)

out.l <- gather(as.data.frame(out), key='population', value='Conc', colnames(out)[-1])
# ggplot(out.l[out.l$population!='R',], aes(x=time, y=Conc))+
ggplot(out.l, aes(x=time, y=Conc))+
  geom_line(aes(color=population))+
  scale_y_log10()+
  theme_bw()+
  scale_colour_brewer(palette='Set1')

# bifuricate over spore retention parameter
i.pDs <- sort(c(0.01,seq(0,1,0.1)))
#define the time to run model
t.model <- 5000
t <- seq (0, t.model, by=1)
#data frame to collect data
d.res <- data.frame(matrix(ncol = 5, nrow = length(i.pDs)))
colnames(d.res) <-  c('pDs','coex.t','N.t','S.t','Nmut.t')
d.res$pDs <- i.pDs
# list to save plots
l.plot <- list()
for (i in 1:length(i.pDs)){
  # spore retention
  pDs <- i.pDs[i]
  # all parameters together
  params1sv <- c(R_in, D, mu, prob_spor, prob_actv, loss, pDs, ads, burst)
  #declare initial conditions for state variables
  current1sv <- c(R = R_in, N = 1e5, Nmut=0, S=0, V=1)
  #solve the model and save the output
  out <- lsoda(current1sv, t, chemostat1sv, params1sv)
  # save data
  #time of equal densities (difference is the lowest)
  d.res$coex.t[i] <-
    out[which.min(abs(out[,'N']+out[,'S']-out[,'Nmut'])),'time']
  #final densities
  d.res$N.t[i] <-
    out[nrow(out),'N']
  d.res$Nmut.t[i] <-
    out[nrow(out),'Nmut']
  d.res$S.t[i] <-
    out[nrow(out),'S']
  # plots
  out.l <- gather(as.data.frame(out), key='population', value='Conc', colnames(out)[-1])
  l.plot[[i]] <- 
    ggplot(out.l, aes(x=time, y=Conc))+
    geom_line(aes(color=population))+
    scale_y_log10()+
    theme_bw()+
    scale_colour_brewer(palette='Set1')+
    ggtitle(paste0('pDs=',pDs))
}

# dev.off()
# plot(d.res$pDs,d.res$coex.t, type='l')

# barplot(log10(as.matrix(d.res[,3:5])), beside = T, ylim = range(pretty(-10:6)), col=rainbow(12), 
#         ylab = paste0("log of final density (t=",t.model,')'),
#         main="final densities in chemostat models\nwith varying spore washout (fraction of dilution rate)")
# legend("bottomright", legend = as.character(d.res$pDs), fill =rainbow(12), ncol=3)

pdf(file = "./figures/spore_washout_densities_phage.pdf", width = 10,paper = 'a4r')
matplot(d.res$pDs,log10(as.matrix(d.res[,3:5])),type='l',col = c(2:4),lty=1, lwd=5,
        xlab = "spore washout (fraction of dilution rate)", ylab = paste0("log of final density (t=",t.model,')'),
        main="final densities in chemostat models with varying spore washout")
legend("right", legend = c('N','S', 'Nmut'),,col = c(2:4),lty=1, lwd=5)
dev.off()

pdf(file = "./figures/spore_washout_plots_phage.pdf", width = 10,paper = 'a4r')
grid.arrange(grobs=l.plot,nrow=4, ncol=3)
dev.off()