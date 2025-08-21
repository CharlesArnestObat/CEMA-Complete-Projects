library(tidyverse)
library(deSolve)
SIS_model=function(t,InitialState,Parameters){
  with(as.list(c(InitialState,Parameters)),{
    #Define the  equations governing your model
    N=S+I
    DS=-beta*S*I/N+gamma*I
    DI=beta*S*I/N-gamma*I
    return(list(c(DS,DI)))#Returns the list of derivatives you what the ode solver to solve
  })
  
}
#Initialise the initial states as a vector
iniatiaConditions=c(S=499,I=1)
#Give values of the parameters as a vector
parametrs=c(beta=0.3,gamma=0.1)
#define the time step/the time the simulation should take using the seq() function eg from day 0 to day365 at an interval of 1
TimeStep=seq(0,365, by =1)
#use ode solver to do the working for you by passing in the required paramemeters
#out <- ode(y = init, times = times, func = sir_model, parms = params)
odesolution2BSimulated <- ode(y = iniatiaConditions, times = TimeStep, func =SIS_model, parms = parametrs)
#  Convert to dataframe
MyDataframe = as.data.frame(odesolution2BSimulated)



