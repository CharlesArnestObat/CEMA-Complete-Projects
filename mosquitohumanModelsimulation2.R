library(tidyverse)
library(deSolve)
vector_Mosquito_human<-function(t,x,parms){
  with(as.list(c(parms,x)),{
    # Total populations
    T_mosq_M=Sm+Im
    TotHuman_H=S+I
    dSm=mu_m*T_mosq_M-beta*I/TotHuman_H*Sm-mu_m*Sm
    dIm=beta*I/TotHuman_H*Sm-mu_m*Im
    #Human Equations
    dS=-alpha*Im/TotHuman_H*Sh+gamma*Ih
    dI=alpha*Im/TotHuman_H*Sh-gamma*Ih
    output=c(dSm,dIm,dS,dI)
    list(output)
  })
}
biting=0.3
seeking=0.01
infectivity=0.03
alpha=biting*seeking*infectivity
k=5 #mosquito to human raio

#define our parameters
parms=c(
  mu_m=0.1,
  alpha=alpha,
  beta=0.5,
  gamma=1/25
)
parms[1]
#ininitial conditions in the compartments
iniatiaConditions=c(
  Sm=4900,
  Im=100,
  Sh=990,
  Ih=10
  
)
iniatiaConditions[1]
#defining the simulation period
times=seq(0,300)
model_solution=ode(
  times = times,
  parms = parms,
  y=iniatiaConditions,
  func = vector_Mosquito_human
)
model_solution
#make the data plotable
df=as.data.frame(model_solution)
df
df=as_tibble(as.data.frame(model_solution))
df
df=as_tibble(as.data.frame(model_solution))%>% pivot_longer(names_to = "variable", cols = !1)
df
