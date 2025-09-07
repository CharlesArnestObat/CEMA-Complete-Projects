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
    dS=-alpha*Im/TotHuman_H*S+gamma*I
    dI=alpha*Im/TotHuman_H*S-gamma*I
    output=c(dSm,dIm,dS,dI)
    list(output)
  })
}
#ininitial conditions in the compartments
iniatiaConditions=c(
  Sm=4000,
  Im=3000,
  S=5000,
  I=1000
)
#define our parameters
parms=c(
  mu_m=1/15,
  alpha=0.12,
  beta=0.25,
  gamma=1/20
)
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
df %>% ggplot()+geom_line(aes(x=time,y=value,colour = variable))+
  theme_minimal()+
  labs(title="Mosquito/Vector-Human Compartments",y=("population"),colour="species")+
  facet_wrap(~variable)