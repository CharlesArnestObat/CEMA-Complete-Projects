library(tidyverse)
library(deSolve)
SIS_Model=function(t,i,param){
  with(as.list(c(i,param)),{
    N=S+I
    dS=-beta*I/N*S+alpha*I
    dI=beta*I/N*S-alpha*I
    output=list(c(dS,dI))
    return(output)
  })
}
#initial Conditions
initiaConditions=c(
  S=499,
  I=1
)
param=c(beta=0.4,alpha=0.1)
#Timestep
time_Step=seq(0,365,by=1)
#solving ODE
solution=ode(times = time_Step,y=initiaConditions,func = SIS_Model,parms = param)
#convert to dataFrame
dF=as.data.frame(solution)%>% tibble()
head(dF)
matplot(dF$time, dF[,2:3], type = "l", lty = 1, col = c("blue", "red"),
        xlab = "Time (days)", ylab = "Number of people")
legend("right", legend = c("Susceptible", "Infected"),
       col = c("blue", "red"), lty = 1)

# Convert to data frame and view
result_df <- as_tibble(solution)
head(solution)


model_plot <- result_df %>% 
  group_by(time) %>% 
  mutate(total_pop = sum(S,I)) %>% 
  pivot_longer(S:I)


ggplot(model_plot)+
  geom_line(aes(time, value, col=name), linewidth=0.2)+
  scale_x_continuous("Time")+
  scale_y_continuous("Population")+
  scale_color_brewer("Compartments", palette = "Set2")+
  ggtitle("SIS model simulation")+
  theme_bw()
