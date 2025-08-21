library(tidyverse)
library(deSolve)

## From Exel loop to deSolve

#Ordinary Differential Equation **solver**
  
#- Eliminates the need to manually choose Delta

#- Automatic Step Size Control: Unlike before for the fixed Delta,
#  ODE solvers automatically adjust step sizes for optimal accuracy



# Define the SIS model function
sis_model <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    # Calculate derivatives
    dS <- -beta * S * I / N + gamma * I # equivalent to S[i+1] = S[i] + Delta*(-(beta*I[i]/N)*S[i] + gamma*I[i])
    dI <- beta * S * I / N - gamma * I  # equivalent to I[t+1] = I[t] + Delta*(beta*I[t]/N*S[t] - gamma*I[t])
    
    # Return derivatives as list
    return(list(c(dS, dI)))
  })
}

# Parameters
parameters <- c(
  beta = 0.3,   # transmission rate
  gamma = 0.1,  # recovery rate  
  N = 500      # total population
)

# Initial conditions
initial_state <- c(
  # N = S +I
  S = 499,  # susceptible
  I = 1     # infected
)

# Time points
times <- seq(0, 365)

# Solve the ODE system
result <- ode(
  y = initial_state,
  times = times,
  func = sis_model,
  parms = parameters
)

# Convert to data frame and view
result_df <- as_tibble(result)
head(result_df)


model_plot <- result_df %>% 
  group_by(time) %>% 
  mutate(total_pop = sum(S,I)) %>% 
  pivot_longer(S:I)


ggplot(model_plot)+
  geom_line(aes(time, value, col=name), linewidth=1.5)+
  scale_x_continuous("Time")+
  scale_y_continuous("Population")+
  scale_color_brewer("Compartments", palette = "Set2")+
  ggtitle("SIS model simulation")+
  theme_bw()

