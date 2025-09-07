
### Objective of the lesson

#Simulate an intervention in the models seen previously treatment,
#extraction parameters, and visualize the impact of the intervention on disease dynamics.
#Know how to analyze intervention impact on disease dynamics


## -----------------------------------------------------------------------------------------------------------------------------
#| label: loadLib

library(tidyverse) # for data manipulation and visualization
library(deSolve)  # for ODE solving
library(plotly)  # interactive plot


## -----------------------------------------------------------------------------------------------------------------------------
#| label: SIS_model
#| code-wrap: true
#| echo: true

#use your knowledge from the previous lesson to complete the code

sis_model <- function(time, compartment, parameters) {
  
  with(as.list(c(compartment, parameters)), {
    # Calculate derivatives
    
    # Return derivatives as list
    return(list(c()))
  })
}





## -----------------------------------------------------------------------------------------------------------------------------
#| label: init_sis

# Parameters
parameters <- c(
  #fill in the parameters
  
)


# Initial conditions
initial_state <- c(
  #fill in the initial conditions
)



## -----------------------------------------------------------------------------------------------------------------------------
#| label: simul_sis
#| code-wrap: true

# Time points
times <- seq(0, 365, by = 1)
# Solve the ODE system
result <- ode(
  y = initial_state,
  times = times,
  func = sis_model,
  parms = parameters
)

# Convert to data frame and view      
result_df <- as_tibble(result)




## -----------------------------------------------------------------------------------------------------------------------------
#| label: viz_sis
#| code-wrap: true
#| echo: true
#| fig-cap: "Dynamics of the SIS model over time"
#| fig-width: 10
#| fig-height: 8

toplot <- result_df %>% 
  group_by(time) %>% 
  mutate(total_pop = sum(S,I)) %>% 
  pivot_longer(S:I)
p1 = ggplot(toplot)+
  geom_line(aes(time, value, col=name), linewidth=1.5)+
  scale_x_continuous("Time")+
  scale_y_continuous("Population")+
  scale_color_brewer("Compartments", palette = "Set1")+
  ggtitle("SIS model simulation")+
  theme_bw()




## -----------------------------------------------------------------------------------------------------------------------------
#| label: obs_sis

plotly::ggplotly(p1)



## -----------------------------------------------------------------------------------------------------------------------------
#| label: R0_sis
#| code-wrap: true
#| echo: true
#| fig-cap: "Basic reproduction number for the SIS model"
#| fig-width: 10
#| fig-height: 8
#| fig-align: center
#| fig-pos: "center"
#| fig-alt: "Basic reproduction number for the SIS model"

R0_sis <- # write the formular 
  cat("The basic reproduction number (R0) for the SIS model is:", R0_sis, "\n")



## -----------------------------------------------------------------------------------------------------------------------------
#| label: SIST_expl
#| code-wrap: true

sits_model <- function(time, compartment, parameters) {
  
  with(as.list(c(compartment, parameters)), {
    # Calculate derivatives
    #fill in the derivatives
    dS <- 
    
    dI <- 
    
    dTrt <- 
      
      # Return derivatives as list
      return(list(c(dS, dI, dTrt )))
  })
}



## -----------------------------------------------------------------------------------------------------------------------------
#| label: init
# Parameters
parameters <- c(
  #fill in the parameters for the SITS model
  
)

# Initial conditions
initial_state <- c(
  #fill in the initial conditions for the SITS model
  
)



## -----------------------------------------------------------------------------------------------------------------------------
#| label: simul


# Time points
times <- seq(0, 365)

# Solve the ODE system
result <- ode(
  y = initial_state,
  times = times,
  func = sits_model,
  parms = parameters
)

# Convert to data frame and view
result_df <- as_tibble(result)
# head(result_df)



## -----------------------------------------------------------------------------------------------------------------------------

#| label: viz1
#| 
toplot <- result_df %>% 
  group_by(time) %>% 
  mutate(total_pop = sum(S,I)) %>% 
  pivot_longer(S:Trt)

# complete this to visualise your plot
p2 = 
  
  
  ## -----------------------------------------------------------------------------------------------------------------------------
ggplotly(p2)


#compute the basic reproduction number for the SITS model
R0_sits <- # write the formular
  


## -----------------------------------------------------------------------------------------------------------------------------
#| label: SIS


# Implementation of a prophylaxis intervention
