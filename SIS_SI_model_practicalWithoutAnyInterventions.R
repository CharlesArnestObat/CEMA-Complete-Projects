##############################
## SIS–SI vector–host model ##
##############################
library(deSolve)
library(tidyverse)
library(plotly)


# -------------------------
# Parameters 
# -------------------------
biting = 0.3
seeking = 0.1
infectivity = 0.3
alpha = biting * seeking * infectivity
k = 5 # mosquito to human ratio

params <- list(
  alpha     = alpha,      # transmission parameter from mosquitoes to humans
  beta = 0.5,      # P(human infection per bite from  an infectious mosquito)
  gamma = 1/25,    # human recovery rate
  mu_m  = 1/10,     # mosquito death rate (1/lifespan)
  Nh    = 1000,     # total number of humans and we assume that the system is closed
  Nm    = 1000*k       # mosquito:human ratio  (Nm = k * Nh)
)


# -------------------------
# Model state variables
# Sh, Ih : susceptible & infectious humans
# Sm, Im : susceptible & infectious mosquitoes
# -------------------------

state0 <- c(
  Sh = params$Nh - 10,  # 1% initially infectious
  Ih = 10,
  Sm = params$Nm - 0.02*params$Nm,  # 2% infectious mosquitoes
  Im = 0.02*params$Nm
)

# -------------------------
# ODE system
# -------------------------
sis_si <- function(t, state, p){
  with(as.list(c(state, p)), {
    Nh <- Sh + Ih # should remain constant 
    Nm <- Sm + Im  # should remain constant 
    
    # Humans 
    dSh <- -(alpha * (Im / Nh)   * Sh) + gamma * Ih
    dIh <-  (alpha * (Im / Nh)   * Sh) - gamma * Ih
    
    # Mosquitoes
    births  <- mu_m * Nm                   
    dSm <- births - (beta * (Ih / Nh)* Sm) - mu_m * Sm
    dIm <-  (beta * (Ih / Nh)* Sm) - mu_m * Im
    
    list(c(dSh, dIh, dSm, dIm))
  })
}

# -------------------------
# Solve
# -------------------------
times <- seq(0, 365*1, by = 1) # 1 year and daily time steps
out <- as.data.frame(ode(y = state0, times = times, func = sis_si, parms = params))



### plotting

subplot(
  plot_ly(out, x = ~time) %>%
    add_lines(y = ~Sh, name = "S_h") %>%
    add_lines(y = ~Ih, name = "I_h") %>%
    layout(title = "Human popn."),
  
  plot_ly(out, x = ~time) %>%
    add_lines(y = ~Sm, name = "S_m") %>%
    add_lines(y = ~Im, name = "I_m") %>%
    layout(title = "Mosquito popn."),
  
  nrows = 1
)
