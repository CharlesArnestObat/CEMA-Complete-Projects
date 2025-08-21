# Load package
library(deSolve)

# 1. Define model function
sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    N <- S + I + R
    dS <- -beta * S * I / N
    dI <- beta * S * I / N - gamma * I
    dR <- gamma * I
    
    list(c(dS, dI, dR))
  })
}

# 2. Initial values
init <- c(S = 999, I = 1, R = 0)

# 3. Parameters
params <- c(beta = 0.3, gamma = 0.1)

# 4. Time steps
times <- seq(0, 160, by = 1)

# 5. Solve ODE
out <- ode(y = init, times = times, func = sir_model, parms = params)

# 6. Convert to dataframe
out <- as.data.frame(out)

# 7. Plot results
matplot(out$time, out[,2:4], type = "l", lty = 1, col = c("blue", "red", "green"),
        xlab = "Time (days)", ylab = "Number of people")
legend("right", legend = c("Susceptible", "Infected", "Recovered"),
       col = c("blue", "red", "green"), lty = 1)