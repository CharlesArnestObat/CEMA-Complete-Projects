library(tidyverse)
library(deSolve)
observed_data=read.csv("sis_si_weekly_cases.csv")
observed_data
observed_data <- read.csv("sis_si_weekly_cases.csv") %>%
  transmute(time = as.integer(week),
            number_infected_individuals = as.numeric(observed_cases))

head(observed_data, 8)
# Quick visualization
#Task 1: Describe the overall pattern. Does the series look stationary over the year?
gg0 <- ggplot() +
  geom_point(data = observed_data,
             aes(x = time, y = number_infected_individuals),
             color = "black", shape = 4, size = 1, stroke = 1.2) +
  geom_line(data = observed_data,
            aes(x = time, y = number_infected_individuals),
            color = "black") +
  labs(x = "Time (weeks)", y = "Number of Infected Individuals") +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text  = element_text(size = 10, color = "black")) +
  scale_x_continuous(limits = c(0, 55), breaks = seq(0, 55, by = 5)) +
  scale_y_continuous(limits = c(0, max(observed_data$number_infected_individuals) * 1.2),
                     breaks = scales::pretty_breaks())

gg0
#3. SIS–SI model with cumulative infections
#We follow the provided structure and add Ch (cumulative human infections) to compute weekly incidence.
# Baseline parameters
biting      <- 0.3
seeking     <- 0.1
infectivity <- 0.3
alpha       <- biting * seeking * infectivity
k           <- 5

params_base <- list(
  alpha = alpha,
  beta  = 0.5,
  gamma = 1/25,   # 0.04 day^-1
  mu_m  = 1/10,   # 0.1 day^-1
  Nh    = 1000,
  Nm    = 1000 * k
)

# Initial state (+ cumulative infections)
state0 <- c(
  Sh = params_base$Nh - 10,
  Ih = 10,
  Sm = params_base$Nm - 0.02 * params_base$Nm,
  Im = 0.02 * params_base$Nm,
  Ch = 0
)

# ODE
sis_si <- function(t, state, p){
  with(as.list(c(state, p)), {
    Nh <- Sh + Ih
    Nm <- Sm + Im
    
    # Humans
    dSh <- -(alpha * (Im / Nh) * Sh) + gamma * Ih
    dIh <-  (alpha * (Im / Nh) * Sh) - gamma * Ih
    
    # Mosquitoes
    births <- mu_m * Nm
    dSm <- births - (beta * (Ih / Nh) * Sm) - mu_m * Sm
    dIm <- (beta * (Ih / Nh) * Sm) - mu_m * Im
    
    # Cumulative human infections
    dCh <- alpha * (Im / Nh) * Sh
    
    list(c(dSh, dIh, dSm, dIm, dCh))
  })
}
#3.1 Simulator → weekly incidence
#We introduce a transmission scale factor alpha_scale so that alpha ← alpha * alpha_scale and compute weekly incidence from Ch.


simulate_weekly_incidence <- function(alpha_scale = 1,
                                      params = params_base,
                                      state = state0){
  p <- params
  p$alpha <- params$alpha * alpha_scale
  
  times <- seq(0, 52*7, by = 1)  # days across 52 weeks
  out <- as.data.frame(ode(y = state, times = times, func = sis_si, parms = p))
  
  idx <- seq(1, nrow(out), by = 7)        # day 0,7,...,364
  weekly_inc <- diff(out$Ch[idx])
  
  tibble(time = 1:52,
         model_weekly_incidence = pmax(weekly_inc, 0))
}

# Quick overlay (alpha_scale=1, rho_scale=1)
alpha_scale=1
rho_scale=1
pred0 <- simulate_weekly_incidence(alpha_scale) %>%
  mutate(pred_cases = rho_scale * model_weekly_incidence)

gg0 + geom_line(data = pred0, aes(x = time, y = pred_cases),
                color = "steelblue", linewidth = 1)
#Task 2: Try alpha_scale = 0.8 and 1.2. What happens to the curve?
#
#4. Least Squares (sse) fit
#Estimate alpha_scale and rho_scale by minimizing the sse between predicted and observed weekly cases.

sse <- function(x, y) sqrt(mean((x - y)^2))

objective_sse <- function(theta_log){
  alpha_scale <- exp(theta_log[1])  # >0
  rho_scale   <- exp(theta_log[2])  # >0
  
  pred <- simulate_weekly_incidence(alpha_scale) %>%
    mutate(pred_cases = rho_scale * model_weekly_incidence)
  
  df <- observed_data %>% left_join(pred, by = "time")
  res <- sse(df$pred_cases, df$number_infected_individuals)
  return(res)
}

theta0_log <- log(c(alpha_scale = 1, rho_scale = 1))
fit_sse <- optim(par = theta0_log, fn = objective_sse, method = "BFGS",
                 control = list(maxit = 1000, reltol = 1e-10))

alpha_scale_sse <- exp(fit_sse$par[1])
rho_scale_sse   <- exp(fit_sse$par[2])
sse_hat       <- fit_sse$value

cat(sprintf("sse fit: alpha_scale = %.3f, rho_scale = %.3f, sse = %.2f\n",
            alpha_scale_sse, rho_scale_sse, sse_hat))
## sse fit: alpha_scale = 0.946, rho_scale = 0.877, sse = 27.75
pred_sse <- simulate_weekly_incidence(alpha_scale_sse) %>%
  mutate(pred_cases = rho_scale_sse * model_weekly_incidence)

gg0 + geom_line(data = pred_sse, aes(x = time, y = pred_cases),
                color = "dodgerblue4", linewidth = 1) +
  labs(title = paste0("sse:   alpha:", round(alpha_scale_sse,3), ",  rho:", round(rho_scale_sse,3), ",  sse:", round(sse_hat,2)))

#4.1. Cumulative comparison
df_plot <- observed_data %>%
  left_join(pred_sse, by = "time") %>%
  transmute(time,
            cum_obs  = cumsum(number_infected_individuals),
            cum_pred = cumsum(pred_cases))

ggplot(df_plot, aes(x = time)) +
  geom_line(aes(y = cum_obs),  linewidth = 1, linetype = "dashed") +
  geom_line(aes(y = cum_pred), linewidth = 1, color = "darkgreen") +
  labs(title = paste0("sse - cumul:   alpha:", round(alpha_scale_sse,3), ",  rho:", round(rho_scale_sse,3), ",  sse:", round(sse_hat,2))) +
  theme_minimal()

#5. Maximum Likelihood (Negative Binomial)
#Assume yw∼NB(μw,size)
#with μw=ρscale×incidencew.
#This minimal MLE uses nlminb on the log-scale with simple bounds.

y <- observed_data$number_infected_individuals
y

# NegBin negative log-likelihood
nll_nb <- function(y, mu, size) {
  mu <- pmax(mu, 1e-12)
  -sum(dnbinom(y, mu = mu, size = size, log = TRUE))
}

# Method-of-moments starting values
inc1  <- simulate_weekly_incidence(alpha_scale = 1)$model_weekly_incidence
rho0  <- sum(y) / sum(inc1)
m     <- mean(y); v <- var(y)
size0 <- if (v > m) (m^2)/(v - m) else 100
size0 <- min(max(size0, 1), 1e3)

theta0_log <- log(c(alpha_scale = 1, rho_scale = rho0, size = size0))

# Log-scale objective
objective_log <- function(theta_log){
  alpha_scale <- exp(theta_log[1])
  rho_scale   <- exp(theta_log[2])
  size        <- exp(theta_log[3])
  
  inc <- simulate_weekly_incidence(alpha_scale)$model_weekly_incidence
  mu  <- rho_scale * inc
  nll_nb(y, mu, size)
}

# Bounded optimization (log-scale)
lower_log <- log(c(0.05, 1e-3, 0.5))
upper_log <- log(c(20,   100,  1e4))

fit_mle <- nlminb(start = theta0_log, objective = objective_log,
                  lower = lower_log, upper = upper_log,
                  control = list(eval.max = 2000, iter.max = 2000))

alpha_scale_mle <- exp(fit_mle$par[1])
rho_scale_mle   <- exp(fit_mle$par[2])
size_mle        <- exp(fit_mle$par[3])
nll_mle         <- fit_mle$objective
AIC_mle         <- 2*3 + 2*nll_mle

cat(sprintf("NB MLE: alpha=%.3f, rho=%.3f, size=%.2f, NLL=%.2f, AIC=%.1f, converged=%s\n",
            alpha_scale_mle, rho_scale_mle, size_mle, nll_mle, AIC_mle, fit_mle$convergence==0))
## NB MLE: alpha=0.877, rho=0.929, size=11.88, NLL=244.89, AIC=495.8, converged=TRUE
pred_mle <- simulate_weekly_incidence(alpha_scale_mle) %>%
  mutate(pred_cases = rho_scale_mle * model_weekly_incidence)

gg0 + geom_line(data = pred_mle, aes(x = time, y = pred_cases),
                color = "dodgerblue4", linewidth = 1) +
  labs(title = paste0("MLE:   alpha:", round(alpha_scale_sse,3), ",  rho:", round(rho_scale_mle,3)))
#Task 3 : Does sse capture the shape well? What does rho_scale control?
  
#  6. To go further
#alpha_scale scales transmission; rho_scale maps model incidence to reported cases; size controls overdispersion.
#With only case counts, alpha_scale and rho_scale can partially trade off; external data on reporting helps.
#Interpretation:
  
 # alpha_scale scales transmission;
#rho_scale maps model incidence to reported cases;
#size controls overdispersion (smaller → noisier counts).
#Identifiability: - With only case counts, alpha_scale and rho_scale can partially trade off; - external data (e.g., test positivity, care-seeking) helps.

#Task 4: Refit the model using MLE method with Poisson likelihood and compare AIC wit
