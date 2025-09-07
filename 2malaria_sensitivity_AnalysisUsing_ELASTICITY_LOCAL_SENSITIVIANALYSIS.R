#############################################
## MALARIA MODEL + SENSITIVITY ANALYSIS   ##
## Tutorial for beginners       ##
##                                         ##
## Learning Objectives:                    ##
## 1. Understand local vs global methods  ##
## 2. Apply One at a time, elasticity, LHS, FAST techniques ##
## 3. Interpret sensitivity results       ##
## 4. Conduct scenario analysis           ##
#############################################

# Clear workspace (good practice at start of scripts)
rm(list = ls())

# -------------------- INSTALL & LOAD PACKAGES --------------------
# First, we check if packages are installed, install if missing


library(deSolve) # for solving differential equations (ODE)
library(lhs) # for latin hypercube sampling
library(sensitivity) # For FAST/Sobol variance-based sensitivity
library(tidyverse) # for plotting (ggplot2)


# if packages are not installed, use install.packages()

# -------------------- UNDERSTANDING THE MALARIA MODEL --------------------
# Before diving into sensitivity analysis, let's understand what we're modeling!
#
# This is a Ross-Macdonald malaria transmission model with:
#
# HUMANS (SIS-like dynamics):
#   - Susceptible humans get infected by mosquito bites
#   - Infectious humans recover (but can be reinfected)
#   - We only track Ih = infectious humans (Sh calculated as Nh - Ih)
#
# MOSQUITOES (SEI dynamics):
#   - Susceptible mosquitoes (Sv) bite infectious humans and become exposed (Ev)
#   - Exposed mosquitoes (Ev) incubate the parasite, then become infectious (Iv)
#   - All mosquitoes die at rate mu_v
#
# KEY INSIGHT: Malaria requires BOTH infected humans AND infected mosquitoes!
# This creates interesting parameter interactions we'll explore.

# -------------------- MODEL DEFINITION --------------------
malaria_ode <- function(t, y, p) {
  # The 'with' function lets us use variable names directly instead of p$a, p$b, etc.
  with(as.list(c(y, p)), {
    
    # STEP 1: Calculate derived quantities
    # Total populations (we assume these are fixed)
    Nh <- Nh                    # total humans (given as parameter)
    Nv <- m * Nh                # total mosquitoes = m times number of humans
    
    # Susceptible populations (calculated from totals minus infected)
    Sh <- Nh - Ih               # susceptible humans
    Sv <- Nv - Ev - Iv          # susceptible mosquitoes
    
    # STEP 2: Calculate forces of infection (transmission rates)
    # These represent the "pressure" of infection in each direction
    
    # Rate humans get infected = biting rate × transmission probability × infected mosquito proportion
    lambda_h <- a * b * (Iv / Nv)
    
    # Rate mosquitoes get infected = biting rate × transmission probability × infected human proportion
    lambda_v <- a * c * (Ih / Nh)
    
    # STEP 3: Write the differential equations
    # Each equation represents "inflows - outflows" for each compartment
    
    # Humans: gain infections from mosquitoes, lose to recovery
    dIh <- lambda_h * Sh - gamma * Ih
    
    # Exposed mosquitoes: gain from biting infected humans, lose to incubation or death
    dEv <- lambda_v * Sv - (sigma + mu_v) * Ev
    
    # Infectious mosquitoes: gain from incubation, lose to death
    dIv <- sigma * Ev - mu_v * Iv
    
    # Return the derivatives (required format for deSolve package)
    return(list(c(dIh, dEv, dIv)))
  })
}

# -------------------- BASELINE PARAMETER VALUES --------------------
# These represent a "typical" malaria setting - we'll test how sensitive
# our results are to changes in these values

par_base <- list(
  a = 0.3,      # mosquito biting rate (bites per mosquito per day)
  b = 0.3,      # probability of transmission from mosquito to human per bite
  c = 0.3,      # probability of transmission from human to mosquito per bite
  m = 5,        # mosquito-to-human ratio (5 mosquitoes per person)
  gamma = 1/14, # human recovery rate (recover in ~14 days on average)
  mu_v = 1/10,  # mosquito death rate (live ~10 days on average)
  sigma = 1/10, # mosquito incubation rate (~10 days to become infectious)
  Nh = 1e5      # total human population (100,000 people)
)

# Print parameters in a readable format
for(i in names(par_base)) {
  cat(sprintf("  %-8s = %g\n", i, par_base[[i]]))
}

# -------------------- INITIAL CONDITIONS --------------------
# Starting values for our state variables at time t=0
# We begin with small numbers of infections in both humans and mosquitoes

y0 <- c(
  Ih = 100,     # 100 infectious humans (0.1% of population)
  Ev = 100,     # 100 exposed mosquitoes
  Iv = 100      # 100 infectious mosquitoes
)


# Time points where we want to evaluate the model (daily for 1 year)
times <- seq(0, 365, by = 1)

# -------------------- HELPER FUNCTION: RUN MODEL --------------------
# This function runs our malaria model and returns a single summary measure:
# the endemic prevalence (% of humans infected) at the end of the simulation.
# We use this as our "outcome of interest" for sensitivity analysis.

run_model <- function(p) {
  # Solve the differential equation using Runge-Kutta method
  out <- ode(y = y0, times = times, func = malaria_ode, parms = p, method = "rk4")
  
  # Convert to data frame for easier handling
  df <- as.data.frame(out)
  
  # Calculate final prevalence = infectious humans / total humans
  prev_end <- tail(df$Ih, 1) / p$Nh
  
  return(prev_end)
}

# Test our baseline model
baseline_prevalence <- run_model(par_base)

#===================================
# ========== 1) LOCAL SENSITIVITY: ONE-AT-A-TIME ==========
# ===========================================================

cat("PART 1: LOCAL SENSITIVITY ANALYSIS (One-At-a-Time)\n")


# LOCAL SENSITIVITY asks: "If I change ONE parameter by a small amount
# while keeping everything else fixed, how much does my outcome change?"
#
# This is like asking: "What happens if mosquito biting rate increases by 10%
# but everything else stays the same?"
#
# ADVANTAGE: Easy to understand and compute
# DISADVANTAGE: Misses parameter interactions, only local around baseline

# Choose which parameters to test (we'll skip b and c to keep it manageable)
params_to_vary <- c("a", "m", "gamma", "mu_v", "sigma")


# For each parameter, calculate outcomes at -10%, baseline, and +10%
oat_results <- lapply(params_to_vary, function(par_name) {
  cat("Testing parameter:", par_name, "\n")
  
  # Start with baseline parameters
  base <- par_base
  baseline_outcome <- run_model(base)
  
  # Create +10% and -10% versions
  plus_10  <- base
  minus_10 <- base
  plus_10[[par_name]]  <- base[[par_name]] * 1.10  # increase by 10%
  minus_10[[par_name]] <- base[[par_name]] * 0.90  # decrease by 10%
  
  # Calculate outcomes for perturbed parameters
  plus_outcome  <- run_model(plus_10)
  minus_outcome <- run_model(minus_10)
  
  cat(sprintf("  Baseline: %.4f, -10%%: %.4f, +10%%: %.4f\n",
              baseline_outcome, minus_outcome, plus_outcome))
  
  # Return results as a named vector (will be combined into data frame)
  c(param = par_name,
    baseline = baseline_outcome,
    minus10 = minus_outcome,
    plus10 = plus_outcome)
})

# Convert list of results into a data frame
oat_df <- as.data.frame(do.call(rbind, oat_results))
# Convert character columns to numeric (except parameter names)
oat_df[, 2:4] <- lapply(oat_df[, 2:4], as.numeric)

# Calculate the "sensitivity measure": difference between +10% and -10% outcomes
# Larger differences = more sensitive to that parameter
oat_df$delta <- oat_df$plus10 - oat_df$minus10
oat_df$abs_delta <- abs(oat_df$delta)

# Sort by absolute sensitivity (most sensitive first)
oat_df <- oat_df[order(oat_df$abs_delta, decreasing = TRUE), ]

cat("\nONE-AT-A-TIME RESULTS (ranked by sensitivity):\n")
print(oat_df[, c("param", "baseline", "delta", "abs_delta")])

# Create a tornado plot (horizontal bar chart)
# This is called "tornado" because sensitive parameters create wide "tornado" shapes
p1 <- ggplot(oat_df, aes(x = reorder(param, abs_delta), y = delta)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  coord_flip() +
  labs(title = "Local Sensitivity Analysis (One-At-a-Time, ±10%)",
       subtitle = "Change in endemic prevalence when parameter changes by ±10%",
       x = "Parameter",
       y = "Δ Prevalence (plus10% − minus10%)",
       caption = "Larger bars = more sensitive parameters") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))

print(p1)

cat(sprintf("\nMOST SENSITIVE PARAMETER: %s (Δ = %.4f)\n",
            oat_df$param[1], oat_df$delta[1]))
cat(sprintf("LEAST SENSITIVE PARAMETER: %s (Δ = %.4f)\n\n",
            tail(oat_df$param, 1), tail(oat_df$delta, 1)))