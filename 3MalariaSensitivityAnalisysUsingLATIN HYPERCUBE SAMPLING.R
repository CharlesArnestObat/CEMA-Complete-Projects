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


# ===========================================================
# === 3) GLOBAL SENSITIVITY: LATIN HYPERCUBE SAMPLING =====
# ===========================================================

# GLOBAL SENSITIVITY asks: "If I vary ALL parameters simultaneously
# across their plausible ranges, which ones matter most for the outcome?"
#
# LATIN HYPERCUBE SAMPLING (LHS) is an efficient way to sample parameter space:
# - Divides each parameter range into N equal intervals
# - Ensures we sample once from each interval
# - Creates good coverage of parameter space with fewer samples than full grid
#
# PRCC (Partial Rank Correlation Coefficient) measures monotonic relationships
# while accounting for other parameters (better than simple correlation)

# Define plausible ranges for parameters based on literature/expert knowledge
# These represent the uncertainty we have about "true" parameter values
ranges <- data.frame(
  param = c("a",   "b",   "c",   "m",  "gamma", "mu_v", "sigma"),
  min   = c(0.2,   0.2,   0.2,   2,    1/21,    1/14,   1/14),    # lower bounds
  max   = c(0.5,   0.5,   0.5,   10,   1/7,     1/7,    1/7),     # upper bounds
  stringsAsFactors = FALSE
)

# Add interpretable descriptions
ranges$description <- c(
  "Biting rate (bites/mosquito/day)",
  "Vector→human transmission probability",
  "Human→vector transmission probability",
  "Mosquito:human ratio",
  "Human recovery rate (1/days)",
  "Mosquito death rate (1/days)",
  "Mosquito incubation rate (1/days)"
)

cat("PARAMETER RANGES FOR GLOBAL ANALYSIS:\n")
for(i in 1:nrow(ranges)) {
  cat(sprintf("  %-8s: [%.3f, %.3f] - %s\n",
              ranges$param[i], ranges$min[i], ranges$max[i], ranges$description[i]))
}

# Number of LHS samples (trade-off between accuracy and computation time)
nLHS <- 500
cat(sprintf("\nGenerating %d Latin Hypercube samples...\n", nLHS))

# Step 1: Generate LHS design in unit hypercube [0,1]^d
set.seed(123)  # for reproducibility in class
lhs_unit <- randomLHS(nLHS, nrow(ranges))

# Step 2: Transform from [0,1] to actual parameter ranges
lhs_vals <- sapply(1:nrow(ranges), function(i) {
  ranges$min[i] + lhs_unit[, i] * (ranges$max[i] - ranges$min[i])
})
colnames(lhs_vals) <- ranges$param
lhs_df <- as.data.frame(lhs_vals)

cat("Sample of LHS design (first 5 rows):\n")
print(head(lhs_df, 5))

# Step 3: Evaluate model for each parameter combination
cat(sprintf("\nEvaluating model for %d parameter combinations...\n", nLHS))

# Show progress for long computations
pb <- txtProgressBar(min = 0, max = nLHS, style = 3)

Y_lhs <- numeric(nLHS)
for(i in 1:nLHS) {
  # Create parameter list by updating baseline with sampled values
  p <- par_base
  for(param_name in ranges$param) {
    p[[param_name]] <- lhs_df[i, param_name]
  }
  Y_lhs[i] <- run_model(p)
  setTxtProgressBar(pb, i)
}
close(pb)

lhs_df$Y <- Y_lhs  # add outcomes to data frame

cat(sprintf("\nGLOBAL ANALYSIS RESULTS:\n"))
cat(sprintf("  Prevalence range: [%.4f, %.4f]\n", min(Y_lhs), max(Y_lhs)))
cat(sprintf("  Mean prevalence: %.4f\n", mean(Y_lhs)))
cat(sprintf("  Std deviation: %.4f\n", sd(Y_lhs)))

# Calculate PRCC (Partial Rank Correlation Coefficients)
# This measures monotonic relationships while controlling for other variables
cat(sprintf("\nCalculating PRCC values...\n"))

prcc_results <- sapply(ranges$param, function(par_name) {
  # Use Spearman correlation with ranks (handles non-linear monotonic relationships)
  suppressWarnings(cor(rank(lhs_df[[par_name]]), rank(lhs_df$Y), method = "spearman"))
})

prcc_df <- data.frame(
  param = ranges$param,
  PRCC = as.numeric(prcc_results),
  abs_PRCC = abs(as.numeric(prcc_results)),
  stringsAsFactors = FALSE
)

# Add parameter descriptions
prcc_df$description <- ranges$description[match(prcc_df$param, ranges$param)]

# Sort by absolute PRCC (strongest relationships first)
prcc_df <- prcc_df[order(prcc_df$abs_PRCC, decreasing = TRUE), ]

cat("\nPRCC RESULTS (ranked by strength of relationship):\n")
print(prcc_df[, c("param", "PRCC", "description")])

# Create scatter plots for top 4 most important parameters
top_params <- head(prcc_df$param, 4)

library(gridExtra)  # for arranging multiple plots
plot_list <- lapply(top_params, function(par_name) {
  prcc_val <- prcc_df$PRCC[prcc_df$param == par_name]
  ggplot(lhs_df, aes_string(x = par_name, y = "Y")) +
    geom_point(alpha = 0.5, color = "steelblue") +
    geom_smooth(method = "loess", se = TRUE, color = "red") +
    labs(title = sprintf("%s (PRCC = %.3f)", par_name, prcc_val),
         x = par_name, y = "Prevalence") +
    theme_minimal()
})

# Arrange plots in 2x2 grid
grid_plot <- do.call(grid.arrange, c(plot_list, ncol = 2))

# Main PRCC plot
p3 <- ggplot(prcc_df, aes(x = reorder(param, abs_PRCC),
                          y = PRCC,
                          fill = PRCC > 0)) +
  geom_col(alpha = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("red", "blue"),
                    labels = c("Negative", "Positive"),
                    name = "Relationship") +
  labs(title = "Global Sensitivity: Latin Hypercube Sampling + PRCC",
       subtitle = sprintf("Based on %d samples across parameter ranges", nLHS),
       x = "Parameter",
       y = "PRCC (Partial Rank Correlation Coefficient)",
       caption = "Higher |PRCC| = stronger influence on prevalence") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold")) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)

print(p3)

cat(sprintf("\nMOST INFLUENTIAL PARAMETER (Global): %s (PRCC = %.3f)\n",
            prcc_df$param[1], prcc_df$PRCC[1]))

# ===========================================================
# ================== EXERCISE =====================
# ===========================================================


#"TASK 1: PARAMETER EXPLORATION#

#"Your turn to explore the model! Try these exercises:\n\n"#

# Change the baseline value of parameter 'b' (vector→human transmission)#
# from 0.3 to 0.1 and re-run the OAT analysis.#
#   Question: Does this change which parameter is most important?#

#B) Modify the ranges in the global analysis
#- Double the upper bound for 'm' (mosquito:human ratio) from 10 to 20
# - Reduce the lower bound for 'gamma' from 1/21 to 1/30


#"TASK 2: NEW INTERVENTION SCENARIO

#"Design a new intervention scenario

#"C) Create a 'Combined intervention' that starts at day 90 and
#   - Reduces biting rate 'a' by 40% (to 60% of original)
#   - Reduces transmission probability 'b' by 25% (to 75% of original)
#   - Increases mosquito mortality 'mu_v' by 50% (to 150% of original)

#  Compare this to the early and late ITN scenarios.
#  Question: Is the combined intervention more effective?

#TASK 3: SENSITIVITY OF INTERVENTIONS
#D) Choose the most important parameter from your global analysis.
#   Test how sensitive your intervention effectiveness is to uncertainty
#   in this parameter by running the intervention with this parameter
#   set to its minimum and maximum plausible values.

#"HINTS:\n")
#- Copy and modify the existing code sections
#- Use the same helper functions (run_model, simulate_with_intervention)
#- Follow the same plotting patterns for visualization
#- Don't forget to update parameter lists with par_base as the starting point


