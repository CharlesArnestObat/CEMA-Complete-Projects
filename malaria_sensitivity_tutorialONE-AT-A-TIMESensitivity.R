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

# ===========================================================
# ======= 2) LOCAL SENSITIVITY: ELASTICITY ANALYSIS =======
# ===========================================================

# ELASTICITY measures proportional sensitivity:
# "If parameter X increases by 1%, by what % does outcome Y change?"
#
# Formula: Elasticity = (∂Y/∂X) × (X/Y)
# where ∂Y/∂X is the derivative of Y with respect to X
#
# INTERPRETATION:
#   |Elasticity| > 1: outcome changes MORE than proportionally (elastic)
#   |Elasticity| < 1: outcome changes LESS than proportionally (inelastic)
#   Elasticity > 0: positive relationship (parameter ↑ → outcome ↑)
#   Elasticity < 0: negative relationship (parameter ↑ → outcome ↓)

cat("Calculating elasticities using numerical derivatives...\n\n")

elasticity_df <- do.call(rbind, lapply(params_to_vary, function(par_name) {
  base <- par_base
  Y <- run_model(base)                               # baseline outcome

  # Calculate numerical derivative using small perturbation
  h <- base[[par_name]] * 1e-4 + 1e-8               # tiny step size (0.01% of parameter)
  p_plus_h <- base
  p_plus_h[[par_name]] <- base[[par_name]] + h      # increase parameter by h

  # Finite-difference approximation: (Y(X+h) - Y(X)) / h
  dYdX <- (run_model(p_plus_h) - Y) / h

  # Elasticity formula: scale derivative by X/Y
  elasticity <- dYdX * (base[[par_name]] / Y)

  cat(sprintf("Parameter %s: baseline Y=%.4f, dY/dX=%.6f, elasticity=%.3f\n",
              par_name, Y, dYdX, elasticity))

  data.frame(param = par_name,
             Y = Y,
             dYdX = dYdX,
             elasticity = elasticity,
             abs_elasticity = abs(elasticity))
}))

# Sort by absolute elasticity
elasticity_df <- elasticity_df[order(elasticity_df$abs_elasticity, decreasing = TRUE), ]

cat("\nELASTICITY RESULTS (ranked by absolute elasticity):\n")
print(elasticity_df[, c("param", "elasticity", "abs_elasticity")])

# Plot elasticities with color coding for positive/negative
p2 <- ggplot(elasticity_df, aes(x = reorder(param, abs_elasticity),
                                y = elasticity,
                                fill = elasticity > 0)) +
  geom_col(alpha = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("red", "blue"),
                    labels = c("Negative", "Positive"),
                    name = "Relationship") +
  labs(title = "Elasticity Analysis: Proportional Sensitivity",
       subtitle = "% change in prevalence per 1% change in parameter",
       x = "Parameter",
       y = "Elasticity",
       caption = "Blue = positive relationship, Red = negative relationship") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold")) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)

print(p2)

cat(sprintf("\nMOST ELASTIC PARAMETER: %s (elasticity = %.3f)\n",
            elasticity_df$param[1], elasticity_df$elasticity[1]))

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
# == 4) VARIANCE-BASED GLOBAL SENSITIVITY: FAST ===========
# ===========================================================
cat("\n", "=" %R% 60, "\n")
cat("PART 4: VARIANCE-BASED SENSITIVITY ANALYSIS (FAST)\n")
cat("=" %R% 60, "\n\n")

# VARIANCE-BASED METHODS ask: "What fraction of output variance
# is caused by uncertainty in each input parameter?"
#
# FAST (Fourier Amplitude Sensitivity Test) decomposes total variance into:
# - First-order effects: variance due to parameter Xi alone
# - Higher-order effects: variance due to interactions between parameters
#
# ADVANTAGE: Captures parameter interactions
# DISADVANTAGE: Requires more model evaluations than correlation methods

cat("Setting up FAST analysis...\n")

# Create wrapper function that accepts parameter matrix
# (required format for sensitivity package)
fast_model <- function(X) {
  # X is a matrix where each row is a parameter combination
  apply(X, 1, function(row) {
    p <- par_base
    names(row) <- ranges$param
    for(k in names(row)) p[[k]] <- row[[k]]
    run_model(p)
  })
}

cat("Generating FAST sampling design...\n")
set.seed(456)  # different seed for FAST

# Create FAST sampling design
# This uses special frequency-based sampling to efficiently estimate variance components
fast_setup <- fast99(
  model = NULL,                    # we'll evaluate separately
  factors = ranges$param,          # parameter names
  n = 200,                         # base sample size (total samples ≈ n × n_params)
  q = "qunif",                     # uniform distribution
  q.arg = lapply(1:nrow(ranges), function(i) {
    list(min = ranges$min[i], max = ranges$max[i])
  })
)

cat(sprintf("Evaluating model at %d FAST design points...\n", nrow(fast_setup$X)))

# Evaluate model at FAST design points
fast_setup$y <- fast_model(fast_setup$X)

cat("Computing FAST sensitivity indices...\n")

# Extract FAST indices
fast_results <- tell(fast_setup)
first_order <- fast_results$D1 / fast_results$V  # first-order indices (normalized)
total_order <- fast_results$Dt / fast_results$V  # total-order indices (normalized)

# Create results data frame
fast_df <- data.frame(
  param = ranges$param,
  first_order = first_order,
  total_order = total_order,
  interaction = total_order - first_order,  # interaction effects
  stringsAsFactors = FALSE
)

# Sort by first-order effects
fast_df <- fast_df[order(fast_df$first_order, decreasing = TRUE), ]

cat(sprintf("\nVARIANCE DECOMPOSITION:\n"))
cat(sprintf("  Total output variance: %.6f\n", fast_results$V))
cat(sprintf("  Sum of first-order effects: %.3f (%.1f%%)\n",
            sum(first_order), sum(first_order) * 100))

cat("\nFAST RESULTS (ranked by first-order effects):\n")
print(fast_df)

# Plot FAST results
library(reshape2)  # for melting data frame

fast_plot_df <- melt(fast_df[, c("param", "first_order", "interaction")],
                     id.vars = "param",
                     variable.name = "effect_type",
                     value.name = "variance_fraction")

p4 <- ggplot(fast_plot_df, aes(x = reorder(param, variance_fraction),
                               y = variance_fraction,
                               fill = effect_type)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("steelblue", "orange"),
                    labels = c("First-order", "Interactions"),
                    name = "Effect Type") +
  labs(title = "Variance-Based Sensitivity Analysis (FAST)",
       subtitle = "Fraction of output variance explained by each parameter",
       x = "Parameter",
       y = "Variance Fraction",
       caption = "First-order = individual parameter effects, Interactions = combined effects") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))

print(p4)

cat(sprintf("\nMOST IMPORTANT PARAMETER (Variance): %s (%.1f%% of variance)\n",
            fast_df$param[1], fast_df$first_order[1] * 100))

# ===========================================================
# ============= 5) SCENARIO ANALYSIS =======================
# ===========================================================
cat("\n", "=" %R% 60, "\n")
cat("PART 5: SCENARIO ANALYSIS - INTERVENTION TIMING\n")
cat("=" %R% 60, "\n\n")

# SCENARIO ANALYSIS asks: "What if we implement a specific intervention?
# How does the timing of implementation affect outcomes?"
#
# Here we simulate distributing Insecticide-Treated bed Nets (ITNs):
# - Reduces mosquito biting rate (parameter 'a')
# - Increases mosquito mortality (parameter 'mu_v')
# We compare early vs late rollout

cat("INTERVENTION: Insecticide-Treated bed Nets (ITNs)\n")
cat("Effects: Reduce biting rate by 30%, increase mosquito mortality by 20%\n")
cat("Scenarios: No intervention, Early rollout (day 60), Late rollout (day 180)\n\n")

# Function to simulate intervention starting at specified time
simulate_with_intervention <- function(t_start = 120,
                                       a_reduction = 0.7,      # reduce to 70% of original
                                       mu_increase = 1.2) {    # increase by 20%

  # Create time-dependent ODE function
  intervention_ode <- function(t, y, p) {
    # Modify parameters after intervention starts
    p_modified <- p
    if (t >= t_start) {
      p_modified$a    <- p$a * a_reduction     # reduce biting
      p_modified$mu_v <- p$mu_v * mu_increase  # increase mortality
    }
    # Call original model with modified parameters
    malaria_ode(t, y, p_modified)
  }

  # Solve ODEs with time-dependent intervention
  ode(y = y0, times = times, func = intervention_ode, parms = par_base, method = "rk4")
}

cat("Running scenario simulations...\n")

# Run three scenarios
baseline_sim  <- as.data.frame(ode(y0, times, malaria_ode, par_base, method = "rk4"))
early_sim     <- as.data.frame(simulate_with_intervention(t_start = 60))
late_sim      <- as.data.frame(simulate_with_intervention(t_start = 180))

# Calculate prevalence time series for each scenario
scenario_df <- rbind(
  data.frame(
    time = baseline_sim$time,
    prevalence = baseline_sim$Ih / par_base$Nh,
    scenario = "No intervention",
    stringsAsFactors = FALSE
  ),
  data.frame(
    time = early_sim$time,
    prevalence = early_sim$Ih / par_base$Nh,
    scenario = "Early ITNs (day 60)",
    stringsAsFactors = FALSE
  ),
  data.frame(
    time = late_sim$time,
    prevalence = late_sim$Ih / par_base$Nh,
    scenario = "Late ITNs (day 180)",
    stringsAsFactors = FALSE
  )
)

# Calculate summary statistics
scenario_summary <- scenario_df %>%
  group_by(scenario) %>%
  summarise(
    final_prevalence = tail(prevalence, 1),
    mean_prevalence = mean(prevalence),
    max_prevalence = max(prevalence),
    .groups = 'drop'
  )

cat("\nSCENARIO RESULTS:\n")
print(scenario_summary)

# Calculate intervention impact
baseline_final <- scenario_summary$final_prevalence[scenario_summary$scenario == "No intervention"]
early_reduction <- (baseline_final - scenario_summary$final_prevalence[scenario_summary$scenario == "Early ITNs (day 60)"]) / baseline_final * 100
late_reduction <- (baseline_final - scenario_summary$final_prevalence[scenario_summary$scenario == "Late ITNs (day 180)"]) / baseline_final * 100

cat(sprintf("\nINTERVENTION IMPACT:\n"))
cat(sprintf("  Early ITNs: %.1f%% reduction in final prevalence\n", early_reduction))
cat(sprintf("  Late ITNs:  %.1f%% reduction in final prevalence\n", late_reduction))
cat(sprintf("  Benefit of early vs late: %.1f percentage points\n", early_reduction - late_reduction))

# Plot scenario comparison
p5 <- ggplot(scenario_df, aes(x = time, y = prevalence, color = scenario)) +
  geom_line(size = 1.2) +
  geom_vline(xintercept = c(60, 180), linetype = "dashed", alpha = 0.5) +
  annotate("text", x = 60, y = max(scenario_df$prevalence) * 0.9,
           label = "Early ITNs", angle = 90, hjust = 1, size = 3) +
  annotate("text", x = 180, y = max(scenario_df$prevalence) * 0.9,
           label = "Late ITNs", angle = 90, hjust = 1, size = 3) +
  scale_color_manual(values = c("red", "blue", "green")) +
  labs(title = "Scenario Analysis: Timing of Intervention Rollout",
       subtitle = "ITNs reduce biting rate by 30% and increase mosquito mortality by 20%",
       x = "Time (days)",
       y = "Infectious Prevalence",
       color = "Scenario",
       caption = "Dashed lines show intervention start times") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom")

print(p5)

# ===========================================================
# ============= SUMMARY AND COMPARISON =====================
# ===========================================================
cat("\n", "=" %R% 70, "\n")
cat("SUMMARY: COMPARISON OF SENSITIVITY METHODS\n")
cat("=" %R% 70, "\n\n")

# Create comparison table of most important parameter from each method
comparison_df <- data.frame(
  Method = c("OAT (±10%)", "Elasticity", "Global (PRCC)", "Variance (FAST)"),
  Most_Important = c(
    oat_df$param[1],
    elasticity_df$param[1],
    prcc_df$param[1],
    fast_df$param[1]
  ),
  Measure = c(
    sprintf("Δ = %.4f", oat_df$delta[1]),
    sprintf("E = %.3f", elasticity_df$elasticity[1]),
    sprintf("PRCC = %.3f", prcc_df$PRCC[1]),
    sprintf("S₁ = %.3f", fast_df$first_order[1])
  ),
  Interpretation = c(
    "Largest change in prevalence",
    "Most elastic response",
    "Strongest rank correlation",
    "Explains most variance"
  ),
  stringsAsFactors = FALSE
)

cat("COMPARISON OF SENSITIVITY METHODS:\n")
print(comparison_df)

cat("\nKEY INSIGHTS:\n")
cat("1. Local vs Global: Local methods may miss important global effects\n")
cat("2. Parameter interactions: Variance-based methods capture interaction effects\n")
cat("3. Method choice: Depends on research question and computational budget\n")
cat("4. Intervention timing: Early intervention is more effective than late\n")

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


