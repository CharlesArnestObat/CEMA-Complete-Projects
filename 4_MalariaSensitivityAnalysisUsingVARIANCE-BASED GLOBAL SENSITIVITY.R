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


