# Mixing Units Practical: Per-week vs Per-day (SIS with standard incidence)
# Goal: show how misinterpreting weekly quantities as daily rates distorts dynamics & inference.

library(deSolve)



# -------------------------
# 0) Helpers
# -------------------------

logit <- function(p) log(p/(1-p))
inv_logit <- function(x) 1/(1+exp(-x))
#inv_logit(2)
clamp01 <- function(x) pmin(pmax(x, 1e-9), 1-1e-9)

sis_rhs <- function(t, state, parms){
  I <- state["I"]
  with(as.list(parms), {
    dI <- beta * I * (N - I) / N - gamma * I
    list(c(dI))
  })
}

simulate_SIS <- function(beta, gamma, I0, N, times){
  out <- ode(y = c(I = I0), times = times, func = sis_rhs,
             parms = list(beta=beta, gamma=gamma, N=N), method="lsoda")
  df <- as.data.frame(out)
  df$p <- clamp01(df$I / N)
  df
}

binom_obs <- function(p, m){ rbinom(length(p), size=m, prob=p) }

# -------------------------
# 1) "Literature" quantity in weeks
# -------------------------
# Suppose a paper reports: weekly recovery PROBABILITY p_week = 0.7.
p_week <- 0.7

# Correct daily RATE (hazard) from weekly probability:
# p_week = 1 - exp(-gamma_day * 7)  =>  gamma_day = -log(1 - p_week)/7
gamma_day_true <- -log(1 - p_week) / 7

# Two WRONG interpretations to demonstrate:
gamma_day_naive_div <- p_week / 7     # linear approx: OK only if p_week << 1 (here it's not)
gamma_day_as_rate   <- p_week         # (very wrong) treating weekly probability as if it were a daily *rate*

# Choose an R0 and compute betas (standard incidence: R0 = beta/gamma)
R0 <- 3
beta_true        <- R0 * gamma_day_true
beta_naive_div   <- R0 * gamma_day_naive_div
beta_as_rate     <- R0 * gamma_day_as_rate

# Early growth rates r = beta - gamma
r_true       <- beta_true - gamma_day_true
r_naive_div  <- beta_naive_div - gamma_day_naive_div
r_as_rate    <- beta_as_rate - gamma_day_as_rate

doubling_time <- function(r){ if(r <= 0) Inf else log(2)/r }

# -------------------------
# 2) Simulate data with the CORRECT parameters
# -------------------------
set.seed(1)
N <- 10000
I0 <- 50
times_daily <- seq(0, 70, by=1)   # 10 weeks daily
times_week  <- seq(0, 70, by=7)
m_week      <- rep(200, length(times_week))

sim_true <- simulate_SIS(beta_true, gamma_day_true, I0, N, times_daily)

# Observations at weekly times
p_week_true <- sim_true$p[match(times_week, sim_true$time)]
y_week <- binom_obs(p_week_true, m_week)
dat <- data.frame(time = times_week, m = m_week, y = y_week)
dat

# -------------------------
# 3) Compute predicted weekly prevalence under three parameter sets
# -------------------------
pred_week <- function(beta, gamma){
  s <- simulate_SIS(beta, gamma, I0, N, times_week)
  s$p
}

p_true      <- pred_week(beta_true,      gamma_day_true)
p_naive_div <- pred_week(beta_naive_div, gamma_day_naive_div)
p_as_rate   <- pred_week(beta_as_rate,   gamma_day_as_rate)

# Negative log-likelihood under Binomial observation model
nll <- function(p, dat){
  p <- clamp01(p)
  -sum(dbinom(dat$y, size=dat$m, prob=p, log=TRUE))
}

nll_true      <- nll(p_true, dat)
nll_naive_div <- nll(p_naive_div, dat)
nll_as_rate   <- nll(p_as_rate, dat)

# -------------------------
# 4) Output & plots
# -------------------------
cat("=== Parameters derived from p_week = 0.7 ===\n")
cat(sprintf("Correct daily rate:        gamma = %.4f /day\n", gamma_day_true))
cat(sprintf("Naive divide-by-7 rate:   gamma = %.4f /day\n", gamma_day_naive_div))
cat(sprintf("Wrong (as a daily rate):  gamma = %.4f /day  <-- WRONG\n\n", gamma_day_as_rate))

cat(sprintf("Early growth r (true):       %.3f /day, doubling time %.2f days\n",
            r_true, doubling_time(r_true)))
cat(sprintf("Early growth r (naive/7):    %.3f /day, doubling time %.2f days\n",
            r_naive_div, doubling_time(r_naive_div)))
cat(sprintf("Early growth r (as rate):    %.3f /day, doubling time %.2f days  <-- WAY TOO FAST\n\n",
            r_as_rate, doubling_time(r_as_rate)))

cat("=== Fit quality on simulated weekly data (Binomial NLL; lower is better) ===\n")
cat(sprintf("NLL (correct):       %.1f\n", nll_true))
cat(sprintf("NLL (naive /7):      %.1f  (Δ=%.1f)\n", nll_naive_div, nll_naive_div - nll_true))
cat(sprintf("NLL (as daily rate): %.1f  (Δ=%.1f)\n\n", nll_as_rate, nll_as_rate - nll_true))

# Plot: weekly prevalence (points) vs three model curves
op <- par(mfrow=c(1,1), mar=c(4,4,2,1))
plot(dat$time, dat$y/dat$m, pch=16, ylim=c(0,1),
     xlab="Time (days)", ylab="Prevalence (observed & predicted)",
     main="Mixing per-week and per-day rates: impact on dynamics")
lines(times_week, p_true,      lwd=2)           # correct
lines(times_week, p_naive_div, lwd=2, lty=2)    # naive /7
lines(times_week, p_as_rate,   lwd=2, lty=3)    # wrong-as-rate
legend("bottomright",
       legend=c("Observed (weekly)", "Correct mapping", "Naive p/7", "Wrong: weekly prob as daily rate"),
       pch=c(16, NA, NA, NA), lty=c(NA,1,2,3), lwd=c(NA,2,2,2), bty="n")
par(op)

# -------------------------
# 5) Exercises
# -------------------------
cat("\nExercises:\n",
    "1) Change p_week (e.g., 0.3, 0.5, 0.9) and rerun. Is the naive p/7 error bigger when p is large?\n",
    "2) Switch R0 (e.g., 1.5, 3, 5) and see how doubling time & trajectories shift.\n",
    "3) Replace SIS ODE with a weekly discrete update using p = 1 - exp(-lambda * 7) vs p ≈ lambda*7.\n",
    "4) Try a 'mix-up': set beta per day but gamma per week (without converting). What happens to R0 and fit?\n", sep="")




