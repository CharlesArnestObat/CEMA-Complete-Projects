#--------------------------------------------------------------------------------------------------------
# PREAMBLE: Model

# Simple SIS (humans) - SI (mosquitoes) malaria model with vaccine + health economics

# Human SIS split into unvaccinated (Su, Iu) and vaccinated (Sv, Iv); mosquito SI (Svec, Ivec).
# Vaccine reduces infection risk for vaccinated susceptibles (VEi), optional waning.
# Runs baseline (no vaccine) vs vaccination campaign (coverage over 30 days starting day 180).
# Computes and discounts over 5 years:
# Cases, deaths, DALYs (YLD + YLL), QALY losses, costs (treatment + vaccination).
# DALYs averted, QALYs gained, ICER ($/DALY, $/QALY).
# Various plots
#--------------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------------
# Helper functions
#--------------------------------------------------------------------------------------------------------

library(deSolve)
library(tidyverse)

# ----------- 1) ODE system -----------
sis_si_vax <- function(t, state, parms){
  with(as.list(c(state, parms)), {
    
    
    # Su — susceptible, unvaccinated (persons)
    # Iu — infected, unvaccinated (persons)
    # Sv — susceptible, vaccinated (persons)
    # Iv — infected, vaccinated (persons)
    # Nh — total humans = Su + Iu + Sv + Iv
    
    # Svec — susceptible mosquitoes
    # Ivec — infected mosquitoes
    # Nv — total mosquitoes = Svec + Ivec 
    
    Nh <- Su + Iu + Sv + Iv
    Nv <- Svec + Ivec
    xh <- (Iu + Iv) / Nh
    xv <- Ivec / Nv

    # Human FOI from mosquitoes (Ross–Macdonald style)
    # a  = (bites / mosquito / day): mosquito biting rate on humans.
    # b  = (unitless, 0–1): per-bite transmission probability from an infectious mosquito to a human.
    # m  = (mosquitoes / human): mosquito density per person
    # xv = (unitless, 0–1): fraction of mosquitoes that are infectious (sporozoite rate).
    # VEi = vaccine efficacy against infection
    
    lambda_h <- a * b * m * xv
    lambda_v <- (1 - VEi) * lambda_h  # vaccinated susceptibles have reduced FOI

    # Vaccination campaign: instantaneous rate nu(t) on Su only during [t_vacc, t_vacc+dur]
    nu_t <- ifelse(t >= t_vacc & t < (t_vacc + campaign_days),
                   -log(1 - coverage) / campaign_days, 0)

    # Flows
    dSu <- - lambda_h * Su + r * Iu + omega * Sv - nu_t * Su
    dIu <-   lambda_h * Su - r * Iu

    dSv <- - lambda_v * Sv + r * Iv - omega * Sv + nu_t * Su
    dIv <-   lambda_v * Sv - r * Iv

    # Mosquito SI
    c_eff <- c * (p ^ n)
    dIvec <-   a * c_eff * xh * Svec  - mu_v * Ivec
    dSvec <- - a * c_eff * xh * Svec  + mu_v * Ivec

    # Incidence and outcomes (per day)
    inc_inf  <- lambda_h * Su + lambda_v * Sv # incidence of (human) infection today
    cases    <- symptomatic_frac * inc_inf
    deaths   <- CFR * cases
    YLD_flow <- cases * (DW * dur_ill_days / 365)   # years
    YLL_flow <- deaths * LYL_years                  # years
    QALY_loss_flow <- cases * (util_loss * dur_ill_days / 365) + deaths * LYL_years

    # Vaccinations per day
    vacc_today <- nu_t * Su

    # Costs (per day)
    cost_flow <- (cases * cost_per_case) + (vacc_today * cost_per_vacc)

    # Discount factor (continuous-time approx)
    disc <- 1 / ((1 + disc_rate) ^ (t / 365))

    # Cumulatives (undiscounted)
    dC_cases <- cases
    dC_deaths <- deaths
    dC_YLD <- YLD_flow
    dC_YLL <- YLL_flow
    dC_DALY <- YLD_flow + YLL_flow
    dC_QALYloss <- QALY_loss_flow
    dC_vacc <- vacc_today
    dC_cost <- cost_flow

    # Cumulatives (discounted)
    dC_cost_disc <- cost_flow * disc
    dC_DALY_disc <- (YLD_flow + YLL_flow) * disc
    dC_QALYloss_disc <- QALY_loss_flow * disc

    list(c(dSu, dIu, dSv, dIv, dSvec, dIvec,
           dC_cases, dC_deaths, dC_YLD, dC_YLL, dC_DALY, dC_QALYloss, dC_vacc, dC_cost,
           dC_cost_disc, dC_DALY_disc, dC_QALYloss_disc),
         c(xh = xh, xv = xv, lambda_h = lambda_h, vacc_today = vacc_today, disc = disc))
  })
}

# ----------- 2) Simulation wrappers -----------------
run_scenario <- function(params, times, init){
  out <- ode(y = init, times = times, func = sis_si_vax, parms = params, method = "lsoda")
  as.data.frame(out)
}
summarise_results <- function(df){
  tail_row <- df[nrow(df), ]
  with(tail_row, {
    list(
      prevalence_end = xh,
      cases = C_cases,
      deaths = C_deaths,
      DALY = C_DALY,
      QALY_loss = C_QALYloss,
      vaccinations = C_vacc,
      cost = C_cost,
      cost_disc = C_cost_disc,
      DALY_disc = C_DALY_disc,
      QALYloss_disc = C_QALYloss_disc
    )
  })
} # this helper pulls end-of-simulation totals from the ODE output and returns them as a tidy list

# ----------- 3) Analysis (baseline vs vaccine) -----------
run_analysis <- function(horizon_years = 5){
  times <- seq(0, horizon_years*365, by=1)

  # Baseline
  base_params <- make_params(vax = FALSE)
  base_init   <- make_init(Nh=100000, m=base_params$m, R0=(base_params$a*base_params$b*base_params$m)/base_params$r)
  base_df     <- run_scenario(base_params, times, base_init)
  base_sum    <- summarise_results(base_df)

  # Vaccination
  vax_params <- make_params(vax = TRUE)
  vax_init   <- make_init(Nh=100000, m=vax_params$m, R0=(vax_params$a*vax_params$b*vax_params$m)/vax_params$r)
  vax_df     <- run_scenario(vax_params, times, vax_init)
  vax_sum    <- summarise_results(vax_df)

  # Differences (discounted)
  DALY_averted_disc  <- base_sum$DALY_disc  - vax_sum$DALY_disc
  QALY_gained_disc   <- base_sum$QALYloss_disc - vax_sum$QALYloss_disc
  dCost_disc         <- vax_sum$cost_disc - base_sum$cost_disc

  ICER_DALY <- dCost_disc / DALY_averted_disc
  ICER_QALY <- dCost_disc / QALY_gained_disc

  list(base=base_sum, vax=vax_sum,
       DALY_averted_disc=DALY_averted_disc,
       QALY_gained_disc=QALY_gained_disc,
       dCost_disc=dCost_disc,
       ICER_DALY=ICER_DALY,
       ICER_QALY=ICER_QALY,
       base_df=base_df, vax_df=vax_df)
}


# ----------- 4) CE plane with PSA cloud (many points) -----------
plot_ce_plane_psa <- function(N = 300, effect = c("QALY","DALY"), k = 1000, seed = 1,
                              alpha = 0.5, save_file = FALSE, filename = "ce_plane_psa.png"){
  effect <- match.arg(effect)
  # Reuse the PSA machinery from build_ceac (does not save files here)
  out <- build_ceac(N = N, effect = effect, wtp_max = 1000, wtp_step = 50, seed = seed, save_files = FALSE)
  psa <- out$psa
  labx <- if (effect == "QALY") "Incremental effect (QALYs gained)" else "Incremental effect (DALYs averted)"
  p <- ggplot(psa, aes(dE, dC)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_abline(intercept = 0, slope = k, linetype = "dashed") +
    geom_point(alpha = alpha, size = 1.2, color = "red") +
    stat_ellipse(type = "norm", level = 0.50, linetype = "solid") +
    stat_ellipse(type = "norm", level = 0.95, linetype = "dotted") +
    labs(x = labx, y = "Incremental cost ($)",
         title = "Cost-effectiveness plane (PSA cloud)",
         subtitle = paste("PSA draws N =", N, "| WTP =", k, if (effect=="QALY") "per QALY" else "per DALY")) +
    theme_minimal(base_size = 12)
  if (save_file) ggsave(filename, p, width = 6.5, height = 4.5, dpi = 200)
  p
}

# ----------- 5) CEAC with PSA -----------

# Helper to run with parameter overrides
run_analysis_override <- function(horizon_years = 5, override = list()){
  times <- seq(0, horizon_years*365, by=1)
  
  base_params <- make_params(vax = FALSE)
  for(nm in intersect(names(override), names(base_params))) base_params[[nm]] <- override[[nm]]
  base_init <- make_init(Nh=100000, m=base_params$m,
                         R0=(base_params$a*base_params$b*base_params$m)/base_params$r)
  base_df  <- run_scenario(base_params, times, base_init)
  base_sum <- summarise_results(base_df)
  
  vax_params <- make_params(vax = TRUE)
  for(nm in intersect(names(override), names(vax_params))) vax_params[[nm]] <- override[[nm]]
  vax_init <- make_init(Nh=100000, m=vax_params$m,
                        R0=(vax_params$a*vax_params$b*vax_params$m)/vax_params$r)
  vax_df  <- run_scenario(vax_params, times, vax_init)
  vax_sum <- summarise_results(vax_df)
  
  list(
    dC = vax_sum$cost_disc - base_sum$cost_disc,
    dQ = base_sum$QALYloss_disc - vax_sum$QALYloss_disc,
    dD = base_sum$DALY_disc - vax_sum$DALY_disc
  )
}

# Simple sampler for uncertain inputs (transparent)
# Change values for the PSA variables here
psa_draw <- function(){
  list(
    VEi = rbeta(1, 20, 80),
    coverage = rbeta(1, 60, 40),
    cost_per_vacc = rlnorm(1, log(10), 0.25),
    cost_per_case = rlnorm(1, log(10), 0.30),
    symptomatic_frac = rbeta(1, 60, 40),
    CFR = rbeta(1, 2, 1998)
  )
}

build_ceac <- function(N = 300, effect = c("QALY","DALY"),
                       horizon_years = 5,
                       wtp_max = 2000, wtp_step = 50, seed = 1,
                       save_files = FALSE){
  effect <- match.arg(effect)
  set.seed(seed)
  
  psa <- replicate(N, {
    o <- psa_draw()
    out <- run_analysis_override(horizon_years, override = o)
    c(dC = out$dC,
      dE = if (effect == "QALY") out$dQ else out$dD)
  })
  psa <- as.data.frame(t(psa))
  
  k_grid <- seq(0, wtp_max, by = wtp_step)
  ce <- sapply(k_grid, function(k) mean((k * psa$dE - psa$dC) > 0))
  se <- sqrt(ce * (1 - ce) / N)
  ceac <- data.frame(WTP = k_grid, CE_Prob = ce, SE = se,
                     CE_L = pmax(0, ce - 1.96 * se),
                     CE_U = pmin(1, ce + 1.96 * se))
  
  xlab <- if (effect == "QALY") "Willingness to pay ($ per QALY gained)"
  else "Willingness to pay ($ per DALY averted)"
  
  p <- ggplot(ceac, aes(WTP, CE_Prob)) +
    geom_ribbon(aes(ymin = CE_L, ymax = CE_U),fill="red", alpha = 0.15) +
    geom_line(linewidth = 1, color = "red") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.2)) +
    labs(x = xlab, y = "Probability cost-effective",
         title = "Cost-Effectiveness Acceptability Curve (CEAC)",
         subtitle = paste0("PSA draws: N = ", N, ", horizon = ", horizon_years, " years; effect = ", effect)) +
    theme_minimal(base_size = 12)
  
  if (save_files){
    ggsave("ceac_plot.png", p, width = 7, height = 4.5, dpi = 200)
    write.csv(ceac, "ceac_table.csv", row.names = FALSE)
  }
  
  list(psa = psa, ceac = ceac, plot = p)
}


# ----------- 6) Defaults ------ -----------

# r = 1/200 ≈ 0.005 (/day) — human recovery rate (mean infectious duration = 200 days).
# a = 0.25 (bites/mosquito/day) — mosquito biting rate on humans.
# b = 0.05 (unitless, 0–1) — per-bite mosquito→human infection probability.
# c = 0.50 (unitless, 0–1) — per-bite human→mosquito infection probability.
# p = 0.90 (unitless, 0–1) — daily mosquito survival probability.
# n = 10 (days) — EIP: extrinsic incubation period in mosquitoes.
#     Often used via ceff=c*p^n (fraction surviving EIP).
# m = 2.0 (mosquitoes/person) — mosquito density per human.
# μ_v (mu_v) = 0.10 (/day) — mosquito death rate (mean life ≈ 10 days).

# VEi = 0.60 (unitless, 0–1 when vax=TRUE, else 0) — vaccine efficacy against infection (reduces FOI for vaccinated susceptibles).
# ω (omega) = 1/(3·365) (/day; when vax=TRUE, else 0) — waning rate from vaccinated back to unvaccinated; mean protection ≈ 3 years.
#   Half-life ≈ ln2/ω≈ 2.1 years.

# t_vacc = 180 (days) — campaign start day.
# campaign_days = 30 (days if vax=TRUE, else 0) — campaign length.
# coverage = 0.60 (fraction if vax=TRUE, else 0) — target fraction of currently unvaccinated susceptibles to vaccinate over the campaign.
# Internally converted to a constant daily hazard: 
#   ν=−ln⁡(1−coverage)/campaign_days

# symptomatic_frac = 0.60 (unitless) — fraction of infections that become symptomatic cases.
# dur_ill_days = 5 (days) — average symptomatic duration per case.
# DW = 0.21 (unitless, 0–1) — disability weight (DALY) for illness episode.
# CFR = 0.001 (unitless) — case-fatality risk among symptomatic cases.
# LYL_years = 30 (years) — life-years lost per death (teaching simplification).
# util_loss = 0.20 (unitless, 0–1) — temporary utility decrement during illness (for QALY losses).

# cost_per_case = 10 (currency per case) — treatment/health-care cost per symptomatic case.
# cost_per_vacc = 10 (currency per person) — vaccine + delivery cost per person vaccinated.
# disc_rate = 0.03 (per year) — annual discount rate applied to costs and health outcomes.

make_params <- function(vax = TRUE){
  list(
    a = 0.25,  b = 0.05, c = 0.50, p = 0.90, n = 10,
    m = 2.0,   mu_v = 0.10,
    r = 1/200,
    VEi = if (vax) 0.20 else 0.0,
    omega = if (vax) 1/(3*365) else 0.0,
    t_vacc = 180, campaign_days = if (vax) 30 else 0, coverage = if (vax) 0.60 else 0.0,
    symptomatic_frac = 0.60, dur_ill_days = 5, DW = 0.21,
    CFR = 0.001, LYL_years = 30, util_loss = 0.20,
    cost_per_case = 10, cost_per_vacc = 10,
    disc_rate = 0.03
  )
}

make_init <- function(Nh = 100000, m = 2.0, R0 = 3.0){
  xstar <- 1 - 1/R0
  I0 <- round(xstar * Nh)
  Su <- Nh - I0
  Iu <- I0
  Sv <- 0
  Iv <- 0
  Nv <- Nh * m
  Ivec0 <- round(0.1 * Nv)
  Svec0 <- Nv - Ivec0
  
  c(Su=Su, Iu=Iu, Sv=Sv, Iv=Iv, Svec=Svec0, Ivec=Ivec0,
    C_cases=0, C_deaths=0, C_YLD=0, C_YLL=0, C_DALY=0, C_QALYloss=0, C_vacc=0, C_cost=0,
    C_cost_disc=0, C_DALY_disc=0, C_QALYloss_disc=0)
}


#--------------------------------------------------------------------------------------------------------
# Outputs
#--------------------------------------------------------------------------------------------------------

# ----------- 5) Main: run & print table + plots -----------
res <- run_analysis(horizon_years = 10)

# Summary TABLE (nicely labeled)
summary_table <- data.frame(
  Metric = c("Discounted cost ($)", "Discounted DALYs", "Discounted QALY loss",
             "Vaccinations (count)", "DALYs averted", "QALYs gained", "ΔCost ($)",
             "ICER ($/DALY)", "ICER ($/QALY)"),
  Baseline = c(res$base$cost_disc, res$base$DALY_disc, res$base$QALYloss_disc,
               res$base$vaccinations, NA, NA, NA, NA, NA),
  Vaccine  = c(res$vax$cost_disc, res$vax$DALY_disc, res$vax$QALYloss_disc,
               res$vax$vaccinations, NA, NA, NA, NA, NA),
  Difference = c(res$dCost_disc, res$base$DALY_disc - res$vax$DALY_disc,
                 res$base$QALYloss_disc - res$vax$QALYloss_disc,
                 res$vax$vaccinations - res$base$vaccinations,
                 res$DALY_averted_disc, res$QALY_gained_disc, res$dCost_disc,
                 res$ICER_DALY, res$ICER_QALY)
)

print(summary_table)

## Prevalence (baseline vs vaccine)
df_prev <- bind_rows(
  transmute(res$base_df, time_year = time/365, value = xh, scenario = "Baseline"),
  transmute(res$vax_df,  time_year = time/365, value = xh, scenario = "Vaccine")
)

p1 <- ggplot(df_prev, aes(time_year, value, colour = scenario, linetype = scenario)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 180/365, linetype = "dotted") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Time (years)", y = "Prevalence (I/N)",
       title = "Prevalence (baseline vs vaccine)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

p1

## Cumulative discounted cost
df_cost <- bind_rows(
  transmute(res$base_df, time_year = time/365, value = C_cost_disc, scenario = "Baseline"),
  transmute(res$vax_df,  time_year = time/365, value = C_cost_disc, scenario = "Vaccine")
)

p2 <- ggplot(df_cost, aes(time_year, value, colour = scenario, linetype = scenario)) +
  geom_line(linewidth = 1) +
  labs(x = "Time (years)", y = "Cum. discounted cost ($)",
       title = "Cumulative discounted cost") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

p2

## Cumulative discounted DALYs
df_daly <- bind_rows(
  transmute(res$base_df, time_year = time/365, value = C_DALY_disc, scenario = "Baseline"),
  transmute(res$vax_df,  time_year = time/365, value = C_DALY_disc, scenario = "Vaccine")
)

p3 <- ggplot(df_daly, aes(time_year, value, colour = scenario, linetype = scenario)) +
  geom_line(linewidth = 1) +
  labs(x = "Time (years)", y = "Cum. discounted DALYs",
       title = "Cumulative discounted DALYs") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

p3

## (Optional) arrange in one row:
# install.packages("patchwork")  # once
library(patchwork)
p1 + p2 + p3 + plot_layout(nrow = 1)


# ----------- 5) Main: run & print table + plots -----------

# CEAC (PSA)
out_ceac <- build_ceac(N = 200, effect = "QALY", wtp_max = 2000, wtp_step = 50, seed = 1, save_files = TRUE)
print(out_ceac$plot)

# PSA CE plane (multi-point)
p_ce_psa <- plot_ce_plane_psa(N = 200, effect = "QALY", k = 1000, seed = 1, filename = "ce_plane_psa.png")
print(p_ce_psa)


# Stuff for you to do 
#--------------------------

cat("\nExercises:\n",
    "1) Change VEi and coverage; try value 0.2 and 0.9 observe table & plots.\n",
    "2) cost_per_vacc; try values 5, 30, 40 \n",
    "3) Extend horizon to 10 years \n",
    "4) Increase CFR or symptomatic_frac and see DALY/cost impact.\n",
    "5) Try a vaccine that also reduces infectiousness (modify FOI for Iv).\n", sep="")


