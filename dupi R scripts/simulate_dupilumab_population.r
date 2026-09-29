## ---------------------------------------------------------------------------
## simulate_dupilumab_population.R
##
## Simulate a virtual population of asthma patients dosed with dupilumab
## 300 mg q2w (600 mg loading dose) using the mrgsolve PopPK model
## (dupilumab_asthma_poppk.cpp), including inter-individual variability (IIV)
## and covariate distributions representative of the asthma PopPK dataset
## (Zhang et al. 2021, CPT Pharmacometrics Syst Pharmacol).
##
## Steady-state exposure metrics (AUCss, Cmax,ss, Cmin,ss) are computed from
## the last dosing interval and compared against the model-derived exposures
## reported in Table 4 of the manuscript for study NCT02414854 (300 mg q2w,
## non-OCS-dependent adult/adolescent asthma patients):
##
##   AUCss   (mg*day/L): 1090 (593)  [54.4% CV]
##   Cmax,ss (mg/L)     : 86.9 (44.8) [51.6% CV]
##   Cmin,ss (mg/L)     : 70.0 (40.9) [58.5% CV]
## ---------------------------------------------------------------------------

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(4783)

## -----------------------------------------------------------------------
## 1. Load model
## -----------------------------------------------------------------------
mod <- mread("dupilumab_asthma_poppk.cpp")

## -----------------------------------------------------------------------
## 2. Define virtual population covariates
##
## Distributions approximate the pooled asthma PopPK dataset (Table 2):
##   Weight (kg): median 78, range 32-186   -> lognormal, truncated
##   Albumin (g/L): median 44, range 30-55  -> normal, truncated
##   CLCRN (mL/min/1.73m2): median 111, range 30-377 -> lognormal, truncated
##   ADA status: ~14.5% non-negative (positive) at baseline/stationary
## -----------------------------------------------------------------------
n_id <- 1000

# Weight: lognormal with median 78 kg, moderate spread, truncated to observed range
wt_cv <- 0.28
weight <- rlnorm(n_id, meanlog = log(78), sdlog = sqrt(log(1 + wt_cv^2)))
weight <- pmin(pmax(weight, 32), 186)

# Albumin: normal with mean/median ~44 g/L, SD ~3.5, truncated to observed range
albumin <- rnorm(n_id, mean = 44, sd = 3.5)
albumin <- pmin(pmax(albumin, 30), 55)

# CLCRN: lognormal with median 111 mL/min/1.73m2, truncated to observed range
clcrn_cv <- 0.30
clcrn <- rlnorm(n_id, meanlog = log(111), sdlog = sqrt(log(1 + clcrn_cv^2)))
clcrn <- pmin(pmax(clcrn, 30), 377)

# ADA status: ~14.5% non-negative (positive) per Table 2 pooled proportion
ada <- rbinom(n_id, size = 1, prob = 0.145)

covset <- tibble(
  ID     = seq_len(n_id),
  WT     = weight,
  ALB    = albumin,
  CLCRN  = clcrn,
  ADA    = ada
)

## -----------------------------------------------------------------------
## 3. Build dosing regimen: 300 mg q2w with 600 mg loading dose
##
## NOTE: loading dose uses addl = 0 (single dose at time 0); maintenance
## dosing begins explicitly at day 14 to avoid duplicate/overlapping doses.
## -----------------------------------------------------------------------
n_doses_maintenance <- 25   # ~52 weeks of q2w maintenance dosing after loading
end_time <- 14 * (n_doses_maintenance + 1) + 14  # simulate one interval past last dose

dose_regimen <- ev(amt = 600, cmt = "GUT", time = 0) |>
  seq(ev(amt = 300, cmt = "GUT", ii = 14, addl = n_doses_maintenance - 1, time = 14))

data_in <- as.data.frame(dose_regimen) |>
  tidyr::crossing(ID = covset$ID) |>
  left_join(covset, by = "ID") |>
  select(ID, time, amt, ii, addl, cmt, evid, WT, ALB, CLCRN, ADA) |>
  arrange(ID, time)

## -----------------------------------------------------------------------
## 4. Simulate with IIV (omega already specified in model)
## -----------------------------------------------------------------------
sim <- mod |>
  data_set(data_in) |>
  Req(CP) |>
  mrgsim(end = end_time, delta = 0.5) |>
  as_tibble()

## -----------------------------------------------------------------------
## 5. Compute steady-state exposure metrics from the last dosing interval
##    (last maintenance dose at time = 14 + 14*(n_doses_maintenance-1))
## -----------------------------------------------------------------------
last_dose_time <- 14 + 14 * (n_doses_maintenance - 1)
interval_start <- last_dose_time
interval_end   <- last_dose_time + 14

ss_window <- sim |>
  filter(time >= interval_start, time <= interval_end)

exposure_metrics <- ss_window |>
  group_by(ID) |>
  summarise(
    Cmax_ss = max(CP),
    Cmin_ss = min(CP),
    AUC_ss  = sum(diff(time) * (head(CP, -1) + tail(CP, -1)) / 2),  # trapezoidal
    .groups = "drop"
  )

## -----------------------------------------------------------------------
## 6. Summary statistics (mean, SD, %CV) and comparison to Table 4
## -----------------------------------------------------------------------
summarise_metric <- function(x) {
  m <- mean(x)
  s <- sd(x)
  cv <- 100 * s / m
  c(mean = m, sd = s, cv_pct = cv)
}

sim_summary <- exposure_metrics |>
  summarise(
    AUCss_mean = mean(AUC_ss),  AUCss_sd  = sd(AUC_ss),  AUCss_cv  = 100 * sd(AUC_ss)  / mean(AUC_ss),
    Cmax_mean  = mean(Cmax_ss), Cmax_sd   = sd(Cmax_ss), Cmax_cv   = 100 * sd(Cmax_ss) / mean(Cmax_ss),
    Cmin_mean  = mean(Cmin_ss), Cmin_sd   = sd(Cmin_ss), Cmin_cv   = 100 * sd(Cmin_ss) / mean(Cmin_ss)
  )

table4_reference <- tibble(
  metric        = c("AUCss (mg*day/L)", "Cmax,ss (mg/L)", "Cmin,ss (mg/L)"),
  reported_mean = c(1090, 86.9, 70.0),
  reported_sd   = c(593, 44.8, 40.9),
  reported_cv   = c(54.4, 51.6, 58.5)
)

comparison <- tibble(
  metric        = c("AUCss (mg*day/L)", "Cmax,ss (mg/L)", "Cmin,ss (mg/L)"),
  sim_mean      = c(sim_summary$AUCss_mean, sim_summary$Cmax_mean, sim_summary$Cmin_mean),
  sim_sd        = c(sim_summary$AUCss_sd, sim_summary$Cmax_sd, sim_summary$Cmin_sd),
  sim_cv        = c(sim_summary$AUCss_cv, sim_summary$Cmax_cv, sim_summary$Cmin_cv)
) |>
  left_join(table4_reference, by = "metric")

print(comparison, width = Inf)

## -----------------------------------------------------------------------
## 7. Visualize simulated exposure distributions vs. reported means
## -----------------------------------------------------------------------
exposure_long <- exposure_metrics |>
  pivot_longer(cols = c(AUC_ss, Cmax_ss, Cmin_ss), names_to = "metric", values_to = "value") |>
  mutate(metric = recode(metric,
                          AUC_ss  = "AUCss (mg*day/L)",
                          Cmax_ss = "Cmax,ss (mg/L)",
                          Cmin_ss = "Cmin,ss (mg/L)"))

ref_long <- table4_reference |>
  rename(value = reported_mean) |>
  select(metric, value)

ggplot(exposure_long, aes(x = value)) +
  geom_histogram(bins = 40) +
  geom_vline(data = ref_long, aes(xintercept = value), linetype = "dashed") +
  facet_wrap(~ metric, scales = "free") +
  labs(
    x = "Exposure",
    y = "Count",
    title = "Simulated steady-state exposures (300 mg q2w) vs. reported mean (dashed line)",
    subtitle = "Reference: Table 4, Zhang et al. 2021 (NCT02414854, non-OCS asthma patients)"
  )
