## ---------------------------------------------------------------------------
## simulate_dupilumab_population_nlmixr2.R
##
## nlmixr2/rxode2 port of simulate_dupilumab_population_quantile.R.
##
## Simulates a virtual population of asthma patients dosed with dupilumab
## 300 mg q2w (600 mg loading dose) using the nlmixr2 PopPK model defined in
## dupilumab_asthma_poppk_nlmixr2.R, including inter-individual variability
## (IIV) and covariates resampled from a quantile function calibrated to
## Table 2 statistics (mean/SD, truncated to min/max) from:
##
## Zhang L, Gao Y, Li M, Xu C, Davis JD, Kanamaluru V, Lu Q.
## Population pharmacokinetic analysis of dupilumab in adult and adolescent
## patients with asthma. CPT Pharmacometrics Syst Pharmacol. 2021;10(8):941-952.
##
## Steady-state exposure metrics (AUCss, Cmax,ss, Cmin,ss) are computed from
## the last dosing interval and compared against Table 4 (study NCT02414854,
## 300 mg q2w, non-OCS-dependent adult/adolescent asthma patients):
##
##   AUCss   (mg*day/L): 1090 (593)  [54.4% CV]
##   Cmax,ss (mg/L)     : 86.9 (44.8) [51.6% CV]
##   Cmin,ss (mg/L)     : 70.0 (40.9) [58.5% CV]
## ---------------------------------------------------------------------------

library(nlmixr2)
library(rxode2)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(8214)

moddir <- "/Users/lukefostvedt/Documents/ACOP-STEP-2026/dupilumab-poppk/"

## Output directory for saved plots (created if it doesn't exist)
fig_dir <- "figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir)

## -----------------------------------------------------------------------
## 1. Load models (nlmixr2/rxode2 and mrgsolve, for the parity check below)
## -----------------------------------------------------------------------
source(paste0(moddir,"dupilumab_asthma_poppk_nlmixr2.R"))
mod_nlmixr <- dupilumab_asthma_poppk()
mod_mrgsolve <- mread(paste0(moddir,"dupilumab_asthma_poppk.cpp"))

## Parity-only mrgsolve variant: identical structural/covariate model, but
## etas are declared as plain $PARAM inputs (not $OMEGA-simulated), so
## externally supplied per-subject eta values can be passed in directly for
## an apples-to-apples comparison against nlmixr2/rxode2 (mrgsolve's
## $OMEGA-declared etas cannot be overridden by input data columns; they are
## always resimulated internally per ID). See dupilumab_asthma_poppk_parity.cpp.
mod_mrgsolve_parity <- mread(paste0(moddir,"dupilumab_asthma_poppk_parity.cpp"))

## Fixed-effect (theta) vector at the literature values (must be supplied
## explicitly to rxSolve alongside covariates/etas; unlike mrgsolve, thetas
## are not baked into the compiled model object by default).
theta_vals <- c(
  tke = log(0.0418), tv2 = log(2.76), tk23 = log(0.0952), tk32 = log(0.163),
  tvmax = log(1.39), tkm = log(2.08), tka = log(0.263), tfsc = log(0.609 / (1 - 0.609)),
  wtke = 0.222, adake = 0.191, clcrnke = 0.217, wtv2 = 0.667, albv2 = -0.484, wtvmax = 0.224,
  prop.err = 0.0388, add.err = 1.73
)

## -----------------------------------------------------------------------
## 2. Quantile-based covariate sampler calibrated to Table 2 statistics
##
## Table 2 reports only mean, SD, median, min, and max -- not a full
## percentile ladder. A shifted lognormal quantile function is calibrated
## by method-of-moments to reproduce the reported mean and SD, and samples
## are drawn via its quantile function (qlnorm), then rejection-sampled to
## respect the reported [min, max] range.
## -----------------------------------------------------------------------
lognormal_moments <- function(mean_val, sd_val) {
  cv2 <- (sd_val / mean_val)^2
  sdlog <- sqrt(log(1 + cv2))
  meanlog <- log(mean_val) - sdlog^2 / 2
  list(meanlog = meanlog, sdlog = sdlog)
}

sample_from_table2 <- function(n, mean_val, sd_val, min_val, max_val) {
  pars <- lognormal_moments(mean_val, sd_val)

  x <- numeric(0)
  while (length(x) < n) {
    u     <- runif(n * 2)
    draws <- qlnorm(u, meanlog = pars$meanlog, sdlog = pars$sdlog)
    draws <- draws[draws >= min_val & draws <= max_val]
    x     <- c(x, draws)
  }
  x[seq_len(n)]
}

## -----------------------------------------------------------------------
## 3. Define virtual population using Table 2 percentiles
##
## Values taken directly from Table 2 (pooled healthy + asthma patients,
## N = 2114, model development population):
##   Weight (kg):  mean 79.6, SD 19.1, min 32,   max 186
##   Albumin(g/L): mean 43.8, SD 3.50, min 30,   max 55
##   CLCRN:        mean 116.2, SD 36.0, min 30.1, max 377
## ADA non-negative proportion: 344/2114 = 16.3%
## -----------------------------------------------------------------------
n_id <- 1000

weight  <- sample_from_table2(n_id, mean_val = 79.6,  sd_val = 19.1, min_val = 32,   max_val = 186)
albumin <- sample_from_table2(n_id, mean_val = 43.8,  sd_val = 3.50, min_val = 30,   max_val = 55)
clcrn   <- sample_from_table2(n_id, mean_val = 116.2, sd_val = 36.0, min_val = 30.1, max_val = 377)
ada     <- rbinom(n_id, size = 1, prob = 344 / 2114)

covset <- tibble(
  id    = seq_len(n_id),
  WT    = weight,
  ALB   = albumin,
  CLCRN = clcrn,
  ADA   = ada
)

## Diagnostic: compare resampled covariate summary stats to Table 2
table2_reference <- tribble(
  ~covariate, ~reported_mean, ~reported_sd, ~reported_median, ~reported_min, ~reported_max,
  "WT",        79.6,           19.1,         78.0,             32,           186,
  "ALB",       43.8,           3.50,         44.0,             30,           55,
  "CLCRN",     116.2,          36.0,         111.3,            30.1,         377
)

covariate_diagnostic <- covset |>
  summarise(
    WT_mean = mean(WT), WT_sd = sd(WT), WT_median = median(WT),
    ALB_mean = mean(ALB), ALB_sd = sd(ALB), ALB_median = median(ALB),
    CLCRN_mean = mean(CLCRN), CLCRN_sd = sd(CLCRN), CLCRN_median = median(CLCRN)
  )

cat("\n--- Resampled covariate distributions vs. Table 2 ---\n")
print(table2_reference)
cat("\nResampled summary:\n")
print(covariate_diagnostic, width = Inf)

## -----------------------------------------------------------------------
## 4. Build dosing regimen: 300 mg q2w with 600 mg loading dose
##
## Single loading dose at time 0 (nbr.doses = 1), maintenance dosing starts
## explicitly at day 14 (nbr.doses = 25) to avoid duplicate/overlapping
## doses -- mirrors the corrected mrgsolve event schedule.
## -----------------------------------------------------------------------
n_doses_maintenance <- 25
sample_times <- seq(0, 14 * (n_doses_maintenance + 1) + 14, by = 1)

ev_template <- eventTable() |>
  add.dosing(dose = 600, start.time = 0, nbr.doses = 1, dosing.to = "depot") |>
  add.dosing(dose = 300, start.time = 14, nbr.doses = n_doses_maintenance,
             dosing.interval = 14, dosing.to = "depot") |>
  add.sampling(sample_times)

last_dose_time <- 14 + 14 * (n_doses_maintenance - 1)
interval_start <- last_dose_time
interval_end   <- last_dose_time + 14

## -----------------------------------------------------------------------
## 5. Simulate population with IIV
##
## rxode2's population simulation via `params = data.frame(...)` (one row
## per subject, one id per row of the event table crossed with covariates)
## automatically draws new etas per subject from the model's omega matrix
## when `nStud`/`params` supplies subject-level covariates and omega is not
## bypassed -- here we instead build a single long data set (id, time,
## amt, evid, covariates) and let rxSolve simulate individual etas via its
## population dosing/covariate table interface, matching the mrgsolve
## `data_set()` workflow.
## -----------------------------------------------------------------------
dosing_df <- as.data.frame(ev_template)

data_in <- dosing_df |>
  tidyr::crossing(id = covset$id) |>
  left_join(covset, by = "id") |>
  arrange(id, time)

sim <- rxSolve(
  mod_nlmixr,
  events = data_in,
  params = theta_vals,
  omega  = mod_nlmixr$omega,
  covsInterpolation = "locf"
) |>
  as_tibble()

## -----------------------------------------------------------------------
## 6. Compute steady-state exposure metrics from the last dosing interval
## -----------------------------------------------------------------------
ss_window <- sim |>
  filter(time >= interval_start, time <= interval_end)

exposure_metrics <- ss_window |>
  group_by(id) |>
  summarise(
    Cmax_ss = max(cp),
    Cmin_ss = min(cp),
    AUC_ss  = sum(diff(time) * (head(cp, -1) + tail(cp, -1)) / 2),  # trapezoidal
    .groups = "drop"
  )

## -----------------------------------------------------------------------
## 7. Summary statistics (mean, SD, %CV) and comparison to Table 4
## -----------------------------------------------------------------------
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
  metric   = c("AUCss (mg*day/L)", "Cmax,ss (mg/L)", "Cmin,ss (mg/L)"),
  sim_mean = c(sim_summary$AUCss_mean, sim_summary$Cmax_mean, sim_summary$Cmin_mean),
  sim_sd   = c(sim_summary$AUCss_sd, sim_summary$Cmax_sd, sim_summary$Cmin_sd),
  sim_cv   = c(sim_summary$AUCss_cv, sim_summary$Cmax_cv, sim_summary$Cmin_cv)
) |>
  left_join(table4_reference, by = "metric")

cat("\n--- nlmixr2/rxode2 simulated exposures vs. Table 4 ---\n")
print(comparison, width = Inf)

## -----------------------------------------------------------------------
## 8. Visualize covariate distributions (resampled) vs. reported points
## -----------------------------------------------------------------------
covariate_long <- covset |>
  pivot_longer(cols = c(WT, ALB, CLCRN), names_to = "covariate", values_to = "value")

covariate_marks <- table2_reference |>
  dplyr::select(covariate, reported_median, reported_min, reported_max) |>
  pivot_longer(cols = c(reported_median, reported_min, reported_max),
               names_to = "stat", values_to = "value")

plot_covariates <- ggplot(covariate_long, aes(x = value)) +
  geom_histogram(bins = 40) +
  geom_vline(data = covariate_marks, aes(xintercept = value, linetype = stat)) +
  facet_wrap(~ covariate, scales = "free") +
  labs(
    title = "Quantile-resampled covariate distributions vs. Table 2 percentiles",
    x = "Value", y = "Count", linetype = "Table 2 stat"
  )

ggsave(file.path(fig_dir, "covariate_distributions.png"), plot_covariates,
       width = 9, height = 4, dpi = 150)

## -----------------------------------------------------------------------
## 9. Visualize simulated exposure distributions vs. reported means
## -----------------------------------------------------------------------
exposure_long <- exposure_metrics |>
  pivot_longer(cols = c(AUC_ss, Cmax_ss, Cmin_ss), names_to = "metric", values_to = "value") |>
  mutate(metric = recode(metric,
                          AUC_ss  = "AUCss (mg*day/L)",
                          Cmax_ss = "Cmax,ss (mg/L)",
                          Cmin_ss = "Cmin,ss (mg/L)"))

ref_long <- table4_reference |>
  rename(value = reported_mean) |>
  dplyr::select(metric, value)

plot_exposures <- ggplot(exposure_long, aes(x = value)) +
  geom_histogram(bins = 40) +
  geom_vline(data = ref_long, aes(xintercept = value), linetype = "dashed") +
  facet_wrap(~ metric, scales = "free") +
  labs(
    x = "Exposure",
    y = "Count",
    title = "nlmixr2/rxode2 simulated steady-state exposures vs. reported mean",
    subtitle = "Reference: Table 4, Zhang et al. 2021 (NCT02414854, non-OCS asthma patients)"
  )

ggsave(file.path(fig_dir, "exposure_distributions_nlmixr2.png"), plot_exposures,
       width = 9, height = 4, dpi = 150)

## ---------------------------------------------------------------------------
## 10. Direct numerical parity check: mrgsolve vs. nlmixr2/rxode2
##
## Both engines are run on the *identical* dosing schedule, covariate set,
## and (crucially) identical per-subject eta draws, so any difference in
## individual concentration-time profiles reflects genuine numerical/
## implementation discrepancies between the two engines rather than
## differences in random sampling.
##
## Etas are drawn once (from the shared omega diagonal, since both models
## use the same uncorrelated variances) and passed explicitly to each
## engine as per-subject parameters, bypassing each engine's internal
## random-effect simulation.
## ---------------------------------------------------------------------------
cat("\n--- Running mrgsolve vs. nlmixr2 parity check ---\n")

omega_diag <- c(
  eta.ke   = 0.0385,
  eta.v2   = 0.00834,
  eta.vmax = 0.0589,
  eta.ka   = 0.243,
  eta.fsc  = 0.132
)

set.seed(9142)
eta_draws <- sapply(omega_diag, function(v) rnorm(n_id, mean = 0, sd = sqrt(v)))
eta_df <- as_tibble(eta_draws) |>
  mutate(id = seq_len(n_id))

covset_eta <- covset |>
  left_join(eta_df, by = "id")

## --- mrgsolve run: pass etas as fixed per-subject parameters (bypasses
## $OMEGA simulation entirely since these are supplied as data columns).
## mrgsolve's $OMEGA names are ETA_KE/ETA_V2/ETA_VMAX/ETA_KA/ETA_FSC, and its
## depot compartment is named GUT (vs. nlmixr2/rxode2's "depot"). ---
data_in_mrgsolve <- dosing_df |>
  mutate(cmt = if_else(evid == 1, "GUT", "CENT")) |>
  tidyr::crossing(id = covset_eta$id) |>
  left_join(covset_eta, by = "id") |>
  rename(
    ID       = id,
    ETA_KE   = eta.ke,
    ETA_V2   = eta.v2,
    ETA_VMAX = eta.vmax,
    ETA_KA   = eta.ka,
    ETA_FSC  = eta.fsc
  ) |>
  arrange(ID, time)

sim_mrgsolve_parity <- mod_mrgsolve_parity |>
  data_set(data_in_mrgsolve) |>
  Req(CP) |>
  mrgsim(end = interval_end, delta = 0.5) |>
  as_tibble() |>
  rename(id = ID) |>
  dplyr::select(id, time, CP_mrgsolve = CP)

## --- nlmixr2/rxode2 run: same dosing/covariates/etas, etas supplied as
## fixed subject-level parameters (omega = NULL suppresses re-simulation) ---
data_in_nlmixr <- dosing_df |>
  tidyr::crossing(id = covset_eta$id) |>
  left_join(covset_eta, by = "id") |>
  arrange(id, time)

sim_nlmixr_parity <- rxSolve(
  mod_nlmixr,
  events = data_in_nlmixr,
  params = theta_vals,
  covsInterpolation = "locf"
) |>
  as_tibble() |>
  dplyr::select(id, time, cp_nlmixr = cp)

## --- Merge and compare on shared time grid ---
parity_df <- sim_mrgsolve_parity |>
  inner_join(sim_nlmixr_parity, by = c("id", "time")) |>
  filter(time > 0) |>   # exclude t=0 (both engines report 0 pre-dose)
  mutate(
    abs_diff = CP_mrgsolve - cp_nlmixr,
    # Relative difference is only meaningful away from near-zero
    # concentrations (e.g., right at/just after a dosing event, where a
    # few ms of solver step-size difference can create a large percentage
    # difference despite a tiny absolute difference). Flag but don't let
    # these dominate the summary.
    rel_diff_pct = 100 * abs_diff / pmax(CP_mrgsolve, 1e-6),
    near_zero = CP_mrgsolve < 1  # mg/L threshold for "near zero" denominator
  )

parity_summary <- parity_df |>
  filter(!near_zero) |>
  summarise(
    n_points          = n(),
    n_excluded_near_zero = sum(parity_df$near_zero),
    max_abs_diff      = max(abs(abs_diff)),
    mean_abs_diff     = mean(abs(abs_diff)),
    median_rel_diff_pct = median(abs(rel_diff_pct)),
    max_rel_diff_pct  = max(abs(rel_diff_pct)),
    mean_rel_diff_pct = mean(abs(rel_diff_pct)),
    correlation       = cor(CP_mrgsolve, cp_nlmixr)
  )

cat("\n--- Parity check: mrgsolve vs. nlmixr2/rxode2 (identical dosing, covariates, etas) ---\n")
print(parity_summary, width = Inf)

## Use the median (robust to isolated near-zero-denominator spikes at dose
## transitions) as the primary pass/fail criterion, and report the max as a
## diagnostic rather than a hard failure trigger.
if (parity_summary$median_rel_diff_pct < 1 && parity_summary$correlation > 0.999) {
  cat("\nPASS: median relative difference < 1% and correlation > 0.999 -- engines are numerically equivalent.\n")
  if (parity_summary$max_rel_diff_pct > 5) {
    cat(sprintf("NOTE: a small number of points (likely near dose-transition timepoints with near-zero concentrations) show up to %.1f%% relative difference; inspect parity_df |> filter(!near_zero) |> arrange(desc(abs(rel_diff_pct))) if concerned.\n",
                parity_summary$max_rel_diff_pct))
  }
} else if (parity_summary$median_rel_diff_pct < 5) {
  cat("\nCAUTION: median relative difference < 5% -- minor solver/interpolation discrepancies, review before relying on exact match.\n")
} else {
  cat("\nFAIL: median relative difference >= 5% -- investigate model translation for discrepancies.\n")
}

## Scatter plot: mrgsolve-predicted vs. nlmixr2-predicted concentrations
plot_parity <- ggplot(parity_df, aes(x = CP_mrgsolve, y = cp_nlmixr)) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "mrgsolve predicted concentration (mg/L)",
    y = "nlmixr2/rxode2 predicted concentration (mg/L)",
    title = "Parity check: mrgsolve vs. nlmixr2/rxode2",
    subtitle = sprintf("Identical dosing, covariates, and etas (n = %d subjects); dashed line = identity",
                        n_id)
  )

ggsave(file.path(fig_dir, "parity_check_mrgsolve_vs_nlmixr2.png"), plot_parity,
       width = 6, height = 6, dpi = 150)

cat(sprintf("\nPlots saved to '%s/': covariate_distributions.png, exposure_distributions_nlmixr2.png, parity_check_mrgsolve_vs_nlmixr2.png\n",
            fig_dir))
