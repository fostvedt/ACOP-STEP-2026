## ---------------------------------------------------------------------------
## dupilumab_figure5_simulation.R
##
## Recreates Figure 5 of Kovalenko et al. 2020 (doi:10.1002/cpdd.780):
## median simulated functional dupilumab concentration over time for
## 16-week induction (300 mg qw or q2w after 600 mg load) followed by
## maintenance to week 52 (qw, q2w, q4w, q8w, placebo), then washout to day 560.
##
## Caveats:
##   - Covariates are sampled to match the popPK phase 3 dataset summary in
##     the FDA Clin Pharm Review (BLA 761055, Table 4.6.1). The paper's
##     own simulated population is not reported. n_sub is kept small (25)
##     for speed; increase (e.g. 1000) for stable medians. IIV on Vc and ke is included; residual error is not.
##   - Model 4 without the ADA term (coefficient unpublished).
##   - Dosing times are inferred from the figure (induction days 0-111,
##     maintenance days 112-363).
## ---------------------------------------------------------------------------

library(nlmixr2)
library(tidyverse)

source("fostvl/dupilumab_ad_phase3_model4_nlmixr2.R")

set.seed(7318)
n_sub <- 250

induction_end   <- 112
maintenance_end <- 364
sim_end         <- 560

dose_times <- function(start, end, interval) {
  if (is.na(interval)) return(numeric(0))
  seq(start, end - 1, by = interval)
}

regimens <- tribble(
  ~induction, ~maintenance,
  "q1w",      "q1w",
  "q1w",      "q4w",
  "q1w",      "q8w",
  "q1w",      "Placebo",
  "q2w",      "q2w",
  "q2w",      "q4w",
  "q2w",      "q8w",
  "q2w",      "Placebo"
) |>
  mutate(
    label = paste0("300 mg ", induction, " SC/",
                   if_else(maintenance == "Placebo", "Placebo",
                           paste0("300 mg ", maintenance, " SC"))),
    label = fct_inorder(label)
  )

interval_days <- c(q1w = 7, q2w = 14, q4w = 28, q8w = 56, Placebo = NA)

make_events <- function(induction, maintenance) {
  maint_doses <- dose_times(induction_end, maintenance_end,
                            interval_days[[maintenance]])
  ind_doses <- dose_times(interval_days[[induction]], induction_end,
                          interval_days[[induction]])
  doses <- tibble(
    time = c(0, ind_doses, maint_doses),
    amt  = c(600, rep(300, length(ind_doses) + length(maint_doses)))
  )
  obs_times <- seq(0, sim_end, by = 1)

  bind_rows(
    doses |> mutate(evid = 1, cmt = "depot"),
    tibble(time = obs_times, amt = 0, evid = 0, cmt = "central")
  ) |>
    arrange(time, desc(evid))
}

# Virtual phase 3 AD population matching the popPK phase 3 dataset (N = 897).
# Source: FDA Clinical Pharmacology Review, BLA 761055 (2017), Tables 4.6.1
# and 4.6.3 (from the Applicant's popPK report, Tables 5 and 9):
#
#   Covariate  Mean  SD    Min   P5    Q1    Median  Q3    P95   Max
#   WT (kg)    75.7  18.6  41.0  50.5  62.3  73.0    85.0  109   175
#   BMI        26.0  5.65  16.3  19.2  22.0  24.9    28.4  37.0  57.3
#   ALB (g/L)  44.6  3.69  22.0  39.0  42.0  45.0    47.0  50.0  57.0
#   EASI       32.1  12.9  0.6   16.8  21.3  29.3    40.9  56.4  72.0
#   White: 615/897 (68.6%)
#
# Distributional choices (assumptions):
#   - WT and BMI: bivariate log-normal fitted to median/IQR; the WT-BMI
#     correlation is not reported and is assumed to be 0.85.
#   - ALB: normal(44.6, 3.69), truncated at the observed range.
#   - EASI: shifted gamma above the entry cutoff of 16 matching mean/SD.
cov_ranges <- list(WT = c(41, 175), BMI = c(16.3, 57.3), ALB = c(22, 57),
                   EASI = c(16, 72))
wt_mu   <- log(73.0); wt_sd  <- log(85.0 / 62.3) / 1.349
bmi_mu  <- log(24.9); bmi_sd <- log(28.4 / 22.0) / 1.349
wt_bmi_rho <- 0.85
easi_excess_mean <- 32.1 - 16
easi_shape <- (easi_excess_mean / 12.9)^2
p_white <- 615 / 897

clamp <- function(x, r) pmin(pmax(x, r[1]), r[2])

sample_covariates <- function(n) {
  z1 <- rnorm(n)
  z2 <- wt_bmi_rho * z1 + sqrt(1 - wt_bmi_rho^2) * rnorm(n)
  tibble(
    id    = seq_len(n),
    WT    = clamp(exp(wt_mu + wt_sd * z1), cov_ranges$WT),
    BMI   = clamp(exp(bmi_mu + bmi_sd * z2), cov_ranges$BMI),
    ALB   = clamp(rnorm(n, 44.6, 3.69), cov_ranges$ALB),
    EASI  = clamp(16 + rgamma(n, easi_shape, easi_shape / easi_excess_mean),
                  cov_ranges$EASI),
    WHITE = rbinom(n, 1, p_white)
  )
}

covs <- sample_covariates(n_sub)

# Callers may supply `fig5_model` (e.g. via source(local = env)); default Model 4
if (!exists("fig5_model", inherits = FALSE)) fig5_model <- dupilumab_ad_model4
mod <- rxode2::zeroRe(nlmixr(fig5_model), which = "sigma")

# Same virtual patients (covariates and etas) across all regimens
sim_one <- function(induction, maintenance, label) {
  ev <- make_events(induction, maintenance) |>
    cross_join(covs) |>
    arrange(id, time, desc(evid))
  set.seed(2961)
  rxode2::rxSolve(mod, ev, returnType = "tibble") |>
    select(id, time, cp) |>
    mutate(label = label)
}

sims <- pmap(select(regimens, induction, maintenance, label), sim_one) |>
  list_rbind()

medians <- sims |>
  group_by(label, time) |>
  summarise(median_cp = median(cp), .groups = "drop")

ggplot(medians, aes(time, median_cp, colour = label)) +
  geom_line() +
  geom_vline(xintercept = c(induction_end, maintenance_end), linetype = "dashed") +
  scale_x_continuous(breaks = seq(0, sim_end, by = 28)) +
  labs(x = "Time (d)", y = "Concentration (mg/L)", colour = NULL) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(ncol = 2))
