## ---------------------------------------------------------------------------
## vpc_nlmixr2.R
##
## Simulation-based "VPC-like" figure for the nlmixr2/rxode2 dupilumab
## asthma PopPK model, for comparison against Figure 2 (visual predictive
## check) of:
##
## Zhang L, Gao Y, Li M, Xu C, Davis JD, Kanamaluru V, Lu Q.
## Population pharmacokinetic analysis of dupilumab in adult and adolescent
## patients with asthma. CPT Pharmacometrics Syst Pharmacol. 2021;10(8):941-952.
##
## IMPORTANT LIMITATION: A true VPC overlays observed patient concentrations
## against simulated prediction intervals to check whether the model
## reproduces the *distribution* of real data. We do not have access to the
## underlying patient-level PK dataset used in the manuscript (it is not
## publicly available), so this script produces a **simulation-only
## prediction interval plot** (median + 90% simulated interval over time)
## using the nlmixr2 model with IIV and residual error, for the same 300 mg
## q2w regimen (with 600 mg loading dose) as Figure 2's phase III studies.
## This can be compared *qualitatively* (median trajectory shape, plateau
## level, width of variability band) against the published VPC, but it is
## not a genuine VPC because no observed data are overlaid.
## ---------------------------------------------------------------------------

library(nlmixr2)
library(rxode2)
library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(3059)

fig_dir <- "figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir)

## -----------------------------------------------------------------------
## 1. Load model
## -----------------------------------------------------------------------
source("dupilumab_asthma_poppk_nlmixr2.R")
mod_nlmixr <- dupilumab_asthma_poppk()

theta_vals <- c(
  tke = log(0.0418), tv2 = log(2.76), tk23 = log(0.0952), tk32 = log(0.163),
  tvmax = log(1.39), tkm = log(2.08), tka = log(0.263), tfsc = log(0.609 / (1 - 0.609)),
  wtke = 0.222, adake = 0.191, clcrnke = 0.217, wtv2 = 0.667, albv2 = -0.484, wtvmax = 0.224,
  prop.err = 0.0388, add.err = 1.73
)

## -----------------------------------------------------------------------
## 2. Quantile-based covariate sampler calibrated to Table 2 statistics
##    (identical to simulate_dupilumab_population_quantile.R /
##    simulate_dupilumab_population_nlmixr2.R)
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

n_id <- 1000

covset <- tibble(
  id    = seq_len(n_id),
  WT    = sample_from_table2(n_id, mean_val = 79.6,  sd_val = 19.1, min_val = 32,   max_val = 186),
  ALB   = sample_from_table2(n_id, mean_val = 43.8,  sd_val = 3.50, min_val = 30,   max_val = 55),
  CLCRN = sample_from_table2(n_id, mean_val = 116.2, sd_val = 36.0, min_val = 30.1, max_val = 377),
  ADA   = rbinom(n_id, size = 1, prob = 344 / 2114)
)

## -----------------------------------------------------------------------
## 3. Dosing regimen: 300 mg q2w with 600 mg loading dose
##    (single loading dose at time 0; maintenance starts explicitly at
##    day 14 to avoid the duplicate-dose/mistimed-sequence bug found
##    earlier in this project)
## -----------------------------------------------------------------------
n_doses_maintenance <- 25
sample_times <- seq(0, 14 * (n_doses_maintenance + 1) + 14, by = 1)

ev_template <- eventTable() |>
  add.dosing(dose = 600, start.time = 0, nbr.doses = 1, dosing.to = "depot") |>
  add.dosing(dose = 300, start.time = 14, nbr.doses = n_doses_maintenance,
             dosing.interval = 14, dosing.to = "depot") |>
  add.sampling(sample_times)

dosing_df <- as.data.frame(ev_template)

data_in <- dosing_df |>
  tidyr::crossing(id = covset$id) |>
  left_join(covset, by = "id") |>
  arrange(id, time)

## -----------------------------------------------------------------------
## 4. Simulate population with IIV + residual error (i.e., simulated
##    "observations", analogous to what a VPC's simulated replicates
##    represent) using rxSolve's built-in residual-error simulation.
## -----------------------------------------------------------------------
sim <- rxSolve(
  mod_nlmixr,
  events = data_in,
  params = theta_vals,
  omega  = mod_nlmixr$omega,
  sigma  = mod_nlmixr$sigma,
  covsInterpolation = "locf"
) |>
  as_tibble()

## -----------------------------------------------------------------------
## 5. Bin nominal time and compute median / 5th / 95th percentiles of the
##    simulated ("observed-like", i.e. including residual error) and
##    typical ("ipred", i.e. individual-prediction without residual error)
##    concentrations at each time bin -- mirroring the median/90% PI bands
##    shown in the manuscript's Figure 2.
##
##    NOTE: because every simulated subject shares the same dosing times
##    (synchronized q2w schedule) and the simulation grid is already dense
##    (daily), we summarize at the *daily* resolution rather than coarser
##    weekly bins -- coarser binning here would straddle peak/trough phases
##    of the q2w cycle inconsistently and distort the band shape (verified:
##    7-day bins produced a spurious sawtooth artifact in the median).
## -----------------------------------------------------------------------
vpc_summary <- sim |>
  group_by(time_bin = time) |>
  summarise(
    median_sim = median(sim, na.rm = TRUE),
    lo90_sim   = quantile(sim, 0.05, na.rm = TRUE),
    hi90_sim   = quantile(sim, 0.95, na.rm = TRUE),
    median_ipred = median(ipredSim, na.rm = TRUE),
    lo90_ipred   = quantile(ipredSim, 0.05, na.rm = TRUE),
    hi90_ipred   = quantile(ipredSim, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(time_bin >= 0)

## -----------------------------------------------------------------------
## 6. Plot: simulation-based VPC-like figure
##
## Solid line + dark ribbon: median and 90% prediction interval of
## simulated "observations" (with residual error) -- the standard VPC
## bands. Dashed line: median individual-prediction trajectory (no
## residual error) for reference.
## -----------------------------------------------------------------------
plot_vpc <- ggplot(vpc_summary, aes(x = time_bin)) +
  geom_ribbon(aes(ymin = pmax(lo90_sim, 0), ymax = hi90_sim), fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = median_sim), color = "steelblue", linewidth = 0.9) +
  geom_line(aes(y = median_ipred), color = "black", linetype = "dashed", linewidth = 0.6) +
  labs(
    x = "Nominal time (days)",
    y = "Dupilumab concentration (mg/L)",
    title = "Simulation-based VPC-like prediction interval (nlmixr2/rxode2 model)",
    subtitle = paste0(
      "300 mg q2w (600 mg loading dose), n = ", n_id, " simulated subjects\n",
      "Shaded band = simulated 5th-95th percentile; solid line = simulated median;\n",
      "dashed line = median individual prediction (no residual error).\n",
      "NOTE: no observed data overlaid -- not a true VPC (see script header)."
    )
  ) +
  theme(plot.subtitle = element_text(size = 8))

ggsave(file.path(fig_dir, "vpc_like_nlmixr2.png"), plot_vpc, width = 8, height = 6, dpi = 150)

print(vpc_summary, n = Inf)
cat(sprintf("\nSteady-state (days 200-378) median simulated trough/peak summary:\n"))
vpc_summary |>
  filter(time_bin >= 200) |>
  summarise(
    median_trough_approx = min(median_sim),
    median_peak_approx   = max(median_sim)
  ) |>
  print()

cat(sprintf("\nPlot saved to '%s/vpc_like_nlmixr2.png'\n", fig_dir))

## -----------------------------------------------------------------------
## 7. Fetch and save the manuscript's actual VPC figure (Figure 2) for
##    direct visual comparison alongside the simulation-based plot above.
## -----------------------------------------------------------------------
fig2_url <- "https://cdn.ncbi.nlm.nih.gov/pmc/blobs/8e14/8376131/d022a091956d/PSP4-10-941-g001.jpg"
fig2_dest <- file.path(fig_dir, "manuscript_figure2_vpc.jpg")
tryCatch(
  download.file(fig2_url, destfile = fig2_dest, mode = "wb", quiet = TRUE),
  error = function(e) message("Could not download manuscript Figure 2: ", conditionMessage(e))
)
