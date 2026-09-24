library(nlmixr2lib)
library(rxode2)
library(dplyr)
library(ggplot2)
library(tibble)
library(purrr)

fig_dir <- "figures"
dir.create(fig_dir, showWarnings = FALSE)

## ---- Steady-state check ----
mod <- readModelDb("Miyano_2022_atopicDermatitis_qsp")
mod_typ <- mod |> rxode2::zeroRe()

ev_ss <- rxode2::et(seq(0, 1000, by = 4)) |>
  rxode2::et(id = 1L)
ev_ss$cmt <- "barrier"

ss <- rxode2::rxSolve(mod_typ, events = ev_ss) |> as.data.frame()
cat("== Steady-state tail ==\n")
print(tail(ss[, c("time", "barrier", "pathogens", "il13", "th2", "EASI")], 5))

p1 <- ss |>
  dplyr::filter(time <= 200) |>
  ggplot(aes(x = time)) +
  geom_line(aes(y = EASI, colour = "EASI")) +
  geom_line(aes(y = 20 * barrier, colour = "20*barrier")) +
  geom_line(aes(y = pathogens, colour = "pathogens")) +
  scale_colour_manual(values = c("EASI" = "black", "20*barrier" = "steelblue", "pathogens" = "firebrick")) +
  labs(x = "Time (weeks)", y = "State value",
       colour = NULL,
       title = "Typical-patient equilibration to steady state",
       caption = "Drift-and-settle from Table-S2 seed IC's to the model's fitted attractor.")
ggsave(file.path(fig_dir, "01_steady_state.png"), p1, width = 7, height = 4.5, dpi = 150)
cat("Saved 01_steady_state.png\n")

## ---- THETA/ETA terminology demo: parameter-level distribution ----
set.seed(123)
n_demo <- 300
lk1 <- 0.537
omega_k1 <- 0.628849  # = 0.793^2; sigma.csv row 1 squared -> etalk1 ~ 0.628849 in ini()

eta_k1 <- rnorm(n_demo, mean = 0, sd = sqrt(omega_k1))
k1_indiv <- exp(lk1 + eta_k1)
k1_typical <- exp(lk1)
cat("Typical-patient k1 =", k1_typical, "\n")

df_theta_eta <- rbind(
  data.frame(quantity = "ETA (etalk1, log scale)", value = eta_k1),
  data.frame(quantity = "individual k1 = exp(lk1 + etalk1)", value = k1_indiv)
)

p_theta_eta <- ggplot(df_theta_eta, aes(x = value)) +
  geom_histogram(bins = 40, fill = "steelblue", alpha = 0.7, colour = "white") +
  facet_wrap(~ quantity, scales = "free", nrow = 1) +
  labs(x = NULL, y = "Count",
       title = "THETA + ETA -> individual parameter (300 simulated virtual patients)",
       subtitle = paste0("ETA ~ N(0, OMEGA = ", round(omega_k1, 3), "); typical-patient k1 = exp(lk1) = ", round(k1_typical, 3)))
ggsave(file.path(fig_dir, "02_theta_eta_distribution.png"), p_theta_eta, width = 8, height = 4, dpi = 150)
cat("Saved 02_theta_eta_distribution.png\n")

## ---- THETA/ETA terminology demo: trajectory-level (spaghetti) ----
set.seed(123)
n_spaghetti <- 12
ev_demo <- rxode2::et(seq(0, 100, by = 2)) |> rxode2::et(id = seq_len(n_spaghetti))
ev_demo$cmt <- "barrier"
ev_demo_typ <- rxode2::et(seq(0, 100, by = 2)) |> rxode2::et(id = 1L)
ev_demo_typ$cmt <- "barrier"

sim_demo <- rxode2::rxSolve(mod, events = ev_demo, nSub = n_spaghetti) |> as.data.frame()
sim_demo_typ <- rxode2::rxSolve(mod_typ, events = ev_demo_typ) |> as.data.frame()

p_spaghetti_demo <- ggplot(sim_demo, aes(x = time, y = EASI, group = id)) +
  geom_line(colour = "grey60", alpha = 0.6) +
  geom_line(data = sim_demo_typ, aes(x = time, y = EASI), inherit.aes = FALSE,
            colour = "black", linewidth = 1.1) +
  labs(x = "Time (weeks)", y = "EASI",
       title = "Same THETAs, different ETAs: 12 virtual patients vs. the typical patient",
       subtitle = "Grey = individual virtual patients (random eta draws); black = typical patient (zeroRe(), eta = 0)",
       caption = "This same eta-sampling mechanism generates between-subject variability in a population PK model.")
ggsave(file.path(fig_dir, "03_spaghetti_patients.png"), p_spaghetti_demo, width = 7, height = 4.5, dpi = 150)
cat("Saved 03_spaghetti_patients.png\n")

## ---- Figure 3B: typical-patient EASI trajectories by drug ----
drug_params <- tibble::tribble(
  ~drug,          ~r_placebo, ~r_il4, ~r_il13, ~r_il17, ~r_il22, ~r_il31, ~r_tslp, ~r_ox40l, ~c_rifng, ~e_a2,
  "Placebo",       1,         0,       0,      0,      0,      0,       0,       0,        0,        1,
  "Dupilumab",     1,         0.99,   0.99,    0,      0,      0,       0,       0,        0,        1,
  "Lebrikizumab",  1,         0,      0.99,    0,      0,      0,       0,       0,        0,        1,
  "Tralokinumab",  1,         0,      0.99,    0,      0,      0,       0,       0,        0,        0.4396,
  "Secukinumab",   1,         0,       0,     0.99,    0,      0,       0,       0,        0,        1,
  "Fezakinumab",   1,         0,       0,      0,     0.99,    0,       0,       0,        0,        1,
  "Nemolizumab",   1,         0,       0,      0,      0,     0.99,     0,       0,        0,        1,
  "Tezepelumab",   1,         0,       0,      0,      0,      0,      0.99,     0,        0,        1,
  "GBR830",        1,         0,       0,      0,      0,      0,       0,      0.99,      0,        1,
  "rIFNg",         1,         0,       0,      0,      0,      0,       0,       0,        210,      1
)

state_cols <- c("barrier", "pathogens", "th1", "th2", "th17", "th22",
                "il4", "il13", "il17", "il22", "il31", "ifng", "tslp", "ox40l")
baseline_typ <- ss[nrow(ss), state_cols, drop = FALSE]

typical_by_drug <- function(row) {
  ev <- rxode2::et(seq(0, 24, by = 0.5)) |> rxode2::et(id = 1L)
  ev$cmt <- "barrier"
  p <- as.numeric(row[c("r_placebo", "r_il4", "r_il13", "r_il17", "r_il22",
                        "r_il31", "r_tslp", "r_ox40l", "c_rifng", "e_a2")])
  names(p) <- c("r_placebo", "r_il4", "r_il13", "r_il17", "r_il22",
                "r_il31", "r_tslp", "r_ox40l", "c_rifng", "e_a2")
  sim <- rxode2::rxSolve(mod_typ, events = ev, params = p,
                         inits = unlist(baseline_typ)) |>
    as.data.frame()
  sim$drug <- row[["drug"]]
  sim
}

typ_all <- purrr::map_dfr(seq_len(nrow(drug_params)), function(i)
  typical_by_drug(drug_params[i, ]))

typ_all <- typ_all |>
  group_by(drug) |>
  mutate(pct_improved_EASI = 100 * (first(EASI) - EASI) / first(EASI)) |>
  ungroup()

p2 <- ggplot(typ_all, aes(x = time, y = pct_improved_EASI, colour = drug)) +
  geom_line(linewidth = 0.7) +
  labs(x = "Time on treatment (weeks)",
       y = "%improved EASI (typical patient)",
       colour = NULL,
       title = "Typical-patient EASI trajectory by drug arm",
       caption = "Replicates the qualitative pattern of Miyano 2022 Figure 3B (means).")
ggsave(file.path(fig_dir, "04_figure3B_drug_trajectories.png"), p2, width = 8, height = 5, dpi = 150)
cat("Saved 04_figure3B_drug_trajectories.png\n")

## ---- EASI-75 virtual cohort ----
n_per_arm <- 50L
mod_stoch <- mod  # keeps eta variability

baseline_typ_named <- as.numeric(baseline_typ[1, ])
names(baseline_typ_named) <- names(baseline_typ)

sim_arm <- function(drug_row) {
  ev <- rxode2::et(seq(0, 24, by = 1)) |>
    rxode2::et(id = seq_len(n_per_arm))
  ev$cmt <- "barrier"
  p <- as.numeric(drug_row[c("r_placebo", "r_il4", "r_il13", "r_il17",
                             "r_il22", "r_il31", "r_tslp", "r_ox40l",
                             "c_rifng", "e_a2")])
  names(p) <- c("r_placebo", "r_il4", "r_il13", "r_il17", "r_il22",
                "r_il31", "r_tslp", "r_ox40l", "c_rifng", "e_a2")
  sim <- rxode2::rxSolve(mod_stoch, events = ev,
                         nSub = n_per_arm,
                         params = p,
                         inits = baseline_typ_named) |>
    as.data.frame()
  sim$drug <- drug_row[["drug"]]
  sim
}

cat("Running cohort simulation (50 patients x 10 arms)...\n")
cohort_all <- purrr::map_dfr(seq_len(nrow(drug_params)), function(i)
  sim_arm(drug_params[i, ]))

baseline_easi <- cohort_all |>
  dplyr::filter(time == 0) |>
  dplyr::select(drug, id, EASI_baseline = EASI)
cohort_all <- cohort_all |>
  dplyr::left_join(baseline_easi, by = c("drug", "id")) |>
  dplyr::mutate(pct_improved = 100 * (EASI_baseline - EASI) / EASI_baseline)

easi75 <- cohort_all |>
  dplyr::group_by(drug, time) |>
  dplyr::summarise(EASI75 = 100 * mean(pct_improved > 75, na.rm = TRUE),
                   .groups = "drop")

p3 <- ggplot(easi75, aes(x = time, y = EASI75, colour = drug)) +
  geom_line(linewidth = 0.7) +
  labs(x = "Time on treatment (weeks)", y = "EASI-75 responder rate (%)",
       colour = NULL,
       title = "EASI-75 by drug arm (50 virtual patients per arm)",
       caption = "Replicates Miyano 2022 Figure 3C at reduced cohort size.")
ggsave(file.path(fig_dir, "05_figure3C_easi75.png"), p3, width = 8, height = 5, dpi = 150)
cat("Saved 05_figure3C_easi75.png\n")

cat("\nAll figures generated successfully in figures/\n")
