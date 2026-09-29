## ---------------------------------------------------------------------------
## dupilumab_figure5_simulation_m2eta.R
##
## Reruns the Figure 5 simulation (dupilumab_figure5_simulation.R) using the
## Model 4 variant with Model 2 IIV on ka and Vm
## (dupilumab_ad_phase3_model4_m2eta_nlmixr2.R). Dosing, virtual population,
## n_sub and seeds are identical to the original script, so differences
## reflect the added etas only.
## ---------------------------------------------------------------------------

library(nlmixr2)
library(tidyverse)

source("fostvl/dupilumab_ad_phase3_model4_m2eta_nlmixr2.R")

# Run in its own environment so objects don't overwrite the Model 4 results
fig5_m2eta <- new.env()
fig5_m2eta$fig5_model <- dupilumab_ad_model4_m2eta
source("fostvl/dupilumab_figure5_simulation.R", local = fig5_m2eta)

fig5_m2eta$medians |>
  ggplot(aes(time, median_cp, colour = label)) +
  geom_line() +
  geom_vline(xintercept = c(fig5_m2eta$induction_end, fig5_m2eta$maintenance_end),
             linetype = "dashed") +
  scale_x_continuous(breaks = seq(0, fig5_m2eta$sim_end, by = 28)) +
  labs(x = "Time (d)", y = "Concentration (mg/L)", colour = NULL,
       title = "Model 4 + Model 2 IIV on ka and Vm") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(ncol = 2))
