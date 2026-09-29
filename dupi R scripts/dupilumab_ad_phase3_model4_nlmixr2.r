## ---------------------------------------------------------------------------
## dupilumab_ad_phase3_model4_nlmixr2.R
##
## nlmixr2 implementation of Model 4 (covariate model parameterized using ke)
## from: Kovalenko P, et al. Base and Covariate Population Pharmacokinetic
## Analyses of Dupilumab Using Phase 3 Data. Clin Pharmacol Drug Dev.
## 2020;9(6):756-767. doi:10.1002/cpdd.780
## (main text Tables 1 and 3; Supplementary Tables S3 and S4).
##
## Structure: 2-compartment, SC absorption through 3 transit compartments
## (MTT) followed by first-order absorption (ka), parallel linear (ke) and
## Michaelis-Menten (Vm in mg/L/day, Km) elimination from the central
## compartment. IV doses go to `central`.
##
## Checked against the paper / supplement:
##   - Estimates, fixed parameters, IIV (SD 0.206 Vc, 0.293 ke, corr -0.450)
##     and residual error (12.5% prop, 6.06 mg/L add) match Suppl. Table S3/S4.
##   - Continuous covariates: (COV/ref)^theta; dichotomous: exp(theta * COV)
##     (Methods, main text).
##
## ASSUMPTIONS / GAPS (not reported in paper or supplement):
##   1. Covariate reference values are assumed to be the phase 3 popPK
##      dataset medians (WT 73 kg, ALB 45 g/L, BMI 24.9, EASI 29.3; FDA Clin
##      Pharm Review BLA 761055, Table 4.6.1). The paper says "median or
##      another selected level" but does not state the values.
##   2. ADA was reported as a significant covariate on ke, but its coefficient
##      is not reported anywhere (authors deliberately withheld it); omitted.
##   3. Transit rate follows the Monolix convention ktr = (n + 1) / MTT, n = 3.
## ---------------------------------------------------------------------------

library(nlmixr2)

dupilumab_ad_model4 <- function() {
  ini({
    # Estimated fixed effects
    lvc  <- log(2.74)    ; label("log Vc (L)")
    lke  <- log(0.0477)  ; label("log ke (1/day)")

    # Covariate effects
    vc_wt    <-  0.817   ; label("Vc ~ weight (power)")
    vc_alb   <- -0.653   ; label("Vc ~ albumin (power)")
    ke_bmi   <-  0.368   ; label("ke ~ BMI (power)")
    ke_easi  <-  0.143   ; label("ke ~ EASI (power)")
    ke_white <- -0.123   ; label("ke ~ race White (exponential)")

    # Fixed parameters (fixed in the paper)
    lkcp  <- fix(log(0.211)) ; label("log kcp (1/day)")
    lkpc  <- fix(log(0.310)) ; label("log kpc (1/day)")
    lka   <- fix(log(0.306)) ; label("log ka (1/day)")
    lmtt  <- fix(log(0.105)) ; label("log MTT (day)")
    lvm   <- fix(log(1.07))  ; label("log Vm (mg/L/day)")
    lkm   <- fix(log(0.01))  ; label("log Km (mg/L)")
    lfsc  <- fix(logit(0.642)) ; label("logit F")

    # IIV: SDs 0.206 (Vc), 0.293 (ke), correlation -0.450
    eta.vc + eta.ke ~ c(0.042436,
                        -0.027161, 0.085849)

    prop.sd <- 0.125     ; label("Proportional residual SD")
    add.sd  <- 6.06      ; label("Additive residual SD (mg/L)")
  })
  model({
    # Reference values are placeholders; hardcoded so rxode2 doesn't treat them as data columns
    vc <- exp(lvc + eta.vc) * (WT / 73)^vc_wt * (ALB / 45)^vc_alb
    ke <- exp(lke + eta.ke) * (BMI / 24.9)^ke_bmi * (EASI / 29.3)^ke_easi *
          exp(ke_white * WHITE)
    kcp <- exp(lkcp)
    kpc <- exp(lkpc)
    ka  <- exp(lka)
    vm  <- exp(lvm)
    km  <- exp(lkm)
    # Monolix convention for 3 transit compartments: ktr = (n + 1) / MTT
    ktr <- 4 / exp(lmtt)

    cp <- central / vc

    d/dt(depot)      <- -ktr * depot
    d/dt(transit1)   <-  ktr * depot    - ktr * transit1
    d/dt(transit2)   <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)   <-  ktr * transit2 - ka  * transit3
    d/dt(central)    <-  ka * transit3 - ke * central - kcp * central +
                         kpc * peripheral - vm * vc * cp / (km + cp)
    d/dt(peripheral) <-  kcp * central - kpc * peripheral

    f(depot) <- expit(lfsc)

    cp ~ add(add.sd) + prop(prop.sd)
  })
}
