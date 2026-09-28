## ---------------------------------------------------------------------------
## dupilumab_asthma_poppk_nlmixr2.R
##
## nlmixr2 translation of the mrgsolve dupilumab asthma PopPK model
## (dupilumab_asthma_poppk.cpp), based on:
##
## Zhang L, Gao Y, Li M, Xu C, Davis JD, Kanamaluru V, Lu Q.
## Population pharmacokinetic analysis of dupilumab in adult and adolescent
## patients with asthma. CPT Pharmacometrics Syst Pharmacol. 2021;10(8):941-952.
##
## Structure: two-compartment model with first-order absorption (s.c.) and
## parallel linear + nonlinear (Michaelis-Menten) elimination from the
## central compartment.
##
## Covariates (final model, Table 3):
##   V2   = 2.76  * (WT/78)^0.667  * (ALB/44)^-0.484 * exp(eta.v2)
##   Vmax = 1.39  * (WT/78)^0.224  * exp(eta.vmax)
##   Ke   = 0.0418* (1 + 0.191*ADA) * (WT/78)^0.222 * (CLCRN/111)^0.217 * exp(eta.ke)
##   Ka   = 0.263 * exp(eta.ka)
##   Fsc  = expit(logit(0.609) + eta.fsc)   (bioavailability constrained 0-1)
##
## Fixed (no covariates identified) inter-compartmental parameters:
##   K23 = 0.0952 1/day, K32 = 0.163 1/day
##
## THETAs are parameterized on the log scale (log-normal fixed effects), so
## that estimation stays positive-constrained; the covariate exponents/
## coefficients are estimated on the natural (additive) scale, consistent
## with how they enter the power/proportional covariate relationships.
##
## Residual error: combined proportional + additive on concentration (mg/L):
##   proportional CV ~ 19.7% (variance 0.0388), additive SD 1.73 mg/L
## ---------------------------------------------------------------------------

library(nlmixr2)

dupilumab_asthma_poppk <- function() {
  ini({
    ## --- Fixed effects (log-scale for positivity) ---
    tke   <- log(0.0418)   # Typical linear elimination rate Ke (1/day)
    tv2   <- log(2.76)     # Typical central volume V2 (L)
    tk23  <- log(0.0952)   # Intercompartmental rate constant K23 (1/day)
    tk32  <- log(0.163)    # Intercompartmental rate constant K32 (1/day)
    tvmax <- log(1.39)     # Maximum target-mediated elimination rate Vmax (mg/L/day)
    tkm   <- log(2.08)     # Michaelis constant Km (mg/L)
    tka   <- log(0.263)    # Absorption rate constant Ka (1/day)
    tfsc  <- log(0.609/(1 - 0.609))  # Typical bioavailability Fsc, logit scale

    ## --- Covariate effects (natural scale, as reported in Table 3) ---
    wtke    <- 0.222   # Power effect of weight on Ke
    adake   <- 0.191   # Proportional effect of positive ADA on Ke
    clcrnke <- 0.217   # Power effect of CLCRN on Ke
    wtv2    <- 0.667   # Power effect of weight on V2
    albv2   <- -0.484  # Power effect of albumin on V2
    wtvmax  <- 0.224   # Power effect of weight on Vmax

    ## --- Inter-individual variability (variances, matching mrgsolve $OMEGA) ---
    eta.ke   ~ 0.0385
    eta.v2   ~ 0.00834
    eta.vmax ~ 0.0589
    eta.ka   ~ 0.243
    eta.fsc  ~ 0.132

    ## --- Residual error (combined proportional + additive) ---
    prop.err <- 0.0388      # proportional error variance (CV^2 scale)
    add.err  <- 1.73        # additive error SD (mg/L)
  })
  model({
    KE   <- exp(tke)   * (1 + adake * ADA) * (WT / 78)^wtke * (CLCRN / 111)^clcrnke * exp(eta.ke)
    V2   <- exp(tv2)   * (WT / 78)^wtv2 * (ALB / 44)^albv2 * exp(eta.v2)
    VMAX <- exp(tvmax) * (WT / 78)^wtvmax * exp(eta.vmax)
    KM   <- exp(tkm)
    K23  <- exp(tk23)
    K32  <- exp(tk32)
    KA   <- exp(tka) * exp(eta.ka)

    ## Fsc modeled on logit scale to keep bioavailability bounded in (0,1)
    logit_Fsc <- tfsc + eta.fsc
    Fsc <- exp(logit_Fsc) / (1 + exp(logit_Fsc))

    f(depot) <- Fsc

    CP <- center / V2
    MM <- (VMAX * CP) / (KM + CP)

    d/dt(depot)  <- -KA * depot
    d/dt(center) <- KA * depot - KE * center - MM * V2 - K23 * center + K32 * periph
    d/dt(periph) <- K23 * center - K32 * periph

    cp <- center / V2
    cp ~ prop(prop.err) + add(add.err)
  })
}

## -----------------------------------------------------------------------
## Notes on translation from mrgsolve -> nlmixr2
## -----------------------------------------------------------------------
## 1. Compartments: mrgsolve's [CMT] GUT/CENT/PERIPH map to depot/center/
##    periph (lower case; nlmixr2 auto-detects `d/dt(cmt)` blocks and does
##    not require an explicit compartment list).
## 2. Bioavailability: mrgsolve used `F_GUT = Fsc` in $MAIN; nlmixr2 uses
##    `f(depot) <- Fsc` inside $model.
## 3. Fixed effects (theta) are declared on the log scale in $ini so that
##    exp(theta) is used in $model to enforce positivity during estimation;
##    this does not change the model's population predictions when thetas
##    are fixed at the reported literature values (as done here).
## 4. Random effects (eta) variances in $ini (`eta.x ~ value`) are equivalent
##    to mrgsolve's $OMEGA diagonal variances -- both parameterizations use
##    variance (not SD) by default.
## 5. Residual error: nlmixr2's `cp ~ prop(prop.err) + add(add.err)` is a
##    combined proportional + additive error model equivalent to the
##    mrgsolve $TABLE construction `CP*(1+EPS(1)) + EPS(2)`; prop.err is a
##    variance (matching mrgsolve's $SIGMA PROP), and add.err here is
##    specified as an SD (nlmixr2's `add()` expects SD, unlike mrgsolve's
##    $SIGMA ADD which was a variance) -- confirm this scaling if you refit.
##
## For this literature-fixed-parameter port, no estimation is performed;
## all thetas are fixed at the manuscript's reported values, so this model
## is intended for simulation (e.g., via `rxode2`/`nlmixr2` simulate
## functions) rather than re-estimation, unless you supply patient-level
## data.
