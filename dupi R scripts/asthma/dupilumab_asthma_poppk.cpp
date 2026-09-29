[PROB]
# Dupilumab Asthma Population PK Model
#
# Reference:
# Zhang L, Gao Y, Li M, Xu C, Davis JD, Kanamaluru V, Lu Q.
# Population pharmacokinetic analysis of dupilumab in adult and adolescent
# patients with asthma. CPT Pharmacometrics Syst Pharmacol. 2021;10(8):941-952.
# https://pmc.ncbi.nlm.nih.gov/articles/PMC8376131/
#
# Structure: two-compartment model with first-order absorption (s.c.) and
# parallel linear + nonlinear (Michaelis-Menten) elimination from the
# central compartment.
#
# Covariates (final model, Table 3):
#   V2   = 2.76  * (WT/78)^0.667  * (ALB/44)^-0.484 * exp(ETA_V2)
#   Vmax = 1.39  * (WT/78)^0.224  * exp(ETA_VMAX)
#   Ke   = 0.0418* (1 + 0.191*ADA) * (WT/78)^0.222 * (CLCRN/111)^0.217 * exp(ETA_KE)
#   Ka   = 0.263 * exp(ETA_KA)
#   Fsc  = 0.609 * exp(ETA_FSC)   (logit-normal on 0-1 scale, see $MAIN)
#
# Fixed (no covariates identified) inter-compartmental parameters:
#   K23 = 0.0952 1/day, K32 = 0.163 1/day
#
# Residual error: combined proportional + additive on concentration (mg/L)
#   prop CV 19.7%(var 0.0388), additive SD 1.73 mg/L (var 2.98 (mg/L)^2 as reported;
#   note: table reports variance-like values for resid; see $SIGMA below.)
#
# Doses: s.c. administration in mg into the GUT (depot) compartment;
# central/peripheral volumes in liters, amounts in mg -> concentration mg/L.

[PARAM] @annotated
WT     : 78   : Body weight (kg)
ALB    : 44   : Albumin (g/L)
CLCRN  : 111  : Creatinine clearance normalized to BSA (mL/min/1.73m2)
ADA    : 0    : Positive ADA status (0 = negative, 1 = positive)

[PARAM] @annotated
TVKE   : 0.0418 : Typical linear elimination rate constant Ke (1/day)
TVV2   : 2.76   : Typical central volume V2 (L)
TVK23  : 0.0952 : Intercompartmental rate constant K23 (1/day)
TVK32  : 0.163  : Intercompartmental rate constant K32 (1/day)
TVVMAX : 1.39   : Maximum target-mediated elimination rate Vmax (mg/L/day)
TVKM   : 2.08   : Michaelis constant Km (mg/L)
TVKA   : 0.263  : Absorption rate constant Ka (1/day)
TVFSC  : 0.609  : Typical bioavailability Fsc (fraction)
WTKE   : 0.222  : Power effect of weight on Ke
ADAKE  : 0.191  : Proportional effect of positive ADA on Ke
CLCRNKE: 0.217  : Power effect of CLCRN on Ke
WTV2   : 0.667  : Power effect of weight on V2
ALBV2  : -0.484 : Power effect of albumin on V2
WTVMAX : 0.224  : Power effect of weight on Vmax

[CMT] @annotated
GUT   : Absorption depot (mg)
CENT  : Central compartment (mg) [ADM, OBS]
PERIPH: Peripheral compartment (mg)

[OMEGA] @annotated @name IIV
ETA_KE   : 0.0385  : IIV on Ke
ETA_V2   : 0.00834 : IIV on V2
ETA_VMAX : 0.0589  : IIV on Vmax
ETA_KA   : 0.243   : IIV on Ka
ETA_FSC  : 0.132   : IIV on Fsc (logit scale)

[SIGMA] @annotated
PROP : 0.0388 : Proportional residual error variance
ADD  : 2.9929 : Additive residual error variance (mg/L)^2 (SD 1.73 mg/L)

[MAIN]
double KE   = TVKE   * (1 + ADAKE * ADA) * pow(WT/78.0, WTKE) * pow(CLCRN/111.0, CLCRNKE) * exp(ETA_KE);
double V2   = TVV2   * pow(WT/78.0, WTV2) * pow(ALB/44.0, ALBV2) * exp(ETA_V2);
double VMAX = TVVMAX * pow(WT/78.0, WTVMAX) * exp(ETA_VMAX);
double KM   = TVKM;
double K23  = TVK23;
double K32  = TVK32;
double KA   = TVKA * exp(ETA_KA);

// Fsc modeled on logit scale to keep bioavailability bounded in (0,1)
double logit_Fsc = log(TVFSC/(1.0 - TVFSC)) + ETA_FSC;
double Fsc = 1.0 / (1.0 + exp(-logit_Fsc));

double V3 = V2 * K23 / K32;   // peripheral volume, derived
F_GUT = Fsc;

[ODE]
double CONC = CENT / V2;
double MM = (VMAX * CONC) / (KM + CONC);

dxdt_GUT    = -KA * GUT;
dxdt_CENT   =  KA * GUT - KE * CENT - MM * V2 - K23 * CENT + K32 * PERIPH;
dxdt_PERIPH =  K23 * CENT - K32 * PERIPH;

[TABLE]
double CP = CENT / V2;
double DV = CP * (1.0 + EPS(1)) + EPS(2);

[CAPTURE] @annotated
CP   : Central compartment concentration (mg/L)
KE   : Individual linear elimination rate (1/day)
V2   : Individual central volume (L)
VMAX : Individual Vmax (mg/L/day)
KA   : Individual Ka (1/day)
Fsc  : Individual bioavailability
DV   : Simulated observed concentration (mg/L)
