[PROB]
# Parity-check variant of dupilumab_asthma_poppk.cpp
#
# Identical structural/covariate model, but etas are declared as ordinary
# $PARAM inputs (not $OMEGA-simulated) so that externally supplied
# per-subject eta values (e.g., shared with an nlmixr2/rxode2 parity run)
# are used directly rather than being resimulated internally by mrgsolve.
#
# See dupilumab_asthma_poppk.cpp for the full annotated/production model.

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

[PARAM] @annotated
ETA_KE   : 0 : Individual eta on Ke (supplied externally for parity check)
ETA_V2   : 0 : Individual eta on V2 (supplied externally for parity check)
ETA_VMAX : 0 : Individual eta on Vmax (supplied externally for parity check)
ETA_KA   : 0 : Individual eta on Ka (supplied externally for parity check)
ETA_FSC  : 0 : Individual eta on Fsc (supplied externally for parity check)

[CMT] @annotated
GUT   : Absorption depot (mg)
CENT  : Central compartment (mg) [ADM, OBS]
PERIPH: Peripheral compartment (mg)

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

F_GUT = Fsc;

[ODE]
double CONC = CENT / V2;
double MM = (VMAX * CONC) / (KM + CONC);

dxdt_GUT    = -KA * GUT;
dxdt_CENT   =  KA * GUT - KE * CENT - MM * V2 - K23 * CENT + K32 * PERIPH;
dxdt_PERIPH =  K23 * CENT - K32 * PERIPH;

[TABLE]
double CP = CENT / V2;

[CAPTURE] @annotated
CP   : Central compartment concentration (mg/L)
