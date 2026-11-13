# ACTG 175 Survival Analysis in SAS (with R validation)

Survival analysis of AIDS Clinical Trials Group Study 175, a randomized trial of 2,139 adults with HIV comparing four antiretroviral regimens: zidovudine (ZDV) only, ZDV + didanosine (ddI), ZDV + zalcitabine (Zal), and ddI only. The analysis is written in SAS, and the key models are reproduced in R.

Project page: https://chunleige.github.io/portfolio/actg175-survival.html

## Endpoint

Composite endpoint: a ≥50% decline in CD4 count, progression to AIDS, or death (`cens = 1`). Time to event or censoring in days (`days`).

## Repository structure

```
actg175-survival-sas/
├── README.md
├── data/
│   └── README.md            # how to download the data (raw CSV not committed)
├── sas/
│   ├── 01_import.sas            # import CSV, structural checks
│   ├── 02_prep_table1.sas       # derived variables, formats, Table 1 (ODS RTF)
│   ├── 03_km_logrank.sas        # Kaplan–Meier curves, log-rank, pairwise tests
│   ├── 04_cox_phreg.sas         # Cox models, Efron vs. Breslow, PH assessment, stratified model
│   ├── 05_subgroup_forest.sas   # subgroup macro, interaction tests, forest plot
│   ├── 06_cd4_change.sas        # ANCOVA of week-20 CD4 change
│   └── 08_rmst.sas              # restricted mean survival time
├── R/
│   ├── 07_r_validation.R        # Cox models, log-rank, ANCOVA reproduced in R
│   └── 08_rmst.R                # RMST cross-check with survRM2
├── figures/
│   ├── km_by_arm.png
│   ├── forest_plot.png
│   ├── cd4_slopes.png
│   ├── ph_cd40.png
│   └── ph_arms1.png
└── output/
    ├── table1.rtf
    └── sas_vs_r_comparison.csv
```

## How to run

**SAS** (SAS OnDemand for Academics or any SAS 9.4 with SAS/STAT 15.1+):

1. Create a folder `~/actg175` and upload the raw CSV as `actg175.csv`.
2. Run the scripts in `sas/` in numeric order. Script 02 stores formats in the `proj` library, so later scripts can be run on their own in a new session.

**R** (4.x, packages `survival` and `survRM2`):

1. Set `path` at the top of each script to the location of the CSV.
2. Run `07_r_validation.R`, then `08_rmst.R`.

## Main results

| Regimen vs. ZDV only | Adjusted HR (95% CI) | RMST gain, τ = 36 mo (95% CI) |
|---|---|---|
| ZDV + ddI | 0.45 (0.36–0.58) | +3.6 months (2.5–4.6) |
| ZDV + Zal | 0.51 (0.41–0.65) | +3.4 months (2.3–4.5) |
| ddI only | 0.57 (0.45–0.71) | +2.8 months (1.7–4.0) |

SAS and R estimates agree to at least five decimal places for all Cox models, the log-rank test, the ANCOVA and RMST.

## Notes

- `PROC PHREG` defaults to Breslow ties and `coxph()` defaults to Efron; both are set explicitly here.
- Baseline CD4 violated proportional hazards and was handled with a model stratified by baseline CD4 quartile.
- Subgroup analyses are exploratory.
- ACTG 175 is a randomized trial; results should not be read as observational registry findings.

## Data source and citation

Hammer SM, Katzenstein DA, Hughes MD, et al. A trial comparing nucleoside monotherapy with combination therapy in HIV-infected adults with CD4 cell counts from 200 to 500 per cubic millimeter. *N Engl J Med*. 1996;335(15):1081–1090.

Dataset: UCI Machine Learning Repository, ID 890 (AIDS Clinical Trials Group Study 175). Please follow the acknowledgment policy on the original dataset website.
