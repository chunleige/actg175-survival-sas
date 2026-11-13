/*-----------------------------------------------------------------
  Project : ACTG 175 survival analysis (portfolio project)
  Script  : 03_km_logrank.sas
  Purpose : Kaplan-Meier curves by treatment arm, log-rank test,
            and pairwise comparisons vs. the ZDV-only reference arm
  Input   : proj.actg
  Output  : KM plot (km_arms.png), survival estimates, test results
  Author  : Ivyann
  Date    : 2026-11-11
  Notes   : - Endpoint: composite (CD4 decline >=50%, AIDS, or death).
            - censor = 1 means censored (see 02_prep_table1.sas).
            - Time in months for plotting; analysis is equivalent in days.
------------------------------------------------------------------*/

libname proj "~/actg175";
options fmtsearch=(proj);

ods listing gpath="~/actg175" image_dpi=300;
ods graphics on / imagename="km_arms" imagefmt=png width=8in height=6in;

/* 1. KM curves with at-risk table and confidence bands */
proc lifetest data=proj.actg notable
              plots=survival(atrisk(outside) cl);
  time time_mo*censor(1);
  strata arms / test=logrank;
  title "Kaplan-Meier: freedom from composite endpoint by arm";
run;

/* 2. Pairwise log-rank vs. ZDV only (no plot) */
proc lifetest data=proj.actg notable plots=none;
  time time_mo*censor(1);
  strata arms / test=logrank adjust=bon diff=all;
run;

/* 3. Survival estimates at fixed time points (12, 24, 36 months) */
proc lifetest data=proj.actg timelist=12 24 36 reduceout
              outsurv=proj.km_est notable;
  time time_mo*censor(1);
  strata arms;
run;

proc print data=proj.km_est noobs;
  title "KM survival probabilities at 12, 24, 36 months";
run;

ods graphics off;
ods listing close;
title;