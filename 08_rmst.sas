/*-----------------------------------------------------------------
  Project : ACTG 175 survival analysis (portfolio project)
  Script  : 08_rmst.sas
  Purpose : Restricted mean survival time (RMST) by treatment arm,
            and RMST differences vs. ZDV only -- a treatment-effect
            summary that does not rely on proportional hazards
  Input   : proj.actg
  Output  : RMST by arm, differences vs. ZDV only, RMST plot
  Author  : Ivyann
  Date    : 2026-11-13
  Notes   : - tau = 36 months (event-free time restricted to 3 years);
              chosen so that every arm still has patients at risk.
            - RMST difference = average additional event-free months
              gained within the first 36 months.
            - R's cox.zph suggested mild non-PH for treatment (p=0.03),
              which motivates this analysis.
------------------------------------------------------------------*/

libname proj "~/actg175";
options fmtsearch=(proj);

ods listing gpath="~/actg175" image_dpi=300;
ods graphics on / imagename="rmst" imagefmt=png;

/* Check risk sets near tau: how many still at risk at 36 months? */
proc lifetest data=proj.actg timelist=30 33 36 reduceout
              outsurv=work.atrisk notable;
  time time_mo*censor(1);
  strata arms;
run;

proc print data=work.atrisk noobs;
  title "Survival and number at risk at 30, 33, 36 months";
run;

/* RMST by arm and differences vs. ZDV only */
proc lifetest data=proj.actg rmst(tau=36) plots=(rmst) notable;
  time time_mo*censor(1);
  strata arms;
  title "RMST (tau = 36 months) by treatment arm";
run;

ods graphics off;
ods listing close;
title;