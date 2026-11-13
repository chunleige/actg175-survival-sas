/*-----------------------------------------------------------------
  Project : ACTG 175 survival analysis (portfolio project)
  Script  : 04_cox_phreg.sas
  Purpose : Cox proportional hazards models: unadjusted and adjusted
            treatment effects, Efron vs. Breslow ties, PH assessment
  Input   : proj.actg
  Output  : hazard ratio tables, PH diagnostic plots, proj.hr_adj
  Author  : Ivyann
  Date    : 2026-11-12
  Notes   : - Endpoint: composite (CD4 decline >=50%, AIDS, or death).
            - Reference arm: ZDV only. Time in days (days).
            - PHREG default ties method is Breslow; Efron is set
              explicitly here and compared against Breslow.
            - Event is cens=1, so censoring value is 0 in the TIME
              statement (do NOT use the 'censor' variable here).
------------------------------------------------------------------*/

libname proj "~/actg175";
options fmtsearch=(proj);

ods listing gpath="~/actg175" image_dpi=300;
ods graphics on / imagename="cox" imagefmt=png;

/* Model 1: unadjusted (Efron ties); reference arm 0 = ZDV only */
proc phreg data=proj.actg;
  class arms(ref="0") / param=ref;
  model days*cens(0) = arms / ties=efron risklimits;
  title "Model 1: Unadjusted Cox model (Efron ties)";
run;

/* Model 2: adjusted for baseline covariates (Efron ties) */
proc phreg data=proj.actg;
  class arms(ref="0") gender(ref="0") race(ref="0")
        hemo(ref="0") homo(ref="0") drugs(ref="0")
        symptom(ref="0") str2(ref="0") / param=ref;
  model days*cens(0) = arms age wtkg karnof cd40 cd80
                       gender race hemo homo drugs symptom str2
        / ties=efron risklimits;
  ods output ParameterEstimates=proj.hr_adj;
  title "Model 2: Adjusted Cox model (Efron ties)";
run;

/* Ties comparison: Breslow vs. Efron on Model 1 */
proc phreg data=proj.actg;
  class arms(ref="0") / param=ref;
  model days*cens(0) = arms / ties=breslow risklimits;
  title "Model 1 with Breslow ties (for comparison)";
run;

/* PH assessment (adjusted model) */
proc phreg data=proj.actg;
  class arms(ref="0") gender(ref="0") race(ref="0")
        hemo(ref="0") homo(ref="0") drugs(ref="0")
        symptom(ref="0") str2(ref="0") / param=ref;
  model days*cens(0) = arms age wtkg karnof cd40 cd80
                       gender race hemo homo drugs symptom str2
        / ties=efron;
  assess ph / resample seed=2026;
  title "PH assessment: cumulative Schoenfeld residuals";
run;

ods graphics off;
ods listing close;
title;

/* Sensitivity analysis: stratify by baseline CD4 quartile
   (cd40 violated the PH assumption, p<.0001) */
proc rank data=proj.actg out=work.actg_q groups=4;
  var cd40;
  ranks cd40_q;
run;

proc phreg data=work.actg_q;
  class arms(ref="0") gender(ref="0") race(ref="0")
        hemo(ref="0") homo(ref="0") drugs(ref="0")
        symptom(ref="0") str2(ref="0") / param=ref;
  model days*cens(0) = arms age wtkg karnof cd80
                       gender race hemo homo drugs symptom str2
        / ties=efron risklimits;
  strata cd40_q;
  title "Model 3: Stratified by baseline CD4 quartile";
run;
title;