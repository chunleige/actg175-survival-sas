/*-----------------------------------------------------------------
  Project : ACTG 175 survival analysis (portfolio project)
  Script  : 06_cd4_change.sas
  Purpose : Compare change in CD4 count from baseline to week 20
            across treatment arms using ANCOVA (adjusted for baseline)
  Input   : proj.actg
  Output  : LS means by arm, differences vs. ZDV only, diagnostic plots
  Author  : Ivyann
  Date    : 2026-11-13
  Notes   : - Outcome: cd4_chg20 = cd420 - cd40 (created in 02).
            - ANCOVA with baseline adjustment; treatment effect is the
              same whether the outcome is the change or the week-20 value.
            - Complete-case analysis: patients missing cd420 are excluded.
------------------------------------------------------------------*/

libname proj "~/actg175";
options fmtsearch=(proj);

ods listing gpath="~/actg175" image_dpi=300;
ods graphics on / imagename="cd4" imagefmt=png;

/* 1. Missingness of week-20 CD4 by arm */
data work.cd4;
  set proj.actg;
  miss20 = missing(cd420);   /* 1 = week-20 CD4 missing */
run;

proc freq data=work.cd4;
  tables arms*miss20 / nocol nopercent chisq;
  title "Missing week-20 CD4 by treatment arm";
run;

/* 2. Descriptive statistics by arm */
proc means data=proj.actg n mean std median maxdec=1;
  class arms;
  var cd40 cd420 cd4_chg20;
  title "CD4 at baseline and week 20, by arm";
run;

/* 3. Check: homogeneity of slopes (arm-by-baseline interaction) */
proc glm data=proj.actg;
  class arms;
  model cd4_chg20 = arms cd40 arms*cd40;
  title "Homogeneity of slopes check";
run;
quit;

/* 4. ANCOVA: adjusted means and comparisons vs. ZDV only */
proc glm data=proj.actg plots=diagnostics order=internal;
  class arms;
  model cd4_chg20 = arms cd40 / solution clparm;
  lsmeans arms / cl pdiff=control adjust=dunnett;
  title "ANCOVA: CD4 change at week 20, adjusted for baseline CD4";
run;
quit;

ods graphics off;
ods listing close;
title;

/* 5. ANCOVA with arm-by-baseline interaction (slopes differ, p=0.0015)
      LS means evaluated at mean baseline CD4 (default) and at 250 / 450 */
proc glm data=proj.actg order=internal;
  class arms;
  model cd4_chg20 = arms cd40 arms*cd40 / solution;
  lsmeans arms / cl pdiff=control adjust=dunnett;
  lsmeans arms / at cd40=250 cl pdiff=control adjust=dunnett;
  lsmeans arms / at cd40=450 cl pdiff=control adjust=dunnett;
  title "ANCOVA with separate slopes: effects at different baseline CD4";
run;
quit;
title;

/* 6. Figure: CD4 change vs. baseline CD4, separate regression line per arm */
ods listing gpath="~/actg175" image_dpi=300;
ods graphics on / imagename="cd4_slopes" imagefmt=png width=7in height=5in;

proc sgplot data=proj.actg;
  reg x=cd40 y=cd4_chg20 / group=arms nomarkers clm clmtransparency=0.8;
  refline 0 / axis=y lineattrs=(pattern=shortdash);
  xaxis label="Baseline CD4 (cells/mm³)";
  yaxis label="Change in CD4 at week 20 (cells/mm³)";
  title "CD4 change by baseline CD4 and treatment arm";
run;

ods graphics off;
ods listing close;
title;