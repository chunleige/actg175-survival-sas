/*-----------------------------------------------------------------
  Project : ACTG 175 survival analysis (portfolio project)
  Script  : 05_subgroup_forest.sas
  Purpose : Subgroup Cox models via a macro; forest plot of HRs
  Input   : proj.actg
  Output  : work.forest, forest_plot.png
  Author  : Ivyann
  Date    : 2026-11-13
  Notes   : - Comparison: any new regimen (treat=1) vs ZDV only (treat=0).
            - Subgroup analyses are exploratory; not powered for
              subgroup-specific inference.
------------------------------------------------------------------*/

libname proj "~/actg175";
options fmtsearch=(proj);

/* Check: treat should be 0 for ZDV only and 1 for the other arms */
proc freq data=proj.actg;
  tables treat*arms / norow nocol nopercent;
run;

/* Add CD4 subgroup: 0 = <350, 1 = >=350 */
data work.sg;
  set proj.actg;
  cd4grp = (cd40 >= 350);
run;

/* Start with an empty results dataset */
proc datasets lib=work nolist; delete forest; quit;

/* Macro: fit a Cox model in one subgroup, append the HR */
%macro sub(var=, val=, label=, order=);
  proc phreg data=work.sg(where=(&var=&val));
    model days*cens(0) = treat / ties=efron risklimits;
    ods output ParameterEstimates=_pe;
  run;

  data _pe;
    length subgroup $40;
    set _pe;
    subgroup = "&label";
    order    = &order;
    keep subgroup order HazardRatio HRLowerCL HRUpperCL ProbChiSq;
  run;

  proc append base=work.forest data=_pe force; run;
%mend sub;

/* Overall (all patients) */
%sub(var=treat, val=treat, label=All patients, order=1);

%sub(var=agegrp,  val=1, label=Age <35,             order=2);
%sub(var=agegrp,  val=2, label=Age >=35,            order=3);
%sub(var=gender,  val=0, label=Female,              order=4);
%sub(var=gender,  val=1, label=Male,                order=5);
%sub(var=race,    val=0, label=White,               order=6);
%sub(var=race,    val=1, label=Non-white,           order=7);
%sub(var=str2,    val=0, label=Antiretroviral-naive,order=8);
%sub(var=str2,    val=1, label=ART-experienced,     order=9);
%sub(var=symptom, val=0, label=Asymptomatic,        order=10);
%sub(var=symptom, val=1, label=Symptomatic,         order=11);
%sub(var=cd4grp,  val=0, label=Baseline CD4 <350,   order=12);
%sub(var=cd4grp,  val=1, label=Baseline CD4 >=350,  order=13);

proc sort data=work.forest; by order; run;

proc print data=work.forest noobs;
  title "Subgroup hazard ratios: new regimen vs ZDV only";
run;

/* Forest plot */
ods listing gpath="~/actg175" image_dpi=300;
ods graphics on / imagename="forest_plot" imagefmt=png width=7in height=6in;

proc sgplot data=work.forest noautolegend;
  scatter y=subgroup x=HazardRatio /
          xerrorlower=HRLowerCL xerrorupper=HRUpperCL
          markerattrs=(symbol=squarefilled size=9);
  refline 1 / axis=x lineattrs=(pattern=shortdash);
  xaxis type=log logbase=2 min=0.125 max=2
        label="Hazard ratio (new regimen vs ZDV only), log scale";
  yaxis discreteorder=data reverse display=(nolabel);
  title "Subgroup analysis: composite endpoint";
run;

ods graphics off;
ods listing close;
title;

/* Interaction tests: does the treatment effect differ across subgroups? */
proc datasets lib=work nolist; delete inter; quit;

%macro inter(var=);
  proc phreg data=work.sg;
    class &var(ref="0") / param=ref;
    model days*cens(0) = treat &var treat*&var / ties=efron;
    ods output Type3=_t3;
  run;
  data _t3;
    set _t3;
    if index(Effect, "*") > 0;   /* keep the interaction row only */
  run;
  proc append base=work.inter data=_t3 force; run;
%mend inter;

%inter(var=gender);
%inter(var=race);
%inter(var=str2);
%inter(var=symptom);
%inter(var=cd4grp);

/* agegrp is coded 1/2, so its reference is 1 */
options mprint;

proc phreg data=work.sg;
  class agegrp(ref="1") / param=ref;
  model days*cens(0) = treat agegrp treat*agegrp / ties=efron;
  ods output Type3=_t3;
run;
  data _t3;
    length Effect $20;
    set _t3;
    if index(Effect, "*") > 0;
  run;
proc append base=work.inter data=_t3 force; run;



options mprint;

/* Remove old results so the new column length takes effect */
proc datasets lib=work nolist; delete inter; quit;

%macro inter(var=, ref=0);
  proc phreg data=work.sg;
    class &var(ref="&ref") / param=ref;
    model days*cens(0) = treat &var treat*&var / ties=efron;
    ods output Type3=_t3;
  run;
  data _t3;
    length Effect $20;
    set _t3;
    if index(Effect, "*") > 0;
  run;
  proc append base=work.inter data=_t3 force; run;
%mend inter;

%inter(var=gender);
%inter(var=race);
%inter(var=str2);
%inter(var=symptom);
%inter(var=cd4grp);
%inter(var=agegrp, ref=1);   /* agegrp is coded 1/2 */

proc print data=work.inter noobs;
  title "Treatment-by-subgroup interaction tests";
run;
title;
proc print data=work.inter noobs;
  title "Treatment-by-subgroup interaction tests";
run;
title;