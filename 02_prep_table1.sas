/*-----------------------------------------------------------------
  Project : ACTG 175 survival analysis (portfolio project)
  Script  : 02_prep_table1.sas
  Purpose : Derive analysis variables, apply formats/labels,
            and produce Table 1 (baseline characteristics by arm)
  Input   : proj.raw
  Output  : proj.actg (analysis dataset), table1.rtf
  Author  : Ivyann
  Date    : 2026-11-11
  Notes   : - cens = 1 indicates the composite endpoint (CD4 decline
              >=50%, progression to AIDS, or death), NOT death alone.
            - Variable coding below follows the dataset documentation;
              verify against the original data description.
            - cd496 is character in the raw file; converted to numeric.
------------------------------------------------------------------*/

libname proj "~/actg175";
options fmtsearch=(proj);

proc format library=proj;
  value armf   0="ZDV only" 1="ZDV + ddI" 2="ZDV + Zal" 3="ddI only";
  value sexf   0="Female" 1="Male";
  value racef  0="White" 1="Non-white";
  value yn     0="No" 1="Yes";
  value expf   0="Antiretroviral-naive" 1="Experienced";
  value stratf 1="Naive" 2="1-52 weeks prior ART" 3=">52 weeks prior ART";
  value evtf   0="Censored" 1="Event";
run;

data proj.actg;
  set proj.raw;
  time_mo   = days / 30.4375;                 /* days -> months */
  censor    = 1 - cens;                       /* 1 = censored (for LIFETEST/PHREG) */
  cd4_chg20 = cd420 - cd40;                   /* CD4 change, baseline to week 20 */
  cd496_n   = input(cd496, ?? best12.);       /* char -> numeric; non-numeric -> missing */
  if age < 35 then agegrp = 1; else agegrp = 2;

  label days      = "Time to event/censoring (days)"
        time_mo   = "Time (months)"
        cens      = "Event indicator (1 = event)"
        arms      = "Treatment arm"
        cd40      = "Baseline CD4"
        cd420     = "CD4 at week 20"
        cd496_n   = "CD4 at week 96"
        cd4_chg20 = "CD4 change, baseline to week 20"
        agegrp    = "Age group (1: <35, 2: >=35)";

  format arms armf. gender sexf. race racef. cens evtf.
         hemo homo drugs symptom yn. str2 expf. strat stratf.;
run;

/* Checks: missingness, event counts */
proc means data=proj.actg n nmiss min max;
  var days cd40 cd420 cd496_n cd4_chg20;
run;

proc freq data=proj.actg;
  tables arms*cens / nocol nopercent;
run;

/* Table 1 -> RTF */
ods rtf file="~/actg175/table1.rtf" style=journal;

proc tabulate data=proj.actg;
  class arms;
  var age wtkg karnof cd40 cd80;
  table (age wtkg karnof cd40 cd80) * (n*f=6.0 mean*f=7.1 std*f=7.1),
        arms="Treatment arm" all="Overall";
  title "Table 1a. Baseline continuous characteristics";
run;

proc tabulate data=proj.actg;
  class arms gender race hemo homo drugs symptom str2 strat / missing;
  table (gender race hemo homo drugs symptom str2 strat) * (n*f=6.0 colpctn*f=7.1),
        arms="Treatment arm" all="Overall";
  title "Table 1b. Baseline categorical characteristics";
run;

ods rtf close;
title;