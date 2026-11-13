/*-----------------------------------------------------------------
  Project : ACTG 175 survival analysis (portfolio project)
  Script  : 01_import.sas
  Purpose : Import the raw ACTG 175 CSV file and run initial checks
            (structure, first records, summary statistics, group sizes)
  Input   : ~/actg175/actg175.csv  (raw data, UCI ML Repository, id=890)
  Output  : proj.raw  (permanent SAS dataset, unmodified raw data)
  Author  : Ivyann
  Date    : 2026-11-11
  Notes   : - Raw data is saved unchanged; all derived variables are
              created in 02_prep_table1.sas.
            - Expected: 2,139 observations.
            - Check the log for ERROR/WARNING messages after each run.
            - Data source: AIDS Clinical Trials Group Study 175
              (Hammer et al., NEJM 1996). Follow the acknowledgment
              policy on the original dataset website when citing.
------------------------------------------------------------------*/

/* Permanent library: datasets stored here persist across sessions */
libname proj "~/actg175";

/* Import CSV; guessingrows=max scans all rows to determine variable types */
proc import datafile="~/actg175/actg175.csv"
    out=proj.raw
    dbms=csv
    replace;
    guessingrows=max;
run;

/* Check 1: structure of the dataset (variable names, types, order) */
proc contents data=proj.raw varnum;
run;

/* Check 2: first 10 records */
proc print data=proj.raw(obs=10);
run;

/* Check 3: summary of numeric variables, including missing counts */
proc means data=proj.raw n nmiss mean std min max;
run;

/* Check 4: treatment arm and event indicator counts */
proc freq data=proj.raw;
    tables arms cens / missing;
run;