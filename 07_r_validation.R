# -----------------------------------------------------------------
# Project : ACTG 175 survival analysis (portfolio project)
# Script  : 07_r_validation.R
# Purpose : Cross-validate the SAS results in R (survival package):
#           log-rank test, Cox models (Efron vs. Breslow ties),
#           adjusted and CD4-stratified models, PH check, ANCOVA
# Input   : the same raw UCI CSV used in SAS (01_import.sas)
# Output  : side-by-side SAS vs. R comparison table
# Author  : Ivyann
# Date    : 2026-11-13
# Notes   : - R's coxph() defaults to Efron ties; SAS PROC PHREG
#             defaults to Breslow. Both are set explicitly here.
#           - SAS values below are copied from the SAS output.
#           - cox.zph() tests PH with scaled Schoenfeld residuals;
#             SAS ASSESS PH uses cumulative martingale residuals,
#             so p-values are not expected to match exactly.
# -----------------------------------------------------------------

library(survival)

# ---- 1. Read data ------------------------------------------------
# Change this path to where your CSV is saved
path <- "C:/Ge/Job/BC-CfE/ACTG175/actg175.csv"
d <- read.csv(path)

stopifnot(nrow(d) == 2139)

d$arms <- factor(d$arms, levels = 0:3,
                 labels = c("ZDV only", "ZDV + ddI", "ZDV + Zal", "ddI only"))
d$cd4_chg20 <- d$cd420 - d$cd40

# Same quartile grouping as SAS PROC RANK GROUPS=4 (TIES=MEAN):
# group = floor(rank * 4 / (n + 1))
n <- nrow(d)
d$cd40_q <- floor(rank(d$cd40, ties.method = "average") * 4 / (n + 1))

# ---- 2. Log-rank test --------------------------------------------
lr <- survdiff(Surv(days, cens) ~ arms, data = d)
lr_chisq <- unname(lr$chisq)

# ---- 3. Cox models -----------------------------------------------
m1_efron   <- coxph(Surv(days, cens) ~ arms, data = d, ties = "efron")
m1_breslow <- coxph(Surv(days, cens) ~ arms, data = d, ties = "breslow")

covars <- "age + wtkg + karnof + cd40 + cd80 + gender + race +
           hemo + homo + drugs + symptom + str2"

m2 <- coxph(as.formula(paste("Surv(days, cens) ~ arms +", covars)),
            data = d, ties = "efron")

m3 <- coxph(Surv(days, cens) ~ arms + age + wtkg + karnof + cd80 +
              gender + race + hemo + homo + drugs + symptom + str2 +
              strata(cd40_q),
            data = d, ties = "efron")

# ---- 4. PH check (Schoenfeld residuals) -----------------------------
zph <- cox.zph(m2)
print(zph)

# ---- 5. ANCOVA: CD4 change at week 20 ------------------------------
# Without the interaction, the arm coefficients equal the
# LS-mean differences vs. ZDV only reported by SAS PROC GLM
anc <- lm(cd4_chg20 ~ arms + cd40, data = d)

# ---- 6. SAS vs. R comparison table -----------------------------------
arm_rows <- c("armsZDV + ddI", "armsZDV + Zal", "armsddI only")

r_coef <- function(fit) unname(coef(fit)[arm_rows])

sas <- list(
  m1_efron   = c(-0.70937, -0.64462, -0.53478),
  m1_breslow = c(-0.70912, -0.64439, -0.53464),
  m2         = c(-0.79223, -0.66523, -0.56784),
  m3         = c(-0.80974, -0.69228, -0.58953),
  ancova     = c(70.236910, 36.206562, 42.287715)
)

cmp <- data.frame(
  analysis = rep(c("Model 1 (Efron)", "Model 1 (Breslow)",
                   "Model 2 adjusted", "Model 3 CD4-stratified",
                   "ANCOVA CD4 change"), each = 3),
  contrast = rep(c("ZDV+ddI vs ZDV", "ZDV+Zal vs ZDV", "ddI vs ZDV"), 5),
  sas = c(sas$m1_efron, sas$m1_breslow, sas$m2, sas$m3, sas$ancova),
  r   = c(r_coef(m1_efron), r_coef(m1_breslow), r_coef(m2),
          r_coef(m3), r_coef(anc))
)
cmp$diff <- cmp$r - cmp$sas

# Hazard ratios for the Cox rows (log-HR -> HR)
cox_rows <- 1:12
cmp$sas_HR <- NA
cmp$r_HR   <- NA
cmp$sas_HR[cox_rows] <- exp(cmp$sas[cox_rows])
cmp$r_HR[cox_rows]   <- exp(cmp$r[cox_rows])

options(width = 120)
print(cmp, digits = 5, row.names = FALSE)

cat("\nLog-rank chi-square:  SAS = 49.1941   R =",
    round(lr_chisq, 4), "\n")

# Efron vs. Breslow: how much does the ties method matter here?
cat("\nNumber of distinct event times:",
    length(unique(d$days[d$cens == 1])),
    "  Total events:", sum(d$cens), "\n")

# Save the comparison table for the website / README
write.csv(cmp, "sas_vs_r_comparison.csv", row.names = FALSE)
