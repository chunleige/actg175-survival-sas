# -----------------------------------------------------------------
# Project : ACTG 175 survival analysis (portfolio project)
# Script  : 08_rmst.R
# Purpose : Cross-validate SAS RMST results with the survRM2 package
#           (each arm vs. ZDV only, tau = 36 months)
# Input   : the same raw UCI CSV used in SAS
# Output  : RMST by arm and RMST differences vs. ZDV only
# Author  : Ivyann
# Date    : 2026-11-13
# Notes   : - survRM2::rmst2() compares two groups at a time, so each
#             new regimen is compared with ZDV only separately.
#           - Point estimates should match SAS; standard errors may
#             differ slightly because variance methods differ.
#           - Run once first: install.packages("survRM2")
# -----------------------------------------------------------------

library(survRM2)

path <- "C:/Ge/Job/BC-CfE/ACTG175/actg175.csv"
d <- read.csv(path)
stopifnot(nrow(d) == 2139)

d$time_mo <- d$days / 30.4375
tau <- 36

arm_labels <- c("0" = "ZDV only", "1" = "ZDV + ddI",
                "2" = "ZDV + Zal", "3" = "ddI only")

res <- data.frame()

for (a in 1:3) {
  sub <- d[d$arms %in% c(0, a), ]
  arm_ind <- as.integer(sub$arms == a)       # 1 = new regimen, 0 = ZDV only

  fit <- rmst2(time = sub$time_mo, status = sub$cens,
               arm = arm_ind, tau = tau)

  est <- fit$unadjusted.result
  res <- rbind(res, data.frame(
    comparison   = paste(arm_labels[as.character(a)], "vs ZDV only"),
    rmst_new     = fit$RMST.arm1$rmst["Est."],
    rmst_zdv     = fit$RMST.arm0$rmst["Est."],
    diff         = est["RMST (arm=1)-(arm=0)", "Est."],
    lower_95     = est["RMST (arm=1)-(arm=0)", "lower .95"],
    upper_95     = est["RMST (arm=1)-(arm=0)", "upper .95"],
    p_unadjusted = est["RMST (arm=1)-(arm=0)", "p"]
  ))
}

res$p_bonferroni <- pmin(res$p_unadjusted * 3, 1)
rownames(res) <- NULL

options(width = 120)
print(res, digits = 4)

write.csv(res, "rmst_r_results.csv", row.names = FALSE)
