# Code for "Assessing State Making Strategies" by D. Lemke and M. Karstens

# Required packages for replication - Please install if needed - VERIFY

library(survival)
library(stargazer)
library(dplyr)
library(ggplot2)
library(pscl)
library(censReg)
library(survminer)
library("xfun")
library(ggthemes)
library(extrafont)
font_import()
loadfonts(device = "win")
library(grid)
library(lintr) #For maintaining coding style
library(foreign)

library(sampleSelection) # Package containing heckit for heckprob substitute
library(lubridate) #Used for survival data formatting
library(reReg) #recurrent events package
library(gmodels)
library(Hmisc)
library(visreg)
library(coefplot)
library(broom)
library(GGally)

# Settings

options(scipen = 50) # bias against scientific notation for convenience

# Loading data and logging area

tc_dat <- read.csv("SMS-Replication-Data.csv")
names(tc_dat)

# Subset datasets

tc_multi <- tc_dat[tc_dat$event_num >= 2, ]

tc23 <- filter(tc_dat, event_num == 2 | event_num == 3)

tc45 <- filter(tc_dat, event_num == 4 | event_num == 5)

tc6_more <- filter(tc_dat, event_num >= 6)

tc1 <- filter(tc_dat, event_num == 1)

tc2 <- filter(tc_dat, event_num == 2)

tc3 <- filter(tc_dat, event_num == 3)

tc4 <- filter(tc_dat, event_num == 4)

tc5 <- filter(tc_dat, event_num == 5)

tc6 <- filter(tc_dat, event_num == 6)

# Survival objects

tc_dat_elapse <- Surv(tc_dat$alt_start, tc_dat$alt_stop, tc_dat$new_tc_dummy)
tc_dat_gap <- Surv(tc_dat$start, tc_dat$stop, tc_dat$new_tc_dummy)

tc_multi_elapse <- Surv(tc_multi$alt_start, tc_multi$alt_stop,
                        tc_multi$new_tc_dummy)
tc_multi_gap <- Surv(tc_multi$start, tc_multi$stop, tc_multi$new_tc_dummy)

tc1_elapse <- Surv(tc1$alt_start, tc1$alt_stop, tc1$new_tc_dummy)
tc1_gap <- Surv(tc1$start, tc1$stop, tc1$new_tc_dummy)

tc2_elapse <- Surv(tc2$alt_start, tc2$alt_stop, tc2$new_tc_dummy)
tc2_gap <- Surv(tc2$start, tc2$stop, tc2$new_tc_dummy)

tc23_elapse <- Surv(tc23$alt_start, tc23$alt_stop, tc23$new_tc_dummy)
tc23_gap <- Surv(tc23$start, tc23$stop, tc23$new_tc_dummy)

tc3_elapse <- Surv(tc3$alt_start, tc3$alt_stop, tc3$new_tc_dummy)
tc3_gap <- Surv(tc3$start, tc3$stop, tc3$new_tc_dummy)

tc4_elapse <- Surv(tc4$alt_start, tc4$alt_stop, tc4$new_tc_dummy)
tc4_gap <- Surv(tc4$start, tc4$stop, tc4$new_tc_dummy)

tc45_elapse <- Surv(tc45$alt_start, tc45$alt_stop, tc45$new_tc_dummy)
tc45_gap <- Surv(tc45$start, tc45$stop, tc45$new_tc_dummy)

tc5_elapse <- Surv(tc5$alt_start, tc5$alt_stop, tc5$new_tc_dummy)
tc5_gap <- Surv(tc5$start, tc5$stop, tc5$new_tc_dummy)

tc6_elapse <- Surv(tc6$alt_start, tc6$alt_stop, tc6$new_tc_dummy)
tc6_gap <- Surv(tc6$start, tc6$stop, tc6$new_tc_dummy)

tc6_more_elapse <- Surv(tc6_more$alt_start, tc6_more$alt_stop,
                          tc6_more$new_tc_dummy)
tc6_more_gap <- Surv(tc6_more$start, tc6_more$stop, tc6_more$new_tc_dummy)

#Linear Decay 20 PWP Models

tc_pwp <- coxph(tc_dat_gap ~ Peace_wtc + Fighting_wtc  + Forceful_LinDecay +
                GoodEnd_LinDecay + polity2 + area_1000_log + lmtnest + ELF +
                tc_tally + strata(event_num) + cluster(ccode), data = tc_dat,
                method = "efron")

tc1_pwp <- coxph(tc1_gap ~ polity2 + area_1000_log + lmtnest + ELF +
                 strata(event_num) + cluster(ccode), data = tc1,
                 method = "efron")

tc23_pwp_lin20 <- coxph(tc23_gap ~ Peace_wtc + Fighting_wtc  +
                        Forceful_LinDecay + GoodEnd_LinDecay + polity2 +
                        area_1000_log + lmtnest + ELF +
                        strata(event_num) + cluster(ccode), data = tc23,
                        method = "efron")

tc45_pwp_lin20 <- coxph(tc45_gap ~ Peace_wtc + Fighting_wtc +
                        Forceful_LinDecay + GoodEnd_LinDecay + polity2 +
                        area_1000_log + lmtnest + ELF +
                        strata(event_num) + cluster(ccode), data = tc45,
                        method = "efron")

tc6_more_pwp_lin20 <- coxph(tc6_more_gap ~ Peace_wtc + Fighting_wtc  +
                            Forceful_LinDecay + GoodEnd_LinDecay + polity2 +
                            area_1000_log + lmtnest + ELF +
                            strata(event_num) + cluster(ccode),
                            data = tc6_more, method = "efron")

summary(tc_pwp)


# Main Results Table

stargazer(tc_pwp, tc1_pwp, tc23_pwp_lin20, tc45_pwp_lin20, tc6_more_pwp_lin20,
          type = "latex",
          title = "PWP Gap Time Model Results",
          #dep.var.labels.include = F,
          model.numbers = F,
          column.labels = c("All TCs",
                             "First TC",
                             " TC 2 or 3",
                             "TC 4 or 5",
                             "TC 6+"),
          dep.var.labels = c("", "Separate PWP Models, By TC Iteration", "",
                             "", ""),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                             "Forceful Reintegration",
                             "Favorable Outcome for TC",
                             "Polity 2 Score",
                             "State Area (logged)",
                             "Mountains",
                             "ELF",
                             "TC Tally"),
          keep.stat = c("n"))

# Coefficient Plot of Main Results - In Progress


#Linear Decay 20 PWP Models without tally

tc_pwp_no_tally <- coxph(tc_dat_gap ~ Peace_wtc + Fighting_wtc  + Forceful_LinDecay +
                        GoodEnd_LinDecay + polity2 + area_1000_log + lmtnest + ELF +
                  strata(event_num) + cluster(ccode), data = tc_dat,
                method = "efron")


stargazer(tc_pwp_no_tally, tc1_pwp, tc23_pwp_lin20, tc45_pwp_lin20, tc6_more_pwp_lin20,
          type = "latex",
          title = "PWP Gap Time Model Results - No TC Tally",
          #dep.var.labels.include = F,
          model.numbers = F,
          column.labels = c("All TCs",
                            "First TC",
                            " TC 2 or 3",
                            "TC 4 or 5",
                            "TC 6+"),
          dep.var.labels = c("", "Separate PWP Models, By TC Iteration", "",
                             "", ""),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Forceful Reintegration",
                               "Favorable Outcome for TC",
                               "Polity 2 Score",
                               "State Area (logged)",
                               "Mountains",
                               "ELF"),
          keep.stat = c("n"))

# State Failure by other tc as a control

tc_failure <- coxph(tc_dat_gap ~ Peace_wtc + Fighting_wtc  +
                     Forceful_LinDecay + GoodEnd_LinDecay +
                     sf + polity2 + area_1000_log +
                     lmtnest + ELF + tc_tally + strata(event_num) +
                     cluster(ccode), data = tc_dat, method = "efron")

tc1_failure <- coxph(tc1_gap ~ sf + polity2 + area_1000_log + lmtnest + ELF +
                    strata(event_num) + cluster(ccode), data = tc1,
                    method = "efron")

tc23_failure <- coxph(tc23_gap ~ Peace_wtc + Fighting_wtc  + Forceful_LinDecay +
                       GoodEnd_LinDecay + sf + polity2 +
                       area_1000_log + lmtnest + ELF + strata(event_num) +
                       cluster(ccode), data = tc23, method = "efron")

tc45_failure <- coxph(tc45_gap ~ Peace_wtc + Fighting_wtc  + Forceful_LinDecay +
                       GoodEnd_LinDecay + sf + polity2 +
                       area_1000_log + lmtnest + ELF + strata(event_num) +
                       cluster(ccode), data = tc45, method = "efron")

tc6_more_failure <- coxph(tc6_more_gap ~ Peace_wtc + Fighting_wtc  +
                           Forceful_LinDecay + GoodEnd_LinDecay +
                           sf + polity2 + area_1000_log +
                           lmtnest + ELF + strata(event_num) +
                           cluster(ccode), data = tc6_more, method = "efron")

summary(tc_failure)

stargazer(tc_failure, tc1_failure, tc23_failure, tc45_failure,
          tc6_more_failure, type = "latex",
          title = "PWP Gap Time - State Failure as a Control",
          dep.var.labels = c("All TCs",
                             "First TC",
                             "TC 2 or 3",
                             "TC 4 or 5",
                             "TC 6+"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Forceful Reintegration",
                               "Favorable Outcome for TC",
                               "State Failure",
                               "Polity 2 Score",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))


# Absorbed by other tc as a control

tc_absorb <- coxph(tc_dat_gap ~ Peace_wtc + Fighting_wtc  +
                   Forceful_LinDecay + GoodEnd_LinDecay +
                   Absorbed_LinDecay + polity2 + area_1000_log +
                   lmtnest + ELF + tc_tally + strata(event_num) +
                   cluster(ccode), data = tc_dat, method = "efron")

tc1_absorb <- coxph(tc1_gap ~ polity2 + area_1000_log + lmtnest + ELF +
                    strata(event_num) + cluster(ccode), data = tc1,
                    method = "efron")

tc23_absorb <- coxph(tc23_gap ~ Peace_wtc + Fighting_wtc  + Forceful_LinDecay +
                     GoodEnd_LinDecay + Absorbed_LinDecay + polity2 +
                     area_1000_log + lmtnest + ELF + strata(event_num) +
                     cluster(ccode), data = tc23, method = "efron")

tc45_absorb <- coxph(tc45_gap ~ Peace_wtc + Fighting_wtc  + Forceful_LinDecay +
                     GoodEnd_LinDecay + Absorbed_LinDecay + polity2 +
                     area_1000_log + lmtnest + ELF + strata(event_num) +
                     cluster(ccode), data = tc45, method = "efron")

tc6_more_absorb <- coxph(tc6_more_gap ~ Peace_wtc + Fighting_wtc  +
                         Forceful_LinDecay + GoodEnd_LinDecay +
                         Absorbed_LinDecay + polity2 + area_1000_log +
                         lmtnest + ELF + strata(event_num) +
                         cluster(ccode), data = tc6_more, method = "efron")

summary(tc_absorb)

stargazer(tc_absorb, tc1_absorb, tc23_absorb, tc45_absorb, tc6_more_absorb,
          type = "latex",
          title = "PWP Gap Time - Absorbed as a Control",
          dep.var.labels = c("All TCs",
                             "First TC",
                             "TC 2 or 3",
                             "TC 4 or 5",
                             "TC 6+"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Forceful Reintegration",
                               "Favorable Outcome for TC",
                               "Absorbed",
                               "Polity 2 Score",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))

# All End Types Disaggregated

tc_disag <- coxph(tc_dat_gap ~ Peace_wtc + Fighting_wtc  + Forceful_LinDecay +
                  Promoted_LinDecay + Peaceful_LinDecay + Absorbed_LinDecay +
                  polity2 + area_1000_log + lmtnest + ELF + tc_tally +
                  strata(event_num) + cluster(ccode), data = tc_dat,
                  method = "efron")

tc1_disag <- coxph(tc1_gap ~ +polity2 + area_1000_log + lmtnest + ELF +
                   strata(event_num) + cluster(ccode), data = tc1,
                   method = "efron")

tc23_disag <- coxph(tc23_gap ~ Peace_wtc + Fighting_wtc  + Forceful_LinDecay +
                    Promoted_LinDecay + Peaceful_LinDecay +
                    Absorbed_LinDecay + polity2 + area_1000_log + lmtnest +
                    ELF + strata(event_num) + cluster(ccode), data = tc23,
                    method = "efron")

tc45_disag <- coxph(tc45_gap ~ Peace_wtc + Fighting_wtc  + Forceful_LinDecay +
                    Promoted_LinDecay + Peaceful_LinDecay + Absorbed_LinDecay +
                    polity2 + area_1000_log + lmtnest + ELF +
                    strata(event_num) + cluster(ccode), data = tc45,
                    method = "efron")

tc6_more_disag <- coxph(tc6_more_gap ~ Peace_wtc + Fighting_wtc  +
                        Forceful_LinDecay + Promoted_LinDecay +
                        Peaceful_LinDecay + Absorbed_LinDecay + polity2 +
                        area_1000_log + lmtnest + ELF + strata(event_num) +
                        cluster(ccode), data = tc6_more, method = "efron")

summary(tc_disag)

stargazer(tc_disag, tc1_disag, tc23_disag, tc45_disag, tc6_more_disag,
          type = "latex",
          title = "PWP Gap Time - Disaggregated",
          dep.var.labels = c("All TCs",
                             "First TC",
                             "TC 2 or 3",
                             "TC 4 or 5",
                             "tc 6+"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Forceful Reintegration",
                               "Promoted",
                               "Peaceful Reintegration",
                               "Absorbed",
                               "Polity 2 Score",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))

# Exponential Decay 20 year half life

tc_exp20 <- coxph(tc_dat_gap ~ Peace_wtc + Fighting_wtc  +
                  Forceful_ExpDecay_20 + GoodEnd_ExpDecay_20 + polity2 +
                  area_1000_log + lmtnest + ELF + tc_tally +
                  strata(event_num) + cluster(ccode), data = tc_dat,
                  method = "efron")

tc1_exp20 <- coxph(tc1_gap ~ polity2 + area_1000_log + lmtnest + ELF +
                   strata(event_num) + cluster(ccode), data = tc1,
                   method = "efron")

tc23_exp20 <- coxph(tc23_gap ~ Peace_wtc + Fighting_wtc  +
                    Forceful_ExpDecay_20 + GoodEnd_ExpDecay_20 + polity2 +
                    area_1000_log + lmtnest + ELF + strata(event_num) +
                    cluster(ccode), data = tc23, method = "efron")

tc45_exp20 <- coxph(tc45_gap ~ Peace_wtc + Fighting_wtc  +
                    Forceful_ExpDecay_20 + GoodEnd_ExpDecay_20 + polity2 +
                    area_1000_log + lmtnest + ELF + strata(event_num) +
                    cluster(ccode), data = tc45, method = "efron")

tc6_more_exp20 <- coxph(tc6_more_gap ~ Peace_wtc + Fighting_wtc  +
                        Forceful_ExpDecay_20 + GoodEnd_ExpDecay_20 + polity2 +
                        area_1000_log + lmtnest + ELF + strata(event_num) +
                        cluster(ccode), data = tc6_more, method = "efron")

summary(tc_exp20)


stargazer(tc_exp20, tc1_exp20, tc23_exp20, tc45_exp20, tc6_more_exp20,
          type = "latex",
          title = "PWP Gap Time - Exponential Decay with 20 year Half Life",
          dep.var.labels = c("All tcs",
                             "First TC",
                             "TC 2 or 3",
                             "TC 4 or 5",
                             "TC 6+"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Forceful Reintegration",
                               "Favorable Outcome for TC",
                               "Polity 2 Score",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))

# Exponential Decay 15 year half life

tc_exp15 <- coxph(tc_dat_gap ~ Peace_wtc + Fighting_wtc  +
                  Forceful_ExpDecay_15 + GoodEnd_ExpDecay_15 + polity2 +
                  area_1000_log + lmtnest + ELF + tc_tally +
                  strata(event_num) + cluster(ccode), data = tc_dat,
                  method = "efron")

tc1_exp15 <- coxph(tc1_gap ~ polity2 + area_1000_log + lmtnest + ELF +
                   strata(event_num) + cluster(ccode), data = tc1,
                   method = "efron")

tc23_exp15 <- coxph(tc23_gap ~ Peace_wtc + Fighting_wtc  +
                      Forceful_ExpDecay_15 + GoodEnd_ExpDecay_15 + polity2 +
                      area_1000_log + lmtnest + ELF + strata(event_num) +
                      cluster(ccode), data = tc23, method = "efron")

tc45_exp15 <- coxph(tc45_gap ~ Peace_wtc + Fighting_wtc  +
                      Forceful_ExpDecay_15 + GoodEnd_ExpDecay_15 + polity2 +
                      area_1000_log + lmtnest + ELF + strata(event_num) +
                      cluster(ccode), data = tc45, method = "efron")

tc6_more_exp15 <- coxph(tc6_more_gap ~ Peace_wtc + Fighting_wtc  +
                          Forceful_ExpDecay_15 + GoodEnd_ExpDecay_15 + polity2 +
                          area_1000_log + lmtnest + ELF + strata(event_num) +
                          cluster(ccode), data = tc6_more, method = "efron")


stargazer(tc_exp15, tc1_exp15, tc23_exp15, tc45_exp15, tc6_more_exp15,
          type = "latex",
          title = "PWP Gap Time - Exponential Decay with 15 year Half Life",
          dep.var.labels = c("All tcs",
                             "First TC",
                             "TC 2 or 3",
                             "TC 4 or 5",
                             "TC 6+"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Forceful Reintegration",
                               "Favorable Outcome for TC",
                               "Polity 2 Score",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))

# Exponential Decay 10 year half life

tc_exp10 <- coxph(tc_dat_gap ~ Peace_wtc + Fighting_wtc  +
                  Forceful_ExpDecay_10 + GoodEnd_ExpDecay_10 + polity2 +
                  area_1000_log + lmtnest + ELF + tc_tally +
                  strata(event_num) + cluster(ccode), data = tc_dat,
                  method = "efron")

tc1_exp10 <- coxph(tc1_gap ~ polity2 + area_1000_log + lmtnest + ELF +
                   strata(event_num) + cluster(ccode), data = tc1,
                   method = "efron")

tc23_exp10 <- coxph(tc23_gap ~ Peace_wtc + Fighting_wtc  +
                    Forceful_ExpDecay_10 + GoodEnd_ExpDecay_10 + polity2 +
                    area_1000_log + lmtnest + ELF + strata(event_num) +
                    cluster(ccode), data = tc23, method = "efron")

tc45_exp10 <- coxph(tc45_gap ~ Peace_wtc + Fighting_wtc  +
                    Forceful_ExpDecay_10 + GoodEnd_ExpDecay_10 + polity2 +
                    area_1000_log + lmtnest + ELF + strata(event_num) +
                    cluster(ccode), data = tc45, method = "efron")

tc6_more_exp10 <- coxph(tc6_more_gap ~ Peace_wtc + Fighting_wtc  +
                        Forceful_ExpDecay_10 + GoodEnd_ExpDecay_10 + polity2 +
                        area_1000_log + lmtnest + ELF + strata(event_num) +
                        cluster(ccode), data = tc6_more, method = "efron")

summary(tc_exp10)


stargazer(tc_exp10, tc1_exp10, tc23_exp10, tc45_exp10, tc6_more_exp10,
          type = "latex",
          title = "PWP Gap Time - Exponential Decay with 10 Year Half Life",
          dep.var.labels = c("All TCs",
                             "First TC",
                             "TC 2 or 3",
                             "TC 4 or 5",
                             "TC 6+"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Forceful Reintegration",
                               "Favorable Outcome for TC",
                               "Polity 2 Score",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))

# Exponential Decay 5 year half life

tc_exp5 <- coxph(tc_dat_gap ~ Peace_wtc + Fighting_wtc  +
                 Forceful_ExpDecay_5 + GoodEnd_ExpDecay_5 + polity2 +
                 area_1000_log + lmtnest + ELF + tc_tally +
                 strata(event_num) + cluster(ccode), data = tc_dat,
                  method = "efron")

tc1_exp5 <- coxph(tc1_gap ~ polity2 + area_1000_log + lmtnest + ELF +
                  strata(event_num) + cluster(ccode), data = tc1,
                  method = "efron")

tc23_exp5 <- coxph(tc23_gap ~ Peace_wtc + Fighting_wtc  +
                   Forceful_ExpDecay_5 + GoodEnd_ExpDecay_5 + polity2 +
                   area_1000_log + lmtnest + ELF + strata(event_num) +
                   cluster(ccode), data = tc23, method = "efron")

tc45_exp5 <- coxph(tc45_gap ~ Peace_wtc + Fighting_wtc  +
                   Forceful_ExpDecay_5 + GoodEnd_ExpDecay_5 + polity2 +
                   area_1000_log + lmtnest + ELF + strata(event_num) +
                   cluster(ccode), data = tc45, method = "efron")

tc6_more_exp5 <- coxph(tc6_more_gap ~ Peace_wtc + Fighting_wtc  +
                       Forceful_ExpDecay_5 + GoodEnd_ExpDecay_5 + polity2 +
                       area_1000_log + lmtnest + ELF + strata(event_num) +
                       cluster(ccode), data = tc6_more, method = "efron")


stargazer(tc_exp5, tc1_exp5, tc23_exp5, tc45_exp5, tc6_more_exp5,
          type = "latex",
          title = "PWP Gap Time - Exponential Decay with 5 year Half Life",
          dep.var.labels = c("All TCs",
                             "First TC",
                             "TC 2 or 3",
                             "TC 4 or 5",
                             "TC 6+"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Forceful Reintegration",
                               "Favorable Outcome for TC",
                               "Polity 2 Score",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))

# 10 Year Binary

tc_binary10 <- coxph(tc_dat_gap ~ Peace_wtc + Fighting_wtc  + Forceful_t10 +
                       GoodEnd_t10 + polity2 + area_1000_log + lmtnest + ELF +
                       tc_tally + strata(event_num) + cluster(ccode),
                     data = tc_dat, method = "efron")

tc1_binary10 <- coxph(tc1_gap ~ polity2 + area_1000_log + lmtnest + ELF +
                        strata(event_num) + cluster(ccode), data = tc1,
                      method = "efron")

tc23_binary10 <- coxph(tc23_gap ~ Peace_wtc + Fighting_wtc  + Forceful_t10 +
                         GoodEnd_t10 + polity2 + area_1000_log + lmtnest + ELF +
                         strata(event_num) + cluster(ccode), data = tc23,
                       method = "efron")

tc45_binary10 <- coxph(tc45_gap ~ Peace_wtc + Fighting_wtc  + Forceful_t10 +
                         GoodEnd_t10 + polity2 + area_1000_log + lmtnest + ELF +
                         strata(event_num) + cluster(ccode), data = tc45,
                       method = "efron")

tc6_more_binary10 <- coxph(tc6_more_gap ~ Peace_wtc + Fighting_wtc +
                             Forceful_t10 + GoodEnd_t10 + polity2 +
                             area_1000_log + lmtnest + ELF + strata(event_num) +
                             cluster(ccode), data = tc6_more, method = "efron")


stargazer(tc_binary10, tc1_binary10, tc23_binary10, tc45_binary10,
          tc6_more_binary10, type = "latex",
          title = "PWP Gap Time - 10 Year Binary Indicator",
          dep.var.labels = c("All TCs",
                             "First TC",
                             " TC 2 or 3",
                             "TC 4 or 5",
                             "TC 6+"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Forceful Reintegration",
                               "Favorable Outcome for TC",
                               "Polity 2 Score",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))

# 20 Year Binary

tc_binary20 <- coxph(tc_dat_gap ~ Peace_wtc + Fighting_wtc  + Forceful_t20 +
                     GoodEnd_t20 + polity2 + area_1000_log + lmtnest + ELF +
                     tc_tally + strata(event_num) + cluster(ccode),
                     data = tc_dat, method = "efron")

tc1_binary20 <- coxph(tc1_gap ~ polity2 + area_1000_log + lmtnest + ELF +
                      strata(event_num) + cluster(ccode), data = tc1,
                      method = "efron")

tc23_binary20 <- coxph(tc23_gap ~ Peace_wtc + Fighting_wtc  + Forceful_t20 +
                       GoodEnd_t20 + polity2 + area_1000_log + lmtnest + ELF +
                       strata(event_num) + cluster(ccode), data = tc23,
                       method = "efron")

tc45_binary20 <- coxph(tc45_gap ~ Peace_wtc + Fighting_wtc  + Forceful_t20 +
                       GoodEnd_t20 + polity2 + area_1000_log + lmtnest + ELF +
                       strata(event_num) + cluster(ccode), data = tc45,
                       method = "efron")

tc6_more_binary20 <- coxph(tc6_more_gap ~ Peace_wtc + Fighting_wtc +
                           Forceful_t20 + GoodEnd_t20 + polity2 +
                           area_1000_log + lmtnest + ELF + strata(event_num) +
                           cluster(ccode), data = tc6_more, method = "efron")

summary(tc_binary20)


stargazer(tc_binary20, tc1_binary20, tc23_binary20, tc45_binary20,
          tc6_more_binary20, type = "latex",
          title = "PWP Gap Time - 20 Year Binary Indicator",
          dep.var.labels = c("All TCs",
                             "First TC",
                             " TC 2 or 3",
                             "TC 4 or 5",
                             "TC 6+"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Forceful Reintegration",
                               "Favorable Outcome for TC",
                               "Polity 2 Score",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))

#Linear Decay 20 PWP  Elapsed Models

tc_pwp_elapse <- coxph(tc_dat_elapse ~ Peace_wtc + Fighting_wtc  +
                       Forceful_LinDecay + GoodEnd_LinDecay + polity2 +
                       area_1000_log + lmtnest + ELF + tc_tally +
                       strata(event_num) + cluster(ccode), data = tc_dat,
                       method = "efron")

tc1_pwp_elapse <- coxph(tc1_elapse ~ polity2 + area_1000_log + lmtnest + ELF +
                        strata(event_num) + cluster(ccode), data = tc1,
                        method = "efron")

tc23_pwp_lin20_elapse <- coxph(tc23_elapse ~ Peace_wtc + Fighting_wtc  +
                               Forceful_LinDecay + GoodEnd_LinDecay + polity2 +
                               area_1000_log + lmtnest + ELF +
                               strata(event_num) + cluster(ccode), data = tc23,
                               method = "efron")

tc45_pwp_lin20_elapse <- coxph(tc45_elapse ~ Peace_wtc + Fighting_wtc +
                               Forceful_LinDecay + GoodEnd_LinDecay + polity2 +
                               area_1000_log + lmtnest + ELF +
                               strata(event_num) + cluster(ccode), data = tc45,
                               method = "efron")

tc6_more_pwp_lin20_elapse <- coxph(tc6_more_elapse ~ Peace_wtc + Fighting_wtc  +
                                   Forceful_LinDecay + GoodEnd_LinDecay + polity2 +
                                   area_1000_log + lmtnest + ELF +
                                   strata(event_num) + cluster(ccode),
                                   data = tc6_more, method = "efron")

summary(tc_pwp_elapse)


stargazer(tc_pwp_elapse, tc1_pwp_elapse, tc23_pwp_lin20_elapse,
          tc45_pwp_lin20_elapse, tc6_more_pwp_lin20_elapse,
          type = "latex",
          title = "PWP Elapsed Time",
          #dep.var.labels.include = F,
          model.numbers = F,
          column.labels = c("All TCs",
                            "First TC",
                            " TC 2 or 3",
                            "TC 4 or 5",
                            "TC 6+"),
          dep.var.labels = c("", "Separate PWP Models, By TC Iteration", "",
                             "", ""),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Forceful Reintegration",
                               "Favorable Outcome for TC",
                               "Polity 2 Score",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))


# Anderson-Gill Model - AG

tc_ag <- coxph(tc_dat_gap ~ Peace_wtc + Fighting_wtc +
               Forceful_LinDecay + GoodEnd_LinDecay + polity2 +
               area_1000_log + lmtnest + ELF +tc_tally +
               cluster(ccode), data = tc_dat, method="efron")

summary(tc_ag)

tc_ag_multi <- coxph(tc_multi_gap ~ Peace_wtc + Fighting_wtc +
                     Forceful_LinDecay + GoodEnd_LinDecay + polity2 +
                     area_1000_log + lmtnest + ELF +tc_tally +
                     cluster(ccode), data = tc_multi, method="efron")

summary(tc_ag_multi)

stargazer(tc_ag, tc_ag_multi,
          type = "latex",
          title = "Anderson Gill Variance Correction Model",
          #dep.var.labels.include = F,
          model.numbers = F,
          column.labels = c("All TCs",
                            "TC 2+"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Forceful Reintegration",
                               "Favorable Outcome for TC",
                               "Polity 2 Score",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))

# WLW Model

library(plyr)
tc_dat_expand <- tc_dat[rep(1:nrow(tc_dat), each =
                        max(tc_dat$event_num)),]

tc_dat_expand$one <- 1

tc_dat_expand <- ddply(tc_dat_expand, c("ccode", "year"), mutate,
                       eventrisk = cumsum(one))

tc_dat_expand$new_tc_dummy <- ifelse(tc_dat_expand$event_num == 
                                     tc_dat_expand$eventrisk &
                                     tc_dat_expand$new_tc_dummy
                                     == 1, 1, 0)

tc_dat_expand_s <- Surv(tc_dat_expand$alt_start,
                        tc_dat_expand$alt_stop,
                        tc_dat_expand$new_tc_dummy)

tc_dat_WLW <- coxph(tc_dat_expand_s ~ Peace_wtc + Fighting_wtc  +
                    Forceful_LinDecay + GoodEnd_LinDecay +
                    polity2 + area_1000_log + lmtnest + ELF +
                    tc_tally + strata(event_num) + cluster(ccode),
                    data = tc_dat_expand, method = "efron")

summary(tc_dat_WLW)

library(dplyr)

stargazer(tc_dat_WLW,
          type = "latex",
          title = "WLW Marginal Elapsed Time Model",
          #dep.var.labels.include = F,
          model.numbers = F,
          column.labels = c("Elapsed Time"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Forceful Reintegration",
                               "Favorable Outcome for TC",
                               "Polity 2 Score",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))

# Hurdle Model - Basic Poisson-Logit

hurdlemod <- hurdle(tc_tally ~ Peace_wtc + Fighting_wtc  +
                    Forceful_LinDecay + GoodEnd_LinDecay + polity2 +
                    area_1000_log + lmtnest + ELF | polity2 +
                    area_1000_log + lmtnest + ELF + Country_Age,
                    data = tc_dat, dist = c("poisson"),
                    zero.dist = c("poisson"), link = c("logit"))

summary(hurdlemod)

# ZINB

tc_zinb <- zeroinfl(tc_tally ~ Peace_wtc + Fighting_wtc  +
                    Forceful_LinDecay + GoodEnd_LinDecay + polity2 +
                    area_1000_log + lmtnest + ELF | polity2 +
                    area_1000_log + lmtnest + ELF, data = tc_dat,
                    dist = "negbin", link = "logit")

summary(tc_zinb)
