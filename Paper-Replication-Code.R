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

# Loading data

tc_dat <- read.csv("Survival-Data-TC.csv")
names(tc_dat)

summary(tc_dat$area_1000)
tc_dat$area_1000_log <- log(tc_dat$area_1000)
summary(tc_dat$area_1000_log)

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

tc_multi_elapse <- Surv(tc_dat$alt_start, tc_dat$alt_stop, tc_dat$new_tc_dummy)
tc_dat_gap <- Surv(tc_dat$start, tc_dat$stop, tc_dat$new_tc_dummy)

tc_multi_elapse <- Surv(tc_multi$alt_start, tc_multi$alt_stop, tc_multi$new_tc_dummy)
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

#Linear Decay PWP Models

tcAll.PWP <- coxph(tc_dat_gap ~ Peace_wtc + Fighting_wtc  + Forceful_LinDecay +
                  GoodEnd_LinDecay + polity2 + area_1000_log + lmtnest + ELF +
                  strata(event_num) + cluster(ccode), data = tc_dat,
                  method = "efron")

tc.1.PWP <- coxph(tc1_gap ~ polity2 + area_1000_log + lmtnest + ELF
                  + strata(event_num) + cluster(ccode), data = tc1,
                  method = "efron")

tc23.PWP.LinDecay <- coxph(tc23_gap ~ Peace_wtc + Fighting_wtc  +
                            Forceful_LinDecay + GoodEnd_LinDecay + polity2
                            + area_1000_log + lmtnest + ELF + strata(event_num)
                            + cluster(ccode), data = tc23, method = "efron")

tc45.PWP.LinDecay <- coxph(tc45_gap ~ Peace_wtc + Fighting_wtc +
                            Forceful_LinDecay + GoodEnd_LinDecay + polity2
                            + area_1000_log + lmtnest + ELF + strata(event_num)
                            + cluster(ccode), data = tc45, method = "efron")

tc6_more.PWP.LinDecay <- coxph(tc6_more_gap ~ Peace_wtc + Fighting_wtc  +
                               Forceful_LinDecay + GoodEnd_LinDecay + polity2
                               + area_1000_log + lmtnest + ELF
                               + strata(event_num) + cluster(ccode),
                               data = tc6_more, method = "efron")

summary(tcAll.PWP)


# Main Results Table

stargazer(tcAll.PWP, tc.1.PWP, tc23.PWP.LinDecay, tc45.PWP.LinDecay,
          tc6_more.PWP.LinDecay, type = "latex",
          title = "PWP Gap Time Model Results",
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
                             "ELF"),
          keep.stat = c("n"))


# Absorbed by other tc as a control

tcAll.AbsorbControl <- coxph(tc_dat_gap ~ Peace_wtc + Fighting_wtc  +
                               Forceful_LinDecay + GoodEnd_LinDecay +
                               Absorbed_LinDecay + polity2 + area_1000_log +
                               lmtnest + ELF + strata(event_num) +
                               cluster(ccode), data = tc_dat, method = "efron")

tc.1.AbsorbControl <- coxph(tc1_gap ~ polity2 + area_1000_log + lmtnest + ELF +
                            strata(event_num) + cluster(ccode), data = tc1,
                            method = "efron")

tc23.AbsorbControl <- coxph(tc23_gap ~ Peace_wtc + Fighting_wtc  +
                             Forceful_LinDecay + GoodEnd_LinDecay +
                             Absorbed_LinDecay + polity2 + area_1000_log +
                             lmtnest + ELF + strata(event_num) +
                             cluster(ccode), data = tc23, method = "efron")

tc45.AbsorbControl <- coxph(tc45_gap ~ Peace_wtc + Fighting_wtc  +
                             Forceful_LinDecay + GoodEnd_LinDecay +
                             Absorbed_LinDecay + polity2 + area_1000_log +
                             lmtnest + ELF + strata(event_num) +
                             cluster(ccode), data = tc45, method = "efron")

tc6_more.AbsorbControl <- coxph(tc6_more_gap ~ Peace_wtc + Fighting_wtc  +
                                Forceful_LinDecay + GoodEnd_LinDecay +
                                Absorbed_LinDecay + polity2 + area_1000_log +
                                lmtnest + ELF + strata(event_num) +
                                cluster(ccode), data = tc6_more,
                                method = "efron")

summary(tcAll.AbsorbControl)

stargazer(tcAll.AbsorbControl, tc.1.AbsorbControl, tc23.AbsorbControl,
          tc45.AbsorbControl, tc6_more.AbsorbControl,
          type = "latex",
          title = "PWP Gap Time Model Results - Absorbed as a Control",
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
                               "ELF"),
          keep.stat = c("n"))

# All End Types Disaggregated

tcAll.Disag <- coxph(tc_dat_gap ~ Peace_wtc + Fighting_wtc  + Forceful_LinDecay +
                    Promoted_LinDecay + Peaceful_LinDecay + Absorbed_LinDecay +
                    polity2 + area_1000_log + lmtnest + ELF + strata(event_num)
                    + cluster(ccode), data = tc_dat, method = "efron")

tc.1.Disag <- coxph(tc1_gap ~ +polity2 + area_1000_log + lmtnest + ELF
                    + strata(event_num) + cluster(ccode), data = tc1,
                    method = "efron")

tc23.Disag <- coxph(tc23_gap ~ Peace_wtc + Fighting_wtc  + Forceful_LinDecay
                     + Promoted_LinDecay + Peaceful_LinDecay +
                     Absorbed_LinDecay + polity2 + area_1000_log + lmtnest +
                     ELF + strata(event_num) + cluster(ccode), data = tc23,
                     method = "efron")

tc45.Disag <- coxph(tc45_gap ~ Peace_wtc + Fighting_wtc  + Forceful_LinDecay
                     + Promoted_LinDecay + Peaceful_LinDecay +
                     Absorbed_LinDecay + polity2 + area_1000_log + lmtnest +
                     ELF + strata(event_num) + cluster(ccode), data = tc45,
                     method = "efron")

tc6_more.Disag <- coxph(tc6_more_gap ~ Peace_wtc + Fighting_wtc  +
                        Forceful_LinDecay + Promoted_LinDecay +
                        Peaceful_LinDecay + Absorbed_LinDecay + polity2 +
                        area_1000_log + lmtnest + ELF + strata(event_num) +
                        cluster(ccode), data = tc6_more, method = "efron")

summary(tcAll.Disag)

stargazer(tcAll.Disag, tc.1.Disag, tc23.Disag, tc45.Disag, tc6_more.Disag,
          type = "latex",
          title = "PWP Gap Time Model Results - Disaggregated",
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
                               "ELF"),
          keep.stat = c("n"))

# Exponential Decay 20 year half life

tcAll.Exponential20 <- coxph(tc_dat_gap ~ Peace_wtc + Fighting_wtc  +
                             Forceful_ExpDecay_20 + GoodEnd_ExpDecay_20 +
                             polity2 + area_1000_log + lmtnest + ELF +
                             strata(event_num) + cluster(ccode), data = tc_dat,
                             method = "efron")

tc.1.Exponential20 <- coxph(tc1_gap ~ polity2 + area_1000_log + lmtnest + ELF +
                            + strata(event_num) + cluster(ccode), data = tc1,
                            method = "efron")

tc23.Exponential20 <- coxph(tc23_gap ~ Peace_wtc + Fighting_wtc  +
                              Forceful_ExpDecay_20 + GoodEnd_ExpDecay_20 +
                              polity2 + area_1000_log + lmtnest + ELF +
                              strata(event_num) + cluster(ccode), data = tc23,
                              method = "efron")

tc45.Exponential20 <- coxph(tc45_gap ~ Peace_wtc + Fighting_wtc  +
                              Forceful_ExpDecay_20 + GoodEnd_ExpDecay_20 +
                              polity2 + area_1000_log + lmtnest + ELF +
                              strata(event_num) + cluster(ccode), data = tc45,
                              method = "efron")

tc6_more.Exponential20 <- coxph(tc6_more_gap ~ Peace_wtc + Fighting_wtc  +
                                 Forceful_ExpDecay_20 + GoodEnd_ExpDecay_20 +
                                 polity2 + area_1000_log + lmtnest + ELF +
                                 strata(event_num) + cluster(ccode),
                                 data = tc6_more, method = "efron")

summary(tcAll.Exponential20)


stargazer(tcAll.Exponential20, tc.1.Exponential20, tc23.Exponential20,
          tc45.Exponential20, tc6_more.Exponential20,
          type = "latex",
          title = "PWP Gap Time Model Results - Exponential Decay with 20 year
          Half Life",
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
                               "ELF"),
          keep.stat = c("n"))

# Exponential Decay 10 year half life

tcAll.Exponential10 <- coxph(tc_dat_gap ~ Peace_wtc + Fighting_wtc  +
                             Forceful_ExpDecay_10 + GoodEnd_ExpDecay_10 +
                             polity2 + area_1000_log + lmtnest + ELF +
                             strata(event_num) + cluster(ccode), data = tc_dat,
                             method = "efron")

tc.1.Exponential10 <- coxph(tc1_gap ~ polity2 + area_1000_log + lmtnest + ELF +
                            strata(event_num) + cluster(ccode), data = tc1,
                            method = "efron")

tc23.Exponential10 <- coxph(tc23_gap ~ Peace_wtc + Fighting_wtc  +
                              Forceful_ExpDecay_10 + GoodEnd_ExpDecay_10 +
                              polity2 + area_1000_log + lmtnest + ELF +
                              strata(event_num) + cluster(ccode), data = tc23,
                              method = "efron")

tc45.Exponential10 <- coxph(tc45_gap ~ Peace_wtc + Fighting_wtc  +
                              Forceful_ExpDecay_10 + GoodEnd_ExpDecay_10 +
                              polity2 + area_1000_log + lmtnest + ELF +
                              strata(event_num) + cluster(ccode), data = tc45,
                              method = "efron")

tc6_more.Exponential10 <- coxph(tc6_more_gap ~ Peace_wtc + Fighting_wtc  +
                                 Forceful_ExpDecay_10 + GoodEnd_ExpDecay_10 +
                                 polity2 + area_1000_log + lmtnest + ELF +
                                 strata(event_num) + cluster(ccode),
                                 data = tc6_more, method = "efron")

summary(tcAll.Exponential10)


stargazer(tcAll.Exponential10, tc.1.Exponential10, tc23.Exponential10,
          tc45.Exponential10, tc6_more.Exponential10,
          type = "latex",
          title = "PWP Gap Time Model Results - Exponential Decay with 10 year
          Half Life",
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
                               "ELF"),
          keep.stat = c("n"))

# 20 Year Binary

tcAll.t20 <- coxph(tc_dat_gap ~ Peace_wtc + Fighting_wtc  + Forceful_t20 +
                   GoodEnd_t20 + polity2 + area_1000_log + lmtnest + ELF +
                   strata(event_num) + cluster(ccode), data = tc_dat,
                   method = "efron")

tc.1.t20 <- coxph(tc1_gap ~ polity2 + area_1000_log + lmtnest + ELF +
                  strata(event_num) + cluster(ccode), data = tc1,
                  method = "efron")

tc23.t20 <- coxph(tc23_gap ~ Peace_wtc + Fighting_wtc  + Forceful_t20 +
                    GoodEnd_t20 + polity2 + area_1000_log + lmtnest + ELF +
                    strata(event_num) + cluster(ccode), data = tc23,
                    method = "efron")

tc45.t20 <- coxph(tc45_gap ~ Peace_wtc + Fighting_wtc  + Forceful_t20 +
                    GoodEnd_t20 + polity2 + area_1000_log + lmtnest + ELF +
                    strata(event_num) + cluster(ccode), data = tc45,
                    method = "efron")

tc6_more.t20 <- coxph(tc6_more_gap ~ Peace_wtc + Fighting_wtc +
                       Forceful_t20 + GoodEnd_t20 + polity2 + area_1000_log +
                       lmtnest + ELF + strata(event_num) + cluster(ccode),
                       data = tc6_more, method = "efron")

summary(tcAll.t20)


stargazer(tcAll.t20, tc.1.t20, tc23.t20, tc45.t20, tc6_more.t20,
          type = "latex",
          title = "PWP Gap Time Model Results - 20 Year Binary Indicator",
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
                               "ELF"),
          keep.stat = c("n"))

# Hurdle Model - This seems so wrong. Its just lumping everything together.

hurdlemod <- hurdle(tc_tally ~ Peace_wtc + Fighting_wtc  + Forceful_LinDecay +
                    GoodEnd_LinDecay + polity2 + area_1000_log + lmtnest +
                    ELF | polity2 + area_1000_log + lmtnest + ELF,
                    data = tc_dat)

summary(hurdlemod)

# Censored Probit

