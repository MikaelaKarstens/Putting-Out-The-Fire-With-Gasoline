# Code for "Assessing State Making Strategies" by D. Lemke and M. Karstens

# Required packages for replication - Please install if needed =================

library(dplyr)
library(survival)
library(stargazer)
library(pscl)

# Settings =====================================================================

options(scipen = 50) # bias against scientific notation for convenience

# Loading data =================================================================

tc_dat <- read.csv("State-Making-Strategies-Replication-Data.csv")
names(tc_dat)


# Subset datasets ==============================================================

tc_multi <- tc_dat[tc_dat$event_num >= 2, ]

tc23 <- filter(tc_dat, event_num == 2 | event_num == 3)

tc45 <- filter(tc_dat, event_num == 4 | event_num == 5)

tc4_more <- filter(tc_dat, event_num >= 4)

tc6_more <- filter(tc_dat, event_num >= 6)

tc1 <- filter(tc_dat, event_num == 1)

tc2 <- filter(tc_dat, event_num == 2)

tc3 <- filter(tc_dat, event_num == 3)

tc4 <- filter(tc_dat, event_num == 4)

tc5 <- filter(tc_dat, event_num == 5)

tc6 <- filter(tc_dat, event_num == 6)

# Survival objects =============================================================

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

tc4_more_elapse <- Surv(tc4_more$alt_start, tc4_more$alt_stop,
                        tc4_more$new_tc_dummy)
tc4_more_gap <- Surv(tc4_more$start, tc4_more$stop, tc4_more$new_tc_dummy)

tc6_more_elapse <- Surv(tc6_more$alt_start, tc6_more$alt_stop,
                          tc6_more$new_tc_dummy)
tc6_more_gap <- Surv(tc6_more$start, tc6_more$stop, tc6_more$new_tc_dummy)

# Code for Tables in the Main Text of the Article ==============================
  # Table 1 - 20 Year Linear Decay PWP Gap Time Models =========================

tc_pwp <- coxph(tc_dat_gap ~ peace_wtc + fighting_wtc  + force_lin_decay +
                favor_tc_lin_decay + absorb_lin_decay + lib_dem_vdem +
                area_1000_log + lmtnest + elf + tc_tally + strata(event_num) +
                cluster(ccode), data = tc_dat, method = "efron")

tc1_pwp <- coxph(tc1_gap ~ lib_dem_vdem + area_1000_log + lmtnest + elf +
                 strata(event_num) + cluster(ccode), data = tc1,
                 method = "efron")

tc23_pwp_lin20 <- coxph(tc23_gap ~ peace_wtc + fighting_wtc  +
                        force_lin_decay + favor_tc_lin_decay +
                        absorb_lin_decay + lib_dem_vdem + area_1000_log +
                        lmtnest + elf + strata(event_num) + cluster(ccode),
                        data = tc23, method = "efron")

tc45_pwp_lin20 <- coxph(tc45_gap ~ peace_wtc + fighting_wtc +
                        force_lin_decay + favor_tc_lin_decay +
                        absorb_lin_decay + lib_dem_vdem + area_1000_log +
                        lmtnest + elf + strata(event_num) + cluster(ccode),
                        data = tc45, method = "efron")

tc6_more_pwp_lin20 <- coxph(tc6_more_gap ~ peace_wtc + fighting_wtc  +
                            force_lin_decay + favor_tc_lin_decay +
                            absorb_lin_decay + lib_dem_vdem + area_1000_log +
                            lmtnest + elf + strata(event_num) +
                            cluster(ccode), data = tc6_more, method = "efron")

stargazer(tc_pwp, tc1_pwp, tc23_pwp_lin20, tc45_pwp_lin20, tc6_more_pwp_lin20,
          type = "latex",
          title = "PWP Gap Time Model Results",
          model.numbers = F,
          column.labels = c("All TCs",
                            "First TC",
                            " TC 2 or 3",
                            "TC 4 or 5",
                            "TC 6+"),
          dep.var.labels = c("(1)", "(2)", "(3)", "(4)", "(5)"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Favorable Outcome for State",
                               "Favorable Outcome for TC",
                               "Absorbed",
                               "Liberal Democracy Index",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))

  # Table 2 - 20 Year Linear Decay PWP Elapsed Time Models =====================

tc_pwp_elapse <- coxph(tc_dat_elapse ~ peace_wtc + fighting_wtc  +
                       force_lin_decay + favor_tc_lin_decay + absorb_lin_decay +
                       lib_dem_vdem + area_1000_log + lmtnest + elf + tc_tally +
                       strata(event_num) + cluster(ccode), data = tc_dat,
                       method = "efron")

tc1_pwp_elapse <- coxph(tc1_elapse ~ lib_dem_vdem + area_1000_log + lmtnest +
                        elf + strata(event_num) + cluster(ccode), data = tc1,
                        method = "efron")

tc23_pwp_elapse <- coxph(tc23_elapse ~ peace_wtc + fighting_wtc  +
                         force_lin_decay + favor_tc_lin_decay +
                         absorb_lin_decay + lib_dem_vdem + area_1000_log +
                         lmtnest + elf + strata(event_num) + cluster(ccode),
                         data = tc23, method = "efron")


tc4_more_pwp_elapse <- coxph(tc4_more_elapse ~ peace_wtc + fighting_wtc +
                             force_lin_decay + favor_tc_lin_decay +
                             absorb_lin_decay + lib_dem_vdem +
                             area_1000_log + lmtnest + elf +
                             strata(event_num) + cluster(ccode),
                             data = tc4_more, method = "efron")

stargazer(tc_pwp_elapse, tc1_pwp_elapse, tc23_pwp_elapse, tc4_more_pwp_elapse,
          type = "latex",
          title = "PWP Elapsed Time Model Results",
          model.numbers = F,
          column.labels = c("All TCs",
                            "First TC",
                            " TC 2 or 3",
                            "TC 4+"),
          dep.var.labels = c("(1)", "(2)", "(3)", "(4)"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Favorable Outcome for State",
                               "Favorable Outcome for TC",
                               "Absorbed",
                               "Liberal Democracy Index",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))

# Appendix B - Alternative Variables ===========================================
  # B1 - PWP Gap Time Model 20 Year Linear Decay - Disaggregated ===============

tc_pwp_disag <- coxph(tc_dat_gap ~ peace_wtc + fighting_wtc  + force_lin_decay +
                      absorb_lin_decay + peaceful_lin_decay +
                      promoted_lin_decay + lib_dem_vdem + area_1000_log +
                      lmtnest + elf + tc_tally + strata(event_num) +
                      cluster(ccode), data = tc_dat, method = "efron")

tc1_pwp_disag <- coxph(tc1_gap ~ lib_dem_vdem + area_1000_log + lmtnest + elf +
                       strata(event_num) + cluster(ccode), data = tc1,
                       method = "efron")

tc23_pwp_disag <- coxph(tc23_gap ~ peace_wtc + fighting_wtc  +
                        force_lin_decay + absorb_lin_decay +
                        peaceful_lin_decay + promoted_lin_decay + lib_dem_vdem +
                        area_1000_log + lmtnest + elf + strata(event_num) +
                        cluster(ccode), data = tc23, method = "efron")

tc45_pwp_disag <- coxph(tc45_gap ~ peace_wtc + fighting_wtc +
                        force_lin_decay + absorb_lin_decay +
                        peaceful_lin_decay + promoted_lin_decay +
                        lib_dem_vdem + area_1000_log + lmtnest + elf +
                        strata(event_num) + cluster(ccode), data = tc45,
                        method = "efron")

tc6_more_pwp_disag <- coxph(tc6_more_gap ~ peace_wtc + fighting_wtc  +
                            force_lin_decay + absorb_lin_decay +
                            peaceful_lin_decay + promoted_lin_decay +
                            lib_dem_vdem + area_1000_log + lmtnest + elf +
                            strata(event_num) + cluster(ccode), data = tc6_more,
                            method = "efron")

stargazer(tc_pwp_disag, tc1_pwp_disag, tc23_pwp_disag, tc45_pwp_disag,
          tc6_more_pwp_disag, type = "latex",
          title = "PWP Gap Time Model Results - Disaggregated",
          model.numbers = F,
          column.labels = c("All TCs",
                            "First TC",
                            " TC 2 or 3",
                            "TC 4 or 5",
                            "TC 6+"),
          dep.var.labels = c("(1)", "(2)", "(3)", "(4)", "(5)"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Forceful Reintegration",
                               "Absorbed",
                               "Peaceful Reintegration",
                               "Promoted to Sovereign",
                               "Liberal Democracy Index",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))


  # B2 - Using Absorb and Forceful for Favorable for Sov. State ================

tc_pwp_ss <- coxph(tc_dat_gap ~ peace_wtc + fighting_wtc  + favor_ss_lin_decay +
                   favor_tc_lin_decay + lib_dem_vdem + area_1000_log + lmtnest +
                   elf + tc_tally + strata(event_num) + cluster(ccode),
                   data = tc_dat, method = "efron")

tc1_pwp_ss <- coxph(tc1_gap ~ lib_dem_vdem + area_1000_log + lmtnest + elf +
                    strata(event_num) + cluster(ccode), data = tc1,
                    method = "efron")

tc23_pwp_ss <- coxph(tc23_gap ~ peace_wtc + favor_ss_lin_decay  +
                     favor_tc_lin_decay + lib_dem_vdem +
                     area_1000_log + lmtnest + elf +
                     strata(event_num) + cluster(ccode), data = tc23,
                     method = "efron")

tc45_pwp_ss <- coxph(tc45_gap ~ peace_wtc + favor_ss_lin_decay +
                     favor_tc_lin_decay + lib_dem_vdem +
                     area_1000_log + lmtnest + elf +
                     strata(event_num) + cluster(ccode), data = tc45,
                     method = "efron")

tc6_more_pwp_ss <- coxph(tc6_more_gap ~ peace_wtc + favor_ss_lin_decay  +
                         favor_tc_lin_decay +
                         lib_dem_vdem + area_1000_log + lmtnest + elf +
                         strata(event_num) + cluster(ccode),
                         data = tc6_more, method = "efron")

stargazer(tc_pwp_ss, tc1_pwp_ss, tc23_pwp_ss, tc45_pwp_ss, tc6_more_pwp_ss,
          type = "latex",
          title = "PWP Gap Time Model Results - Favor. Sovereign States",
          model.numbers = F,
          column.labels = c("All TCs",
                            "First TC",
                            " TC 2 or 3",
                            "TC 4 or 5",
                            "TC 6+"),
          dep.var.labels = c("(1)", "(2)", "(3)", "(4)", "(5)"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Favorable Outcome for State",
                               "Favorable Outcome for TC",
                               "Liberal Democracy Index",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))

  # B3 - PWP Gap Time Model 20 Year Linear Decay - State Failure Control =======

tc_pwp_sf <- coxph(tc_dat_gap ~ peace_wtc + fighting_wtc  + force_lin_decay +
                  favor_tc_lin_decay + absorb_lin_decay + lib_dem_vdem +
                  area_1000_log + lmtnest + elf + state_fail + tc_tally +
                  strata(event_num) + cluster(ccode), data = tc_dat,
                  method = "efron")

tc1_pwp_sf <- coxph(tc1_gap ~ lib_dem_vdem + area_1000_log + lmtnest + elf +
                   strata(event_num) + cluster(ccode), data = tc1,
                 method = "efron")

tc23_pwp_sf <- coxph(tc23_gap ~ peace_wtc + fighting_wtc  +
                     force_lin_decay + favor_tc_lin_decay +
                     absorb_lin_decay + lib_dem_vdem + area_1000_log +
                     lmtnest + elf + state_fail + strata(event_num) +
                     cluster(ccode), data = tc23, method = "efron")

tc45_pwp_sf <- coxph(tc45_gap ~ peace_wtc + fighting_wtc +
                     force_lin_decay + favor_tc_lin_decay +
                     absorb_lin_decay + lib_dem_vdem + area_1000_log +
                     lmtnest + elf + state_fail + strata(event_num) +
                     cluster(ccode), data = tc45, method = "efron")

tc6_more_pwp_sf <- coxph(tc6_more_gap ~ peace_wtc + fighting_wtc  +
                         force_lin_decay + favor_tc_lin_decay +
                         absorb_lin_decay + lib_dem_vdem + area_1000_log +
                         lmtnest + elf + state_fail + strata(event_num) +
                         cluster(ccode), data = tc6_more, method = "efron")

stargazer(tc_pwp_sf, tc1_pwp_sf, tc23_pwp_sf, tc45_pwp_sf, tc6_more_pwp_sf,
          type = "latex",
          title = "PWP Gap Time Model Results - State Failure",
          model.numbers = F,
          column.labels = c("All TCs",
                            "First TC",
                            " TC 2 or 3",
                            "TC 4 or 5",
                            "TC 6+"),
          dep.var.labels = c("(1)", "(2)", "(3)", "(4)", "(5)"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Favorable Outcome for State",
                               "Favorable Outcome for TC",
                               "Absorbed",
                               "Liberal Democracy Index",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "State Failure",
                               "TC Tally"),
          keep.stat = c("n"))

  # B4 - 20 Year Linear Decay PWP Gap Time Models - Polity 2 ===================

tc_pwp_polity <- coxph(tc_dat_gap ~ peace_wtc + fighting_wtc  +
                       force_lin_decay + favor_tc_lin_decay + absorb_lin_decay +
                       polity2 + area_1000_log + lmtnest + elf + tc_tally +
                       strata(event_num) + cluster(ccode), data = tc_dat,
                       method = "efron")

tc1_pwp_polity <- coxph(tc1_gap ~ polity2 + area_1000_log + lmtnest + elf +
                        strata(event_num) + cluster(ccode), data = tc1,
                        method = "efron")

tc23_pwp_polity <- coxph(tc23_gap ~ peace_wtc + fighting_wtc  +
                         force_lin_decay + favor_tc_lin_decay +
                         absorb_lin_decay + polity2 + area_1000_log +
                         lmtnest + elf + strata(event_num) + cluster(ccode),
                         data = tc23, method = "efron")

tc45_pwp_polity <- coxph(tc45_gap ~ peace_wtc + fighting_wtc +
                         force_lin_decay + favor_tc_lin_decay +
                         absorb_lin_decay + polity2 + area_1000_log +
                         lmtnest + elf + strata(event_num) + cluster(ccode),
                         data = tc45, method = "efron")

tc6_more_pwp_polity <- coxph(tc6_more_gap ~ peace_wtc + fighting_wtc  +
                              force_lin_decay + favor_tc_lin_decay +
                              absorb_lin_decay + polity2 + area_1000_log +
                              lmtnest + elf + strata(event_num) +
                              cluster(ccode), data = tc6_more, method = "efron")

stargazer(tc_pwp_polity, tc1_pwp_polity, tc23_pwp_polity, tc45_pwp_polity,
          tc6_more_pwp_polity, type = "latex",
          title = "PWP Gap Time Model Results - Polity 2",
          model.numbers = F,
          column.labels = c("All TCs",
                            "First TC",
                            " TC 2 or 3",
                            "TC 4 or 5",
                            "TC 6+"),
          dep.var.labels = c("(1)", "(2)", "(3)", "(4)", "(5)"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Favorable Outcome for State",
                               "Favorable Outcome for TC",
                               "Absorbed",
                               "Polity Score",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))

  # B5 - PWP Gap Time Model 20 Year Linear Decay - Electoral Democracy =========

tc_pwp_elec <- coxph(tc_dat_gap ~ peace_wtc + fighting_wtc  + force_lin_decay +
                     favor_tc_lin_decay + absorb_lin_decay + elec_dem_vdem +
                     area_1000_log + lmtnest + elf + tc_tally +
                     strata(event_num) + cluster(ccode), data = tc_dat,
                     method = "efron")

tc1_pwp_elec <- coxph(tc1_gap ~ elec_dem_vdem + area_1000_log + lmtnest + elf +
                      strata(event_num) + cluster(ccode), data = tc1,
                      method = "efron")

tc23_pwp_elec <- coxph(tc23_gap ~ peace_wtc + fighting_wtc  + force_lin_decay +
                       favor_tc_lin_decay + absorb_lin_decay + elec_dem_vdem +
                       area_1000_log + lmtnest + elf + strata(event_num) +
                       cluster(ccode), data = tc23, method = "efron")

tc45_pwp_elec <- coxph(tc45_gap ~ peace_wtc + fighting_wtc + force_lin_decay +
                       favor_tc_lin_decay + absorb_lin_decay + elec_dem_vdem +
                       area_1000_log + lmtnest + elf + strata(event_num) +
                       cluster(ccode), data = tc45, method = "efron")

tc6_more_pwp_elec <- coxph(tc6_more_gap ~ peace_wtc + fighting_wtc  +
                             force_lin_decay + favor_tc_lin_decay +
                             absorb_lin_decay + elec_dem_vdem + area_1000_log +
                             lmtnest + elf + strata(event_num) +
                             cluster(ccode), data = tc6_more, method = "efron")

stargazer(tc_pwp_elec, tc1_pwp_elec, tc23_pwp_elec, tc45_pwp_elec,
          tc6_more_pwp_elec, type = "latex",
          title = "PWP Gap Time Model Results - Electoral Democracy",
          model.numbers = F,
          column.labels = c("All TCs",
                            "First TC",
                            " TC 2 or 3",
                            "TC 4 or 5",
                            "TC 6+"),
          dep.var.labels = c("(1)", "(2)", "(3)", "(4)", "(5)"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Favorable Outcome for State",
                               "Favorable Outcome for TC",
                               "Absorbed",
                               "Electoral Democracy Index",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))


# Appendix C - Decay Function Alternatives =====================================
  # C1 - PWP Gap Time Models - 10 Year Binary ==================================

tc_pwp_t10 <- coxph(tc_dat_gap ~ peace_wtc + fighting_wtc  + force_t10 +
                    favor_tc_t10 + absorb_t10 + lib_dem_vdem +
                    area_1000_log + lmtnest + elf + tc_tally +
                    strata(event_num) + cluster(ccode), data = tc_dat,
                    method = "efron")

tc1_pwp_t10 <- coxph(tc1_gap ~ lib_dem_vdem + area_1000_log + lmtnest + elf +
                     strata(event_num) + cluster(ccode), data = tc1,
                     method = "efron")

tc23_pwp_t10 <- coxph(tc23_gap ~ peace_wtc + fighting_wtc  + force_t10 +
                      favor_tc_t10 + absorb_t10 + lib_dem_vdem + area_1000_log +
                      lmtnest + elf + strata(event_num) + cluster(ccode),
                      data = tc23, method = "efron")

tc45_pwp_t10 <- coxph(tc45_gap ~ peace_wtc + fighting_wtc + force_t10 +
                      favor_tc_t10 + absorb_t10 + lib_dem_vdem +
                      area_1000_log + lmtnest + elf + strata(event_num) +
                      cluster(ccode), data = tc45, method = "efron")

tc6_more_pwp_t10 <- coxph(tc6_more_gap ~ peace_wtc + fighting_wtc  + force_t10 +
                          favor_tc_t10 + absorb_t10 + lib_dem_vdem +
                          area_1000_log + lmtnest + elf + strata(event_num) +
                          cluster(ccode), data = tc6_more, method = "efron")

stargazer(tc_pwp_t10, tc1_pwp_t10, tc23_pwp_t10, tc45_pwp_t10, tc6_more_pwp_t10,
          type = "latex",
          title = "PWP Gap Time Model Results - 10 Year Binary",
          model.numbers = F,
          column.labels = c("All TCs",
                            "First TC",
                            " TC 2 or 3",
                            "TC 4 or 5",
                            "TC 6+"),
          dep.var.labels = c("(1)", "(2)", "(3)", "(4)", "(5)"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Favorable Outcome for State",
                               "Favorable Outcome for TC",
                               "Absorbed",
                               "Liberal Democracy Index",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))

  # C2 - PWP Gap Time Models - 20 Year Binary ==================================

tc_pwp_t20 <- coxph(tc_dat_gap ~ peace_wtc + fighting_wtc  + force_t20 +
                    favor_tc_t20 + absorb_t20 + lib_dem_vdem +
                    area_1000_log + lmtnest + elf + tc_tally +
                    strata(event_num) + cluster(ccode), data = tc_dat,
                    method = "efron")

tc1_pwp_t20 <- coxph(tc1_gap ~ lib_dem_vdem + area_1000_log + lmtnest + elf +
                     strata(event_num) + cluster(ccode), data = tc1,
                     method = "efron")

tc23_pwp_t20 <- coxph(tc23_gap ~ peace_wtc + fighting_wtc  + force_t20 +
                      favor_tc_t20 + absorb_t20 + lib_dem_vdem +
                      area_1000_log + lmtnest + elf + strata(event_num) +
                      cluster(ccode), data = tc23, method = "efron")

tc45_pwp_t20 <- coxph(tc45_gap ~ peace_wtc + fighting_wtc + force_t20 +
                      favor_tc_t20 + absorb_t20 + lib_dem_vdem +
                      area_1000_log + lmtnest + elf + strata(event_num) +
                      cluster(ccode), data = tc45, method = "efron")

tc6_more_pwp_t20 <- coxph(tc6_more_gap ~ peace_wtc + fighting_wtc + force_t20 +
                          favor_tc_t20 + absorb_t20 + lib_dem_vdem +
                          area_1000_log + lmtnest + elf + strata(event_num) +
                          cluster(ccode), data = tc6_more, method = "efron")

stargazer(tc_pwp_t20, tc1_pwp_t20, tc23_pwp_t20, tc45_pwp_t20, tc6_more_pwp_t20,
          type = "latex",
          title = "PWP Gap Time Model Results - 20 Year Binary",
          model.numbers = F,
          column.labels = c("All TCs",
                            "First TC",
                            " TC 2 or 3",
                            "TC 4 or 5",
                            "TC 6+"),
          dep.var.labels = c("(1)", "(2)", "(3)", "(4)", "(5)"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Favorable Outcome for State",
                               "Favorable Outcome for TC",
                               "Absorbed",
                               "Liberal Democracy Index",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))

  # C3 - PWP Gap Time Models - Exponential 5 Year Halflife =====================

tc_pwp_exp5 <- coxph(tc_dat_gap ~ peace_wtc + fighting_wtc  +
                     force_exp_decay_5 + favor_tc_exp_decay_5 + absorb_exp_5 +
                     lib_dem_vdem + area_1000_log + lmtnest + elf + tc_tally +
                     strata(event_num) + cluster(ccode), data = tc_dat,
                     method = "efron")

tc1_pwp_exp5 <- coxph(tc1_gap ~ lib_dem_vdem + area_1000_log + lmtnest + elf +
                      strata(event_num) + cluster(ccode), data = tc1,
                      method = "efron")

tc23_pwp_exp5 <- coxph(tc23_gap ~ peace_wtc + fighting_wtc  +
                       force_exp_decay_5 + favor_tc_exp_decay_5 +
                       absorb_exp_5 + lib_dem_vdem + area_1000_log +
                       lmtnest + elf + strata(event_num) + cluster(ccode),
                       data = tc23, method = "efron")

tc45_pwp_exp5 <- coxph(tc45_gap ~ peace_wtc + fighting_wtc +
                       force_exp_decay_5 + favor_tc_exp_decay_5 +
                       absorb_exp_5 + lib_dem_vdem + area_1000_log +
                       lmtnest + elf + strata(event_num) + cluster(ccode),
                       data = tc45, method = "efron")

tc6_more_pwp_exp5 <- coxph(tc6_more_gap ~ peace_wtc + fighting_wtc  +
                           force_exp_decay_5 + favor_tc_exp_decay_5 +
                           absorb_exp_5 + lib_dem_vdem + area_1000_log +
                           lmtnest + elf + strata(event_num) +
                           cluster(ccode), data = tc6_more, method = "efron")

stargazer(tc_pwp_exp5, tc1_pwp_exp5, tc23_pwp_exp5, tc45_pwp_exp5,
          tc6_more_pwp_exp5, type = "latex",
          title = "PWP Gap Time Model - Exponential Decay 5-Year Halflife",
          model.numbers = F,
          column.labels = c("All TCs",
                            "First TC",
                            " TC 2 or 3",
                            "TC 4 or 5",
                            "TC 6+"),
          dep.var.labels = c("(1)", "(2)", "(3)", "(4)", "(5)"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Favorable Outcome for State",
                               "Favorable Outcome for TC",
                               "Absorbed",
                               "Liberal Democracy Index",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))

  # C4 - PWP Gap Time Models - Exponential 10 Year Halflife ====================

tc_pwp_exp10 <- coxph(tc_dat_gap ~ peace_wtc + fighting_wtc  +
                      force_exp_decay_10 + favor_tc_exp_decay_10 +
                      absorb_exp_10 + lib_dem_vdem + area_1000_log + lmtnest +
                      elf + tc_tally + strata(event_num) + cluster(ccode),
                      data = tc_dat, method = "efron")

tc1_pwp_exp10 <- coxph(tc1_gap ~ lib_dem_vdem + area_1000_log + lmtnest + elf +
                       strata(event_num) + cluster(ccode), data = tc1,
                       method = "efron")

tc23_pwp_exp10 <- coxph(tc23_gap ~ peace_wtc + fighting_wtc  +
                        force_exp_decay_10 + favor_tc_exp_decay_10 +
                        absorb_exp_10 + lib_dem_vdem + area_1000_log +
                        lmtnest + elf + strata(event_num) + cluster(ccode),
                        data = tc23, method = "efron")

tc45_pwp_exp10 <- coxph(tc45_gap ~ peace_wtc + fighting_wtc +
                        force_exp_decay_10 + favor_tc_exp_decay_10 +
                        absorb_exp_10 + lib_dem_vdem + area_1000_log +
                        lmtnest + elf + strata(event_num) + cluster(ccode),
                        data = tc45, method = "efron")

tc6_more_pwp_exp10 <- coxph(tc6_more_gap ~ peace_wtc + fighting_wtc  +
                            force_exp_decay_10 + favor_tc_exp_decay_10 +
                            absorb_exp_10 + lib_dem_vdem + area_1000_log +
                            lmtnest + elf + strata(event_num) +
                            cluster(ccode), data = tc6_more, method = "efron")

stargazer(tc_pwp_exp10, tc1_pwp_exp10, tc23_pwp_exp10, tc45_pwp_exp10,
          tc6_more_pwp_exp10, type = "latex",
          title = "PWP Gap Time Model - Exponential Decay 10 Year Halflife",
          model.numbers = F,
          column.labels = c("All TCs",
                            "First TC",
                            " TC 2 or 3",
                            "TC 4 or 5",
                            "TC 6+"),
          dep.var.labels = c("(1)", "(2)", "(3)", "(4)", "(5)"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Favorable Outcome for State",
                               "Favorable Outcome for TC",
                               "Absorbed",
                               "Liberal Democracy Index",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))

  # C5 - PWP Gap Time Models - Exponential 15 Year Halflife ====================

tc_pwp_exp15 <- coxph(tc_dat_gap ~ peace_wtc + fighting_wtc  +
                      force_exp_decay_15 + favor_tc_exp_decay_15 +
                      absorb_exp_15 + lib_dem_vdem + area_1000_log + lmtnest +
                      elf + tc_tally + strata(event_num) +
                      cluster(ccode), data = tc_dat, method = "efron")

tc1_pwp_exp15 <- coxph(tc1_gap ~ lib_dem_vdem + area_1000_log + lmtnest + elf +
                       strata(event_num) + cluster(ccode), data = tc1,
                       method = "efron")

tc23_pwp_exp15 <- coxph(tc23_gap ~ peace_wtc + fighting_wtc  +
                        force_exp_decay_15 + favor_tc_exp_decay_15 +
                        absorb_exp_15 + lib_dem_vdem + area_1000_log +
                        lmtnest + elf + strata(event_num) + cluster(ccode),
                        data = tc23, method = "efron")

tc45_pwp_exp15 <- coxph(tc45_gap ~ peace_wtc + fighting_wtc +
                        force_exp_decay_15 + favor_tc_exp_decay_15 +
                        absorb_exp_15 + lib_dem_vdem + area_1000_log +
                        lmtnest + elf + strata(event_num) + cluster(ccode),
                        data = tc45, method = "efron")

tc6_more_pwp_exp15 <- coxph(tc6_more_gap ~ peace_wtc + fighting_wtc  +
                            force_exp_decay_15 + favor_tc_exp_decay_15 +
                            absorb_exp_15 + lib_dem_vdem + area_1000_log +
                            lmtnest + elf + strata(event_num) +
                            cluster(ccode), data = tc6_more, method = "efron")

stargazer(tc_pwp_exp15, tc1_pwp_exp15, tc23_pwp_exp15, tc45_pwp_exp15,
          tc6_more_pwp_exp15, type = "latex",
          title = "PWP Gap Time Model - Exponential Decay 15 Year Halflife",
          model.numbers = F,
          column.labels = c("All TCs",
                            "First TC",
                            " TC 2 or 3",
                            "TC 4 or 5",
                            "TC 6+"),
          dep.var.labels = c("(1)", "(2)", "(3)", "(4)", "(5)"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Favorable Outcome for State",
                               "Favorable Outcome for TC",
                               "Absorbed",
                               "Liberal Democracy Index",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))

  # C6 - PWP Gap Time Models - Exponential 20 Year Halflife ====================

tc_pwp_exp20 <- coxph(tc_dat_gap ~ peace_wtc + fighting_wtc  +
                      force_exp_decay_20 + favor_tc_exp_decay_20 +
                      absorb_exp_20 + lib_dem_vdem + area_1000_log + lmtnest +
                      elf + tc_tally + strata(event_num) + cluster(ccode),
                      data = tc_dat, method = "efron")

tc1_pwp_exp20 <- coxph(tc1_gap ~ lib_dem_vdem + area_1000_log + lmtnest + elf +
                       strata(event_num) + cluster(ccode), data = tc1,
                       method = "efron")

tc23_pwp_exp20 <- coxph(tc23_gap ~ peace_wtc + fighting_wtc  +
                        force_exp_decay_20 + favor_tc_exp_decay_20 +
                        absorb_exp_20 + lib_dem_vdem + area_1000_log +
                        lmtnest + elf + strata(event_num) + cluster(ccode),
                        data = tc23, method = "efron")

tc45_pwp_exp20 <- coxph(tc45_gap ~ peace_wtc + fighting_wtc +
                        force_exp_decay_20 + favor_tc_exp_decay_20 +
                        absorb_exp_20 + lib_dem_vdem + area_1000_log +
                        lmtnest + elf + strata(event_num) + cluster(ccode),
                        data = tc45, method = "efron")

tc6_more_pwp_exp20 <- coxph(tc6_more_gap ~ peace_wtc + fighting_wtc  +
                            force_exp_decay_20 + favor_tc_exp_decay_20 +
                            absorb_exp_20 + lib_dem_vdem + area_1000_log +
                            lmtnest + elf + strata(event_num) +
                            cluster(ccode), data = tc6_more, method = "efron")

stargazer(tc_pwp_exp20, tc1_pwp_exp20, tc23_pwp_exp20, tc45_pwp_exp20,
          tc6_more_pwp_exp20, type = "latex",
          title = "PWP Gap Time Model - Exponential Decay 20 Year Halflife",
          model.numbers = F,
          column.labels = c("All TCs",
                            "First TC",
                            " TC 2 or 3",
                            "TC 4 or 5",
                            "TC 6+"),
          dep.var.labels = c("(1)", "(2)", "(3)", "(4)", "(5)"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Favorable Outcome for State",
                               "Favorable Outcome for TC",
                               "Absorbed",
                               "Liberal Democracy Index",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))

# Appendix D - Alternative Modeling Strategies =================================
  # D1 - Anderson-Gill Model - AG ==============================================

tc_ag <- coxph(tc_dat_gap ~ peace_wtc + fighting_wtc +
               force_lin_decay + favor_tc_lin_decay + absorb_lin_decay +
               lib_dem_vdem +
               area_1000_log + lmtnest + elf + tc_tally +
               cluster(ccode), data = tc_dat, method = "efron")

summary(tc_ag)

tc_ag_multi <- coxph(tc_multi_gap ~ peace_wtc + fighting_wtc +
                     force_lin_decay + favor_tc_lin_decay + absorb_lin_decay +
                     lib_dem_vdem + area_1000_log + lmtnest + elf + tc_tally +
                     cluster(ccode), data = tc_multi, method = "efron")

summary(tc_ag_multi)

stargazer(tc_ag, tc_ag_multi,
          type = "latex",
          title = "Anderson Gill Variance Correction Model",
          dep.var.labels.include = F,
          model.numbers = F,
          column.labels = c("All TCs",
                            "TC 2+"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Forceful Reintegration",
                               "Favorable Outcome for TC",
                               "Absorbed",
                               "Liberal Democracy Index",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))

  # D2 - WLW Model =============================================================

library(plyr)

tc_dat_expand <- tc_dat[rep(seq_len(nrow(tc_dat)), each =
                              max(tc_dat$event_num)), ]

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

tc_dat_wlw <- coxph(tc_dat_expand_s ~ peace_wtc + fighting_wtc  +
                    force_lin_decay + favor_tc_lin_decay + absorb_lin_decay +
                    lib_dem_vdem + area_1000_log + lmtnest + elf +
                    tc_tally + strata(event_num) + cluster(ccode),
                    data = tc_dat_expand, method = "efron")

summary(tc_dat_wlw)

library(dplyr)

stargazer(tc_dat_wlw,
          type = "latex",
          title = "WLW Marginal Elapsed Time Model",
          dep.var.labels.include = F,
          model.numbers = F,
          column.labels = c("Elapsed Time"),
          covariate.labels = c("Peace w/ TC",
                               "Fighting w/ TC",
                               "Forceful Reintegration",
                               "Favorable Outcome for TC",
                               "Absorbed",
                               "Liberal Democracy Index",
                               "State Area (logged)",
                               "Mountains",
                               "ELF",
                               "TC Tally"),
          keep.stat = c("n"))

  # D3 - Conditional Frailty Model - Gamma =====================================

tc_dat_gamma <- coxph(tc_dat_gap ~ peace_wtc + fighting_wtc  +
                      force_lin_decay + favor_tc_lin_decay +
                      absorb_lin_decay + lib_dem_vdem + area_1000_log +
                      lmtnest + elf + frailty.gamma(ccode) +
                      strata(event_num), data = tc_dat, method = "efron")

tc1_frail_gamma <- coxph(tc1_gap ~ lib_dem_vdem + area_1000_log + lmtnest +
                         elf + frailty.gamma(ccode) + strata(event_num),
                         data = tc1, method = "efron")

tc23_frail_gamma <- coxph(tc23_gap ~ peace_wtc + fighting_wtc  +
                          force_lin_decay + favor_tc_lin_decay +
                          absorb_lin_decay + lib_dem_vdem + area_1000_log +
                          lmtnest + elf + frailty.gamma(ccode) +
                          strata(event_num), data = tc23,
                          method = "efron")

tc45_frail_gamma <- coxph(tc45_gap ~ peace_wtc + fighting_wtc +
                          force_lin_decay + favor_tc_lin_decay +
                          absorb_lin_decay + lib_dem_vdem + area_1000_log +
                          lmtnest + elf + frailty.gamma(ccode) +
                          strata(event_num),
                          data = tc45, method = "efron")

tc6_more_frail_gamma <- coxph(tc6_more_gap ~ peace_wtc + fighting_wtc  +
                              force_lin_decay + favor_tc_lin_decay +
                              absorb_lin_decay + lib_dem_vdem + area_1000_log +
                              lmtnest + elf + + frailty.gamma(ccode) +
                              strata(event_num),
                              data = tc6_more, method = "efron")

summary(tc_dat_gamma)
summary(tc1_frail_gamma)
summary(tc23_frail_gamma)
summary(tc45_frail_gamma)
summary(tc6_more_frail_gamma)

  # D4 - Conditional Frailty Model - Gaussian ==================================

tc_dat_frail_gaus <- coxph(tc_dat_gap ~ peace_wtc + fighting_wtc  +
                           force_lin_decay + favor_tc_lin_decay +
                           absorb_lin_decay + lib_dem_vdem + area_1000_log +
                           lmtnest + elf + frailty.gaussian(ccode) +
                           strata(event_num), data = tc_dat, method = "efron")

tc1_frail_gaus <- coxph(tc1_gap ~ lib_dem_vdem + area_1000_log + lmtnest +
                        elf + frailty.gaussian(ccode) + strata(event_num),
                        data = tc1, method = "efron")

tc23_frail_gaus <- coxph(tc23_gap ~ peace_wtc + fighting_wtc  +
                         force_lin_decay + favor_tc_lin_decay +
                         absorb_lin_decay + lib_dem_vdem + area_1000_log +
                         lmtnest + elf + frailty.gaussian(ccode) +
                         strata(event_num), data = tc23,
                         method = "efron")

tc45_frail_gaus <- coxph(tc45_gap ~ peace_wtc + fighting_wtc +
                         force_lin_decay + favor_tc_lin_decay +
                         absorb_lin_decay + lib_dem_vdem + area_1000_log +
                         lmtnest + elf + frailty.gaussian(ccode) +
                         strata(event_num),
                         data = tc45, method = "efron")

tc6_more_frail_gaus <- coxph(tc6_more_gap ~ peace_wtc + fighting_wtc  +
                             force_lin_decay + favor_tc_lin_decay +
                             absorb_lin_decay + lib_dem_vdem + area_1000_log +
                             lmtnest + elf + frailty.gaussian(ccode) +
                             strata(event_num), data = tc6_more,
                             method = "efron")

summary(tc_dat_frail_gaus)
summary(tc1_frail_gaus)
summary(tc23_frail_gaus)
summary(tc45_frail_gaus)
summary(tc6_more_frail_gaus)
  # D5 - Censored Probit - In Stata do-file ====================================
  # D6 - Zero-Inflated Negative Binomial - In Stata do-file ====================
  # D7 - Hurdle Hurdle Model - In Stata do-file ==============================