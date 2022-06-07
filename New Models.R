# New Models for possible inclusion

dat <- read.csv("Putting-Out-The-Fire-With-Gasoline.csv")

HS <- read.csv("Hanson-Sigman-Capacity-Data.csv")

library(dplyr)

merg2 <- select(HS, unique_id, hs_capacity)

dat <- left_join(dat, merg2)


# correlate HS with Devel

test <- filter(dat, year > 1959)

test <- select(dat, hs_capacity, Devel)

cor(test$hs_capacity, test$Devel, use = "complete.obs", method = c("pearson"))

library(peacesciencer)

# Rerunning Model with Thom. Rivalry and Capacity


tc_pwp_all <- coxph(tc_dat_gap ~ peace_wtc + fighting_wtc  + force_lin_decay +
                      favor_tc_lin_decay + absorb_lin_decay + lib_dem_vdem +
                      area_1000_log + lmtnest + elf + rivals_thom + Devel +
                      tc_tally +
                      strata(event_num) + cluster(ccode), data = tc_dat,
                    method = "efron")

tc1_pwp_all <- coxph(tc1_gap ~ lib_dem_vdem + area_1000_log + lmtnest + elf +
                       rivals_thom + Devel + strata(event_num) +
                       cluster(ccode), data = tc1, method = "efron")

tc23_pwp_all <- coxph(tc23_gap ~ peace_wtc + fighting_wtc  +
                        force_lin_decay + favor_tc_lin_decay +
                        absorb_lin_decay + lib_dem_vdem + area_1000_log +
                        lmtnest + elf + rivals_thom + Devel + strata(event_num) +
                        cluster(ccode), data = tc23, method = "efron")

tc45_pwp_all <- coxph(tc45_gap ~ peace_wtc + fighting_wtc +
                        force_lin_decay + favor_tc_lin_decay +
                        absorb_lin_decay + lib_dem_vdem + area_1000_log +
                        lmtnest + elf + rivals_thom + Devel + strata(event_num) +
                        cluster(ccode), data = tc45, method = "efron")

tc6_more_pwp_all <- coxph(tc6_more_gap ~ peace_wtc + fighting_wtc  +
                            force_lin_decay + favor_tc_lin_decay +
                            absorb_lin_decay + lib_dem_vdem + area_1000_log +
                            lmtnest + elf + rivals_thom + Devel + strata(event_num) +
                            cluster(ccode), data = tc6_more, method = "efron")

stargazer(tc_pwp_all, tc1_pwp_all, tc23_pwp_all, tc45_pwp_all, tc6_more_pwp_all,
          type = "latex",
          title = "PWP Gap Time Model Results - Including Rivalry and Capacity",
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
                               "Number of Rivals",
                               "Capacity",
                               "TC Tally"),
          keep.stat = c("n"))


