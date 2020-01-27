# Clean Code #

# Packages, Settings, and Data #

library(foreign) # download .dta files for Stata < 13
library(haven) # download .dta files for Stata 13+
library(sampleSelection) # Package containing heckit for heckprob substitute
library(dplyr) #Used for dataset laterations and survival time calculations
library(lubridate) #Used for survival data formatting
library(survival)
library(stargazer)
library(reReg) #recurrent events package
library(gmodels)
library(Hmisc)
library(survminer)
library(visreg)



options(scipen = 50) # no scientific notation

# setwd("C:/Users/mjw65/Dropbox/State Making Strategies/State-Making-Strategies") #209 Pond Lab

TCDat <- read.csv("Survival-Data-TC.csv")
names(TCDat)

# Making Variables

#TCDat$BadEnd <- ifelse(TCDat$Forceful == 1 | TCDat$Absorbed == 1, 1, 0)

#TCDat$GoodEnd <- ifelse(TCDat$Peaceful == 1 | TCDat$Promoted == 1, 1, 0)

#TCDat$temp <- TCDat$Forceful_t10 + TCDat$Absorbed_t10
#TCDat$BadEnd_t10 <- ifelse(TCDat$temp >= 1,1,0)

#TCDat$temp <- TCDat$Forceful_t20 + TCDat$Absorbed_t20
#TCDat$BadEnd_t20 <- ifelse(TCDat$temp >= 1,1,0)

#TCDat$temp <- TCDat$Forceful_end + TCDat$Absorbed_end
#TCDat$BadEnd_tever <- ifelse(TCDat$temp >= 1,1,0)

#TCDat$temp <- TCDat$Peaceful_t10 + TCDat$Promoted_t10
#TCDat$GoodEnd_t10 <- ifelse(TCDat$temp >= 1,1,0)

#TCDat$temp <- TCDat$Peaceful_t20 + TCDat$Promoted_t20
#TCDat$GoodEnd_t20 <- ifelse(TCDat$temp >= 1,1,0)

#TCDat$temp <- TCDat$Peaceful_end + TCDat$Promoted_end
#TCDat$GoodEnd_tever <- ifelse(TCDat$temp >= 1,1,0)

#write.csv(TCDat, "Survival-Data-TC.csv")

# Summary Statistics & Crosstabs

CrossTable(TCDat$new_TC_dummy, TCDat$Absorbed_t20)

CrossTable(TCDat$new_TC_dummy, TCDat$Forceful_t20)

CrossTable(TCDat$new_TC_dummy, TCDat$Peaceful_t20)

CrossTable(TCDat$new_TC_dummy, TCDat$Promoted_t20)

CrossTable(TCDat$new_TC_dummy, TCDat$TC_persists_alt)

CrossTable(TCDat$new_TC_dummy, TCDat$Peace_wTC)

CrossTable(TCDat$new_TC_dummy, TCDat$Fighting_wTC)

CrossTable(TCDat$new_TC_dummy, TCDat$TempExisting)

CrossTable(TCDat$new_TC_dummy, TCDat$polity2)

plot(TCDat$new_TC_dummy, TCDat$area_1000) 

plot(TCDat$new_TC_dummy, TCDat$lmtnest)

plot(TCDat$new_TC_dummy, TCDat$ELF)

plot(TCDat$new_TC_dummy, TCDat$Country_Age)

# Creating subset datasets

TC2On <- TCDat[TCDat$event_num>=2,]

TC2or3 <- filter(TCDat, event_num == 2 | event_num == 3)

TC3or4 <- filter(TCDat, event_num == 3 | event_num == 4)

TC4or5 <- filter(TCDat, event_num == 4 | event_num == 5)

TC5or6 <- filter(TCDat, event_num == 5 | event_num == 6)

TC5orMore <- filter(TCDat, event_num >= 5)

TC6orMore <- filter(TCDat, event_num >= 6)

TC7orMore <- filter(TCDat, event_num >= 7)

TC1 <- filter(TCDat, event_num == 1)

TC2 <- filter(TCDat, event_num == 2)

TC3 <- filter(TCDat, event_num == 3)

TC4 <- filter(TCDat, event_num == 4)

TC5 <- filter(TCDat, event_num == 5)

TC6 <- filter(TCDat, event_num == 6)

TC7 <- filter(TCDat, event_num == 7)

TC8 <- filter(TCDat, event_num == 8)

TC9 <- filter(TCDat, event_num == 9)

TC10 <- filter(TCDat, event_num == 10)


# Creating all survival objects

TCDat.Elapsed<-Surv(TCDat$alt_start,TCDat$alt_stop,TCDat$new_TC_dummy)
TCDat.Gap<-Surv(TCDat$start,TCDat$stop,TCDat$new_TC_dummy)

TC2On.Elapsed<-Surv(TC2On$alt_start,TC2On$alt_stop,TC2On$new_TC_dummy)
TC2On.Gap<-Surv(TC2On$start,TC2On$stop,TC2On$new_TC_dummy)

TC1.Elapsed<-Surv(TC1$alt_start,TC1$alt_stop,TC1$new_TC_dummy)
TC1.Gap<-Surv(TC1$start,TC1$stop,TC1$new_TC_dummy)

TC2.Elapsed<-Surv(TC2$alt_start,TC2$alt_stop,TC2$new_TC_dummy)
TC2.Gap<-Surv(TC2$start,TC2$stop,TC2$new_TC_dummy)

TC2or3.Elapsed<-Surv(TC2or3$alt_start,TC2or3$alt_stop,TC2or3$new_TC_dummy)
TC2or3.Gap<-Surv(TC2or3$start,TC2or3$stop,TC2or3$new_TC_dummy)

TC3.Elapsed<-Surv(TC3$alt_start,TC3$alt_stop,TC3$new_TC_dummy)
TC3.Gap<-Surv(TC3$start,TC3$stop,TC3$new_TC_dummy)

TC3or4.Elapsed<-Surv(TC3or4$alt_start,TC3or4$alt_stop,TC3or4$new_TC_dummy)
TC3or4.Gap<-Surv(TC3or4$start,TC3or4$stop,TC3or4$new_TC_dummy)

TC4.Elapsed<-Surv(TC4$alt_start,TC4$alt_stop,TC4$new_TC_dummy)
TC4.Gap<-Surv(TC4$start,TC4$stop,TC4$new_TC_dummy)

TC4or5.Elapsed<-Surv(TC4or5$alt_start,TC4or5$alt_stop,TC4or5$new_TC_dummy)
TC4or5.Gap<-Surv(TC4or5$start,TC4or5$stop,TC4or5$new_TC_dummy)

TC5.Elapsed<-Surv(TC5$alt_start,TC5$alt_stop,TC5$new_TC_dummy)
TC5.Gap<-Surv(TC5$start,TC5$stop,TC5$new_TC_dummy)

TC5orMore.Elapsed<-Surv(TC5orMore$alt_start,TC5orMore$alt_stop,TC5orMore$new_TC_dummy)
TC5orMore.Gap<-Surv(TC5orMore$start,TC5orMore$stop,TC5orMore$new_TC_dummy)

TC5or6.Elapsed<-Surv(TC5or6$alt_start,TC5or6$alt_stop,TC5or6$new_TC_dummy)
TC5or6.Gap<-Surv(TC5or6$start,TC5or6$stop,TC5or6$new_TC_dummy)

TC6.Elapsed<-Surv(TC6$alt_start,TC6$alt_stop,TC6$new_TC_dummy)
TC6.Gap<-Surv(TC6$start,TC6$stop,TC6$new_TC_dummy)

TC6orMore.Elapsed<-Surv(TC6orMore$alt_start,TC6orMore$alt_stop,TC6orMore$new_TC_dummy)
TC6orMore.Gap<-Surv(TC6orMore$start,TC6orMore$stop,TC6orMore$new_TC_dummy)

TC7orMore.Elapsed<-Surv(TC7orMore$alt_start,TC7orMore$alt_stop,TC7orMore$new_TC_dummy)
TC7orMore.Gap<-Surv(TC7orMore$start,TC7orMore$stop,TC7orMore$new_TC_dummy)

TC7.Elapsed<-Surv(TC7$alt_start,TC7$alt_stop,TC7$new_TC_dummy)
TC7.Gap<-Surv(TC7$start,TC7$stop,TC7$new_TC_dummy)

TC8.Elapsed<-Surv(TC8$alt_start,TC8$alt_stop,TC8$new_TC_dummy)
TC8.Gap<-Surv(TC8$start,TC8$stop,TC8$new_TC_dummy)

TC9.Elapsed<-Surv(TC9$alt_start,TC9$alt_stop,TC9$new_TC_dummy)
TC9.Gap<-Surv(TC9$start,TC9$stop,TC9$new_TC_dummy)

TC10.Elapsed<-Surv(TC10$alt_start,TC10$alt_stop,TC10$new_TC_dummy)
TC10.Gap<-Surv(TC10$start,TC10$stop,TC10$new_TC_dummy)

# Anderson-Gill Model - Variance Correction

TC.AG.2On.1 <- coxph(TC2On.Gap ~ Peace_wTC + Fighting_wTC + GoodEnd_t10 + BadEnd_t10
                     +polity2 + area_1000 + lmtnest + ELF   
                    +cluster(ccode),data=TC2On, method="efron")
summary(TC.AG.2On.1)

TC.AG.2On.2 <- coxph(TC2On.Elapsed ~ Peace_wTC + Fighting_wTC + Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20
                     +polity2 + area_1000 + lmtnest + ELF   
                     +cluster(ccode),data=TC2On, method="efron")
summary(TC.AG.2On.2)

# Conditional PWP Model

TC.PWP.1 <-coxph(TCDat.Gap ~ Peace_wTC + GoodEnd_t10 + BadEnd_t10
                    +polity2 + area_1000 + lmtnest + ELF
                    +strata(event_num) + cluster(ccode),data=TCDat, method="efron")

summary(TC.PWP.1)

TC.PWP.2 <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC + GoodEnd_t20 + BadEnd_t20
                     +polity2 + area_1000 + lmtnest + ELF + Country_Age
                     +strata(event_num) + cluster(ccode),data=TCDat, method="efron")

summary(TC.PWP.2)

# Cox Gap Time with Gaussian Frailty

TC.GausFrail<-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC + GoodEnd_t20 + BadEnd_t20
             +polity2 + area_1000 + lmtnest + ELF + Country_Age
             +frailty.gaussian(ccode),data=TCDat)

summary(TC.GausFrail) 

# Cox Gap Time with Gamma Frailty

TC.GammaFrail<-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC + Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20
                    +polity2 + area_1000 + lmtnest + ELF + Country_Age
                    +frailty.gamma(ccode),data=TCDat)

summary(TC.GammaFrail) 
#################################################### Unified PWP - All
TC.1.PWP <- coxph(TC1.Gap ~ +polity2 + area_1000 + lmtnest + ELF 
                  +strata(event_num) + cluster(ccode),data=TC1, method="efron")

TCAll.PWP.t10 <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC + GoodEnd_t10 + BadEnd_t10
                               +polity2 + area_1000 + lmtnest + ELF 
                               +strata(event_num) + cluster(ccode),data=TCDat, method="efron") 

TCAll.PWP.t20 <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_t20 + BadEnd_t20
                               +polity2 + area_1000 + lmtnest + ELF 
                               +strata(event_num) + cluster(ccode),data=TCDat, method="efron") 

TCAll.PWP.tever <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_end + Forceful_end + Peaceful_end + Promoted_end
                                 +polity2 + area_1000 + lmtnest + ELF 
                                 +strata(event_num) + cluster(ccode),data=TCDat, method="efron") 

TCAll.PWP.LinDecay <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_LinDecay + BadEnd_LinDecay
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +strata(event_num) + cluster(ccode),data=TCDat, method="efron") 

TCAll.PWP.LinDecay.disag <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_LinDecay + Forceful_LinDecay + Peaceful_LinDecay + Promoted_LinDecay
                                          +polity2 + area_1000 + lmtnest + ELF 
                                          +strata(event_num) + cluster(ccode),data=TCDat, method="efron") 

TCAll.PWP.ExpDecay5 <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_5 + BadEnd_ExpDecay_5
                                     +polity2 + area_1000 + lmtnest + ELF 
                                     +strata(event_num) + cluster(ccode),data=TCDat, method="efron") 

TCAll.PWP.ExpDecay5.disag <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_5 + Forceful_ExpDecay_5 + Peaceful_ExpDecay_5 + Promoted_ExpDecay_5
                                           +polity2 + area_1000 + lmtnest + ELF 
                                           +strata(event_num) + cluster(ccode),data=TCDat, method="efron") 

TCAll.PWP.ExpDecay10 <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_10 + BadEnd_ExpDecay_10
                                      +polity2 + area_1000 + lmtnest + ELF 
                                      +strata(event_num) + cluster(ccode),data=TCDat, method="efron") 

TCAll.PWP.ExpDecay10.disag <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_10 + Forceful_ExpDecay_10 + Peaceful_ExpDecay_10 + Promoted_ExpDecay_10
                                            +polity2 + area_1000 + lmtnest + ELF 
                                            +strata(event_num) + cluster(ccode),data=TCDat, method="efron") 

TCAll.PWP.ExpDecay15 <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_15 + BadEnd_ExpDecay_15
                                      +polity2 + area_1000 + lmtnest + ELF 
                                      +strata(event_num) + cluster(ccode),data=TCDat, method="efron") 

TCAll.PWP.ExpDecay15.disag <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_15 + Forceful_ExpDecay_15 + Peaceful_ExpDecay_15 + Promoted_ExpDecay_15
                                            +polity2 + area_1000 + lmtnest + ELF 
                                            +strata(event_num) + cluster(ccode),data=TCDat, method="efron") 

TCAll.PWP.ExpDecay20 <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_20 + BadEnd_ExpDecay_20
                                      +polity2 + area_1000 + lmtnest + ELF 
                                      +strata(event_num) + cluster(ccode),data=TCDat, method="efron") 

TCAll.PWP.ExpDecay20.disag <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_20 + Forceful_ExpDecay_20 + Peaceful_ExpDecay_20 + Promoted_ExpDecay_20
                                            +polity2 + area_1000 + lmtnest + ELF 
                                            +strata(event_num) + cluster(ccode),data=TCDat, method="efron") 

summary(TCAll.PWP.t10)
summary(TCAll.PWP.t20)
summary(TCAll.PWP.LinDecay)
summary(TCAll.PWP.ExpDecay5)
summary(TCAll.PWP.ExpDecay10)
summary(TCAll.PWP.ExpDecay15)
summary(TCAll.PWP.ExpDecay20)

stargazer(TCAll.PWP.t10, TCAll.PWP.t20, TCAll.PWP.LinDecay, TCAll.PWP.ExpDecay5, TCAll.PWP.ExpDecay10, TCAll.PWP.ExpDecay15, TCAll.PWP.ExpDecay20, 
          covariate.labels = c("Peace w/ TC", "Fighting w/ TC", "Good - t10", 
                               "Bad - t10", "Good - t20", "Bad - t20", 
                               "Good - LinDec", "Bad - LinDec", "Good - Exp5",
                               "Bad - Exp5", "Good Exp 10", "Bad- Exp 10", 
                               "Good - Exp 15", "Bad - Exp 15", "Good - Exp20",
                               "Bad - Exp20", "Polity", "Area", "Mountain", "ELF"),
          dep.var.caption = "", type = "latex" ,title = "PWP Models - all variations")

#################################################### Unified PWP - 2 or 3


TC2or3.PWP.t10 <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC + GoodEnd_t10 + BadEnd_t10
                               +polity2 + area_1000 + lmtnest + ELF 
                               +strata(event_num) + cluster(ccode),data=TC2or3, method="efron") 

TC2or3.PWP.t20 <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_t20 + BadEnd_t20
                               +polity2 + area_1000 + lmtnest + ELF 
                               +strata(event_num) + cluster(ccode),data=TC2or3, method="efron") 

TC2or3.PWP.tever <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_end + Forceful_end + Peaceful_end + Promoted_end
                                 +polity2 + area_1000 + lmtnest + ELF 
                                 +strata(event_num) + cluster(ccode),data=TC2or3, method="efron") 

TC2or3.PWP.LinDecay <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_LinDecay + BadEnd_LinDecay
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +strata(event_num) + cluster(ccode),data=TC2or3, method="efron") 

TC2or3.PWP.LinDecay.disag <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_LinDecay + Forceful_LinDecay + Peaceful_LinDecay + Promoted_LinDecay
                                          +polity2 + area_1000 + lmtnest + ELF 
                                          +strata(event_num) + cluster(ccode),data=TC2or3, method="efron") 

TC2or3.PWP.ExpDecay5 <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_5 + BadEnd_ExpDecay_5
                                     +polity2 + area_1000 + lmtnest + ELF 
                                     +strata(event_num) + cluster(ccode),data=TC2or3, method="efron") 

TC2or3.PWP.ExpDecay5.disag <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_5 + Forceful_ExpDecay_5 + Peaceful_ExpDecay_5 + Promoted_ExpDecay_5
                                           +polity2 + area_1000 + lmtnest + ELF 
                                           +strata(event_num) + cluster(ccode),data=TC2or3, method="efron") 

TC2or3.PWP.ExpDecay10 <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_10 + BadEnd_ExpDecay_10
                                      +polity2 + area_1000 + lmtnest + ELF 
                                      +strata(event_num) + cluster(ccode),data=TC2or3, method="efron") 

TC2or3.PWP.ExpDecay10.disag <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_10 + Forceful_ExpDecay_10 + Peaceful_ExpDecay_10 + Promoted_ExpDecay_10
                                            +polity2 + area_1000 + lmtnest + ELF 
                                            +strata(event_num) + cluster(ccode),data=TC2or3, method="efron") 

TC2or3.PWP.ExpDecay15 <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_15 + BadEnd_ExpDecay_15
                                      +polity2 + area_1000 + lmtnest + ELF 
                                      +strata(event_num) + cluster(ccode),data=TC2or3, method="efron") 

TC2or3.PWP.ExpDecay15.disag <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_15 + Forceful_ExpDecay_15 + Peaceful_ExpDecay_15 + Promoted_ExpDecay_15
                                            +polity2 + area_1000 + lmtnest + ELF 
                                            +strata(event_num) + cluster(ccode),data=TC2or3, method="efron") 

TC2or3.PWP.ExpDecay20 <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_20 + BadEnd_ExpDecay_20
                                      +polity2 + area_1000 + lmtnest + ELF 
                                      +strata(event_num) + cluster(ccode),data=TC2or3, method="efron") 

TC2or3.PWP.ExpDecay20.disag <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_20 + Forceful_ExpDecay_20 + Peaceful_ExpDecay_20 + Promoted_ExpDecay_20
                                            +polity2 + area_1000 + lmtnest + ELF 
                                            +strata(event_num) + cluster(ccode),data=TC2or3, method="efron") 

summary(TC2or3.PWP.LinDecay)

#################################################### Unified PWP - 4 or 5


TC4or5.PWP.t10 <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC + GoodEnd_t10 + BadEnd_t10
                               +polity2 + area_1000 + lmtnest + ELF 
                               +strata(event_num) + cluster(ccode),data=TC4or5, method="efron") 

TC4or5.PWP.t20 <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_t20 + BadEnd_t20
                               +polity2 + area_1000 + lmtnest + ELF 
                               +strata(event_num) + cluster(ccode),data=TC4or5, method="efron") 

TC4or5.PWP.tever <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_end + Forceful_end + Peaceful_end + Promoted_end
                                 +polity2 + area_1000 + lmtnest + ELF 
                                 +strata(event_num) + cluster(ccode),data=TC4or5, method="efron") 

TC4or5.PWP.LinDecay <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_LinDecay + BadEnd_LinDecay
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +strata(event_num) + cluster(ccode),data=TC4or5, method="efron") 

TC4or5.PWP.LinDecay.disag <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_LinDecay + Forceful_LinDecay + Peaceful_LinDecay + Promoted_LinDecay
                                          +polity2 + area_1000 + lmtnest + ELF 
                                          +strata(event_num) + cluster(ccode),data=TC4or5, method="efron") 

TC4or5.PWP.ExpDecay5 <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_5 + BadEnd_ExpDecay_5
                                     +polity2 + area_1000 + lmtnest + ELF 
                                     +strata(event_num) + cluster(ccode),data=TC4or5, method="efron") 

TC4or5.PWP.ExpDecay5.disag <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_5 + Forceful_ExpDecay_5 + Peaceful_ExpDecay_5 + Promoted_ExpDecay_5
                                           +polity2 + area_1000 + lmtnest + ELF 
                                           +strata(event_num) + cluster(ccode),data=TC4or5, method="efron") 

TC4or5.PWP.ExpDecay10 <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_10 + BadEnd_ExpDecay_10
                                      +polity2 + area_1000 + lmtnest + ELF 
                                      +strata(event_num) + cluster(ccode),data=TC4or5, method="efron") 

TC4or5.PWP.ExpDecay10.disag <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_10 + Forceful_ExpDecay_10 + Peaceful_ExpDecay_10 + Promoted_ExpDecay_10
                                            +polity2 + area_1000 + lmtnest + ELF 
                                            +strata(event_num) + cluster(ccode),data=TC4or5, method="efron") 

TC4or5.PWP.ExpDecay15 <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_15 + BadEnd_ExpDecay_15
                                      +polity2 + area_1000 + lmtnest + ELF 
                                      +strata(event_num) + cluster(ccode),data=TC4or5, method="efron") 

TC4or5.PWP.ExpDecay15.disag <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_15 + Forceful_ExpDecay_15 + Peaceful_ExpDecay_15 + Promoted_ExpDecay_15
                                            +polity2 + area_1000 + lmtnest + ELF 
                                            +strata(event_num) + cluster(ccode),data=TC4or5, method="efron") 

TC4or5.PWP.ExpDecay20 <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_20 + BadEnd_ExpDecay_20
                                      +polity2 + area_1000 + lmtnest + ELF 
                                      +strata(event_num) + cluster(ccode),data=TC4or5, method="efron") 

TC4or5.PWP.ExpDecay20.disag <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_20 + Forceful_ExpDecay_20 + Peaceful_ExpDecay_20 + Promoted_ExpDecay_20
                                            +polity2 + area_1000 + lmtnest + ELF 
                                            +strata(event_num) + cluster(ccode),data=TC4or5, method="efron") 

#################################################### Unified PWP - 6 or more


TC6orMore.PWP.t10 <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC + GoodEnd_t10 + BadEnd_t10
                               +polity2 + area_1000 + lmtnest + ELF 
                               +strata(event_num) + cluster(ccode),data=TC6orMore, method="efron") 

TC6orMore.PWP.t20 <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_t20 + BadEnd_t20
                               +polity2 + area_1000 + lmtnest + ELF 
                               +strata(event_num) + cluster(ccode),data=TC6orMore, method="efron") 

TC6orMore.PWP.tever <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_end + Forceful_end + Peaceful_end + Promoted_end
                                 +polity2 + area_1000 + lmtnest + ELF 
                                 +strata(event_num) + cluster(ccode),data=TC6orMore, method="efron") 

TC6orMore.PWP.LinDecay <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_LinDecay + BadEnd_LinDecay
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +strata(event_num) + cluster(ccode),data=TC6orMore, method="efron") 

TC6orMore.PWP.LinDecay.disag <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_LinDecay + Forceful_LinDecay + Peaceful_LinDecay + Promoted_LinDecay
                                          +polity2 + area_1000 + lmtnest + ELF 
                                          +strata(event_num) + cluster(ccode),data=TC6orMore, method="efron") 

TC6orMore.PWP.ExpDecay5 <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_5 + BadEnd_ExpDecay_5
                                     +polity2 + area_1000 + lmtnest + ELF 
                                     +strata(event_num) + cluster(ccode),data=TC6orMore, method="efron") 

TC6orMore.PWP.ExpDecay5.disag <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_5 + Forceful_ExpDecay_5 + Peaceful_ExpDecay_5 + Promoted_ExpDecay_5
                                           +polity2 + area_1000 + lmtnest + ELF 
                                           +strata(event_num) + cluster(ccode),data=TC6orMore, method="efron") 

TC6orMore.PWP.ExpDecay10 <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_10 + BadEnd_ExpDecay_10
                                      +polity2 + area_1000 + lmtnest + ELF 
                                      +strata(event_num) + cluster(ccode),data=TC6orMore, method="efron") 

TC6orMore.PWP.ExpDecay10.disag <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_10 + Forceful_ExpDecay_10 + Peaceful_ExpDecay_10 + Promoted_ExpDecay_10
                                            +polity2 + area_1000 + lmtnest + ELF 
                                            +strata(event_num) + cluster(ccode),data=TC6orMore, method="efron") 

TC6orMore.PWP.ExpDecay15 <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_15 + BadEnd_ExpDecay_15
                                      +polity2 + area_1000 + lmtnest + ELF 
                                      +strata(event_num) + cluster(ccode),data=TC6orMore, method="efron") 

TC6orMore.PWP.ExpDecay15.disag <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_15 + Forceful_ExpDecay_15 + Peaceful_ExpDecay_15 + Promoted_ExpDecay_15
                                            +polity2 + area_1000 + lmtnest + ELF 
                                            +strata(event_num) + cluster(ccode),data=TC6orMore, method="efron") 

TC6orMore.PWP.ExpDecay20 <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_20 + BadEnd_ExpDecay_20
                                      +polity2 + area_1000 + lmtnest + ELF 
                                      +strata(event_num) + cluster(ccode),data=TC6orMore, method="efron") 

TC6orMore.PWP.ExpDecay20.disag <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_20 + Forceful_ExpDecay_20 + Peaceful_ExpDecay_20 + Promoted_ExpDecay_20
                                            +polity2 + area_1000 + lmtnest + ELF 
                                          +strata(event_num) + cluster(ccode),data=TC6orMore, method="efron") 


########################################################### Gaussian Frailty - all

TCAll.Frailty.Gaus.t10 <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC + GoodEnd_t10 + BadEnd_t10
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gaussian(ccode),data=TCDat) 

TCAll.Frailty.Gaus.t20 <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_t20 + BadEnd_t20
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gaussian(ccode),data=TCDat) 

TCAll.Frailty.Gaus.tever <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_end + Forceful_end + Peaceful_end + Promoted_end
                               +polity2 + area_1000 + lmtnest + ELF 
                               +frailty.gaussian(ccode),data=TCDat) 

TCAll.Frailty.Gaus.LinDecay <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_LinDecay + BadEnd_LinDecay
                                  +polity2 + area_1000 + lmtnest + ELF 
                                  +frailty.gaussian(ccode),data=TCDat) 

TCAll.Frailty.Gaus.LinDecay.disag <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_LinDecay + Forceful_LinDecay + Peaceful_LinDecay + Promoted_LinDecay
                                        +polity2 + area_1000 + lmtnest + ELF 
                                        +frailty.gaussian(ccode),data=TCDat) 

TCAll.Frailty.Gaus.ExpDecay5 <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_5 + BadEnd_ExpDecay_5
                                   +polity2 + area_1000 + lmtnest + ELF 
                                   +frailty.gaussian(ccode),data=TCDat) 

TCAll.Frailty.Gaus.ExpDecay5.disag <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_5 + Forceful_ExpDecay_5 + Peaceful_ExpDecay_5 + Promoted_ExpDecay_5
                                         +polity2 + area_1000 + lmtnest + ELF 
                                         +frailty.gaussian(ccode),data=TCDat) 

TCAll.Frailty.Gaus.ExpDecay10 <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_10 + BadEnd_ExpDecay_10
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +frailty.gaussian(ccode),data=TCDat) 

TCAll.Frailty.Gaus.ExpDecay10.disag <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_10 + Forceful_ExpDecay_10 + Peaceful_ExpDecay_10 + Promoted_ExpDecay_10
                                          +polity2 + area_1000 + lmtnest + ELF 
                                          +frailty.gaussian(ccode),data=TCDat) 

TCAll.Frailty.Gaus.ExpDecay15 <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_15 + BadEnd_ExpDecay_15
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +frailty.gaussian(ccode),data=TCDat) 

TCAll.Frailty.Gaus.ExpDecay15.disag <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_15 + Forceful_ExpDecay_15 + Peaceful_ExpDecay_15 + Promoted_ExpDecay_15
                                          +polity2 + area_1000 + lmtnest + ELF 
                                          +frailty.gaussian(ccode),data=TCDat) 

TCAll.Frailty.Gaus.ExpDecay20 <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_20 + BadEnd_ExpDecay_20
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +frailty.gaussian(ccode),data=TCDat) 

TCAll.Frailty.Gaus.ExpDecay20.disag <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_20 + Forceful_ExpDecay_20 + Peaceful_ExpDecay_20 + Promoted_ExpDecay_20
                                          +polity2 + area_1000 + lmtnest + ELF 
                                          +frailty.gaussian(ccode),data=TCDat) 


############################################################################################

# Stratified models for gaussian frailty

TC1.Frailty.Gaus <-coxph(TC1.Gap ~ polity2 + area_1000 + lmtnest + ELF 
                    +frailty.gaussian(ccode),data=TC1) 

###########################################################

TC2.Frailty.Gaus.t10 <-coxph(TC2.Gap ~ Peace_wTC + Fighting_wTC + GoodEnd_t10 + BadEnd_t10
                    +polity2 + area_1000 + lmtnest + ELF 
                    +frailty.gaussian(ccode),data=TC2) 

TC2.Frailty.Gaus.t20 <-coxph(TC2.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_t20 + BadEnd_t20
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gaussian(ccode),data=TC2) 

TC2.Frailty.Gaus.tever <-coxph(TC2.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_end + Forceful_end + Peaceful_end + Promoted_end
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gaussian(ccode),data=TC2) 

TC2.Frailty.Gaus.LinDecay <-coxph(TC2.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_LinDecay + BadEnd_LinDecay
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gaussian(ccode),data=TC2) 

TC2.Frailty.Gaus.LinDecay.disag <-coxph(TC2.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_LinDecay + Forceful_LinDecay + Peaceful_LinDecay + Promoted_LinDecay
                                  +polity2 + area_1000 + lmtnest + ELF 
                                  +frailty.gaussian(ccode),data=TC2) 

TC2.Frailty.Gaus.ExpDecay5 <-coxph(TC2.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_5 + BadEnd_ExpDecay_5
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gaussian(ccode),data=TC2) 

TC2.Frailty.Gaus.ExpDecay5.disag <-coxph(TC2.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_5 + Forceful_ExpDecay_5 + Peaceful_ExpDecay_5 + Promoted_ExpDecay_5
                                   +polity2 + area_1000 + lmtnest + ELF 
                                   +frailty.gaussian(ccode),data=TC2) 

TC2.Frailty.Gaus.ExpDecay10 <-coxph(TC2.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_10 + BadEnd_ExpDecay_10
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gaussian(ccode),data=TC2) 

TC2.Frailty.Gaus.ExpDecay10.disag <-coxph(TC2.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_10 + Forceful_ExpDecay_10 + Peaceful_ExpDecay_10 + Promoted_ExpDecay_10
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +frailty.gaussian(ccode),data=TC2) 

TC2.Frailty.Gaus.ExpDecay15 <-coxph(TC2.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_15 + BadEnd_ExpDecay_15
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gaussian(ccode),data=TC2) 

TC2.Frailty.Gaus.ExpDecay15.disag <-coxph(TC2.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_15 + Forceful_ExpDecay_15 + Peaceful_ExpDecay_15 + Promoted_ExpDecay_15
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +frailty.gaussian(ccode),data=TC2) 

TC2.Frailty.Gaus.ExpDecay20 <-coxph(TC2.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_20 + BadEnd_ExpDecay_20
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +frailty.gaussian(ccode),data=TC2) 

TC2.Frailty.Gaus.ExpDecay20.disag <-coxph(TC2.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_20 + Forceful_ExpDecay_20 + Peaceful_ExpDecay_20 + Promoted_ExpDecay_20
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +frailty.gaussian(ccode),data=TC2) 

##############################################

TC2or3.Frailty.Gaus.t10 <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC   + GoodEnd_t10 + BadEnd_t10
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gaussian(ccode),data=TC2or3) 

TC2or3.Frailty.Gaus.t20 <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_t20 + BadEnd_t20
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gaussian(ccode),data=TC2or3) 

TC2or3.Frailty.Gaus.tever <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_end + Forceful_end + Peaceful_end + Promoted_end
                               +polity2 + area_1000 + lmtnest + ELF 
                               +frailty.gaussian(ccode),data=TC2or3)

TC2or3.Frailty.Gaus.LinDecay <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC + GoodEnd_LinDecay + BadEnd_LinDecay
                                  +polity2 + area_1000 + lmtnest + ELF 
                                  +frailty.gaussian(ccode),data=TC2or3)

TC2or3.Frailty.Gaus.LinDecay.disag <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_LinDecay + Forceful_LinDecay + Peaceful_LinDecay + Promoted_LinDecay
                                           +polity2 + area_1000 + lmtnest + ELF 
                                           +frailty.gaussian(ccode),data=TC2or3)

TC2or3.Frailty.Gaus.ExpDecay5 <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_5 + BadEnd_ExpDecay_5
                                   +polity2 + area_1000 + lmtnest + ELF 
                                   +frailty.gaussian(ccode),data=TC2or3) 

TC2or3.Frailty.Gaus.ExpDecay5.disag <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_5 + Forceful_ExpDecay_5 + Peaceful_ExpDecay_5 + Promoted_ExpDecay_5
                                            +polity2 + area_1000 + lmtnest + ELF 
                                            +frailty.gaussian(ccode),data=TC2or3) 

TC2or3.Frailty.Gaus.ExpDecay10 <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_10 + BadEnd_ExpDecay_10
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +frailty.gaussian(ccode),data=TC2or3) 

TC2or3.Frailty.Gaus.ExpDecay10.disag <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_10 + Forceful_ExpDecay_10 + Peaceful_ExpDecay_10 + Promoted_ExpDecay_10
                                             +polity2 + area_1000 + lmtnest + ELF 
                                             +frailty.gaussian(ccode),data=TC2or3) 

TC2or3.Frailty.Gaus.ExpDecay15 <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_15 + BadEnd_ExpDecay_15
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +frailty.gaussian(ccode),data=TC2or3) 

TC2or3.Frailty.Gaus.ExpDecay15.disag <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_15 + Forceful_ExpDecay_15 + Peaceful_ExpDecay_15 + Promoted_ExpDecay_15
                                             +polity2 + area_1000 + lmtnest + ELF 
                                             +frailty.gaussian(ccode),data=TC2or3) 

TC2or3.Frailty.Gaus.ExpDecay20 <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_20 + BadEnd_ExpDecay_20
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +frailty.gaussian(ccode),data=TC2or3) 

TC2or3.Frailty.Gaus.ExpDecay20.disag <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_20 + Forceful_ExpDecay_20 + Peaceful_ExpDecay_20 + Promoted_ExpDecay_20
                                             +polity2 + area_1000 + lmtnest + ELF 
                                             +frailty.gaussian(ccode),data=TC2or3) 

##############################################

TC3.Frailty.Gaus.t10 <-coxph(TC3.Gap ~ Peace_wTC + Fighting_wTC   + GoodEnd_t10 + BadEnd_t10
                    +polity2 + area_1000 + lmtnest + ELF 
                    +frailty.gaussian(ccode),data=TC3) 

TC3.Frailty.Gaus.t20 <-coxph(TC3.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_t20 + BadEnd_t20
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gaussian(ccode),data=TC3) 

TC3.Frailty.Gaus.tever <-coxph(TC3.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_end + Forceful_end + Peaceful_end + Promoted_end
                               +polity2 + area_1000 + lmtnest + ELF 
                               +frailty.gaussian(ccode),data=TC3) 

TC3.Frailty.Gaus.LinDecay <-coxph(TC3.Gap ~ Peace_wTC + Fighting_wTC + GoodEnd_LinDecay + BadEnd_LinDecay
                                  +polity2 + area_1000 + lmtnest + ELF 
                                  +frailty.gaussian(ccode),data=TC3) 

TC3.Frailty.Gaus.LinDecay.disag <-coxph(TC3.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_LinDecay + Forceful_LinDecay + Peaceful_LinDecay + Promoted_LinDecay
                                  +polity2 + area_1000 + lmtnest + ELF 
                                  +frailty.gaussian(ccode),data=TC3)

TC3.Frailty.Gaus.ExpDecay5 <-coxph(TC3.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_5 + BadEnd_ExpDecay_5
                                   +polity2 + area_1000 + lmtnest + ELF 
                                   +frailty.gaussian(ccode),data=TC3)

TC3.Frailty.Gaus.ExpDecay5.disag <-coxph(TC3.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_5 + Forceful_ExpDecay_5 + Peaceful_ExpDecay_5 + Promoted_ExpDecay_5
                                         +polity2 + area_1000 + lmtnest + ELF 
                                         +frailty.gaussian(ccode),data=TC3) 

TC3.Frailty.Gaus.ExpDecay10 <-coxph(TC3.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_10 + BadEnd_ExpDecay_10
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +frailty.gaussian(ccode),data=TC3) 

TC3.Frailty.Gaus.ExpDecay10.disag <-coxph(TC3.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_10 + Forceful_ExpDecay_10 + Peaceful_ExpDecay_10 + Promoted_ExpDecay_10
                                          +polity2 + area_1000 + lmtnest + ELF 
                                          +frailty.gaussian(ccode),data=TC3) 

TC3.Frailty.Gaus.ExpDecay15 <-coxph(TC3.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_15 + BadEnd_ExpDecay_15
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +frailty.gaussian(ccode),data=TC3) 

TC3.Frailty.Gaus.ExpDecay15.disag <-coxph(TC3.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_15 + Forceful_ExpDecay_15 + Peaceful_ExpDecay_15 + Promoted_ExpDecay_15
                                          +polity2 + area_1000 + lmtnest + ELF 
                                          +frailty.gaussian(ccode),data=TC3) 

TC3.Frailty.Gaus.ExpDecay20 <-coxph(TC3.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_20 + BadEnd_ExpDecay_20
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +frailty.gaussian(ccode),data=TC3) 

TC3.Frailty.Gaus.ExpDecay20.disag <-coxph(TC3.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_20 + Forceful_ExpDecay_20 + Peaceful_ExpDecay_20 + Promoted_ExpDecay_20
                                          +polity2 + area_1000 + lmtnest + ELF 
                                          +frailty.gaussian(ccode),data=TC3) 

##############################################

TC3or4.Frailty.Gaus.t10 <-coxph(TC3or4.Gap ~ Peace_wTC + Fighting_wTC   + GoodEnd_t10 + BadEnd_t10
                                +polity2 + area_1000 + lmtnest + ELF 
                                +frailty.gaussian(ccode),data=TC3or4) 

TC3or4.Frailty.Gaus.t20 <-coxph(TC3or4.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_t20 + BadEnd_t20
                                +polity2 + area_1000 + lmtnest + ELF 
                                +frailty.gaussian(ccode),data=TC3or4) 

TC3or4.Frailty.Gaus.tever <-coxph(TC3or4.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_end + Forceful_end + Peaceful_end + Promoted_end
                                  +polity2 + area_1000 + lmtnest + ELF 
                                  +frailty.gaussian(ccode),data=TC3or4) 

TC3or4.Frailty.Gaus.LinDecay <-coxph(TC3or4.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_LinDecay + BadEnd_LinDecay
                                     +polity2 + area_1000 + lmtnest + ELF 
                                     +frailty.gaussian(ccode),data=TC3or4) 

TC3or4.Frailty.Gaus.LinDecay.disag <-coxph(TC3or4.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_LinDecay + Forceful_LinDecay + Peaceful_LinDecay + Promoted_LinDecay
                                     +polity2 + area_1000 + lmtnest + ELF 
                                     +frailty.gaussian(ccode),data=TC3or4) 

TC3or4.Frailty.Gaus.ExpDecay5 <-coxph(TC3or4.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_5 + BadEnd_ExpDecay_5
                                      +polity2 + area_1000 + lmtnest + ELF 
                                      +frailty.gaussian(ccode),data=TC3or4) 

TC3or4.Frailty.Gaus.ExpDecay5.disag <-coxph(TC3or4.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_5 + Forceful_ExpDecay_5 + Peaceful_ExpDecay_5 + Promoted_ExpDecay_5
                                      +polity2 + area_1000 + lmtnest + ELF 
                                      +frailty.gaussian(ccode),data=TC3or4) 

TC3or4.Frailty.Gaus.ExpDecay10 <-coxph(TC3or4.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_10 + BadEnd_ExpDecay_10
                                       +polity2 + area_1000 + lmtnest + ELF 
                                       +frailty.gaussian(ccode),data=TC3or4) 

TC3or4.Frailty.Gaus.ExpDecay10.disag <-coxph(TC3or4.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_10 + Forceful_ExpDecay_10 + Peaceful_ExpDecay_10 + Promoted_ExpDecay_10
                                       +polity2 + area_1000 + lmtnest + ELF 
                                       +frailty.gaussian(ccode),data=TC3or4) 

TC3or4.Frailty.Gaus.ExpDecay15 <-coxph(TC3or4.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_15 + BadEnd_ExpDecay_15
                                       +polity2 + area_1000 + lmtnest + ELF 
                                       +frailty.gaussian(ccode),data=TC3or4) 

TC3or4.Frailty.Gaus.ExpDecay15.disag <-coxph(TC3or4.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_15 + Forceful_ExpDecay_15 + Peaceful_ExpDecay_15 + Promoted_ExpDecay_15
                                       +polity2 + area_1000 + lmtnest + ELF 
                                       +frailty.gaussian(ccode),data=TC3or4) 

TC3or4.Frailty.Gaus.ExpDecay20 <-coxph(TC3or4.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_20 + BadEnd_ExpDecay_20
                                       +polity2 + area_1000 + lmtnest + ELF 
                                       +frailty.gaussian(ccode),data=TC3or4) 

TC3or4.Frailty.Gaus.ExpDecay20.disag <-coxph(TC3or4.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_20 + Forceful_ExpDecay_20 + Peaceful_ExpDecay_20 + Promoted_ExpDecay_20
                                       +polity2 + area_1000 + lmtnest + ELF 
                                       +frailty.gaussian(ccode),data=TC3or4) 

##############################################

TC4.Frailty.Gaus.t10 <-coxph(TC4.Gap ~ Peace_wTC + Fighting_wTC   + GoodEnd_t10 + BadEnd_t10
                    +polity2 + area_1000 + lmtnest + ELF 
                    +frailty.gaussian(ccode),data=TC4) 

TC4.Frailty.Gaus.t20 <-coxph(TC4.Gap ~ Peace_wTC + Fighting_wTC   + GoodEnd_t20 + BadEnd_t20
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gaussian(ccode),data=TC4)

TC4.Frailty.Gaus.tever <-coxph(TC4.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_end + Forceful_end + Peaceful_end + Promoted_end
                               +polity2 + area_1000 + lmtnest + ELF 
                               +frailty.gaussian(ccode),data=TC4) 

TC4.Frailty.Gaus.LinDecay <-coxph(TC4.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_LinDecay + BadEnd_LinDecay
                                  +polity2 + area_1000 + lmtnest + ELF 
                                  +frailty.gaussian(ccode),data=TC4) 

TC4.Frailty.Gaus.LinDecay.disag <-coxph(TC4.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_LinDecay + Forceful_LinDecay + Peaceful_LinDecay + Promoted_LinDecay
                                        +polity2 + area_1000 + lmtnest + ELF 
                                        +frailty.gaussian(ccode),data=TC4) 

TC4.Frailty.Gaus.ExpDecay5 <-coxph(TC4.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_5 + BadEnd_ExpDecay_5
                                   +polity2 + area_1000 + lmtnest + ELF 
                                   +frailty.gaussian(ccode),data=TC4) 

TC4.Frailty.Gaus.ExpDecay5.disag <-coxph(TC4.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_5 + Forceful_ExpDecay_5 + Peaceful_ExpDecay_5 + Promoted_ExpDecay_5
                                         +polity2 + area_1000 + lmtnest + ELF 
                                         +frailty.gaussian(ccode),data=TC4) 

TC4.Frailty.Gaus.ExpDecay10 <-coxph(TC4.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_10 + BadEnd_ExpDecay_10
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +frailty.gaussian(ccode),data=TC4) 

TC4.Frailty.Gaus.ExpDecay10.disag <-coxph(TC4.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_10 + Forceful_ExpDecay_10 + Peaceful_ExpDecay_10 + Promoted_ExpDecay_10
                                          +polity2 + area_1000 + lmtnest + ELF 
                                          +frailty.gaussian(ccode),data=TC4) 

TC4.Frailty.Gaus.ExpDecay15 <-coxph(TC4.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_15 + BadEnd_ExpDecay_15
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +frailty.gaussian(ccode),data=TC4) 

TC4.Frailty.Gaus.ExpDecay15.disag <-coxph(TC4.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_15 + Forceful_ExpDecay_15 + Peaceful_ExpDecay_15 + Promoted_ExpDecay_15
                                          +polity2 + area_1000 + lmtnest + ELF 
                                          +frailty.gaussian(ccode),data=TC4) 

TC4.Frailty.Gaus.ExpDecay20 <-coxph(TC4.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_20 + BadEnd_ExpDecay_20
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +frailty.gaussian(ccode),data=TC4) 

TC4.Frailty.Gaus.ExpDecay20.disag <-coxph(TC4.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_20 + Forceful_ExpDecay_20 + Peaceful_ExpDecay_20 + Promoted_ExpDecay_20
                                          +polity2 + area_1000 + lmtnest + ELF 
                                          +frailty.gaussian(ccode),data=TC4) 

##############################################

TC4or5.Frailty.Gaus.t10 <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC   + GoodEnd_t10 + BadEnd_t10
                                +polity2 + area_1000 + lmtnest + ELF 
                                +frailty.gaussian(ccode),data=TC4or5) 

TC4or5.Frailty.Gaus.t20 <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_t20 + BadEnd_t20
                                +polity2 + area_1000 + lmtnest + ELF 
                                +frailty.gaussian(ccode),data=TC4or5) 

TC4or5.Frailty.Gaus.tever <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_end + Forceful_end + Peaceful_end + Promoted_end
                                  +polity2 + area_1000 + lmtnest + ELF 
                                  +frailty.gaussian(ccode),data=TC4or5) 

TC4or5.Frailty.Gaus.LinDecay <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_LinDecay + BadEnd_LinDecay
                                     +polity2 + area_1000 + lmtnest + ELF 
                                     +frailty.gaussian(ccode),data=TC4or5)

TC4or5.Frailty.Gaus.LinDecay.disag <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_LinDecay + Forceful_LinDecay + Peaceful_LinDecay + Promoted_LinDecay
                                           +polity2 + area_1000 + lmtnest + ELF 
                                           +frailty.gaussian(ccode),data=TC4or5)

TC4or5.Frailty.Gaus.ExpDecay5 <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_5 + BadEnd_ExpDecay_5
                                      +polity2 + area_1000 + lmtnest + ELF 
                                      +frailty.gaussian(ccode),data=TC4or5) 

TC4or5.Frailty.Gaus.ExpDecay5.disag <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_5 + Forceful_ExpDecay_5 + Peaceful_ExpDecay_5 + Promoted_ExpDecay_5
                                            +polity2 + area_1000 + lmtnest + ELF 
                                            +frailty.gaussian(ccode),data=TC4or5)

TC4or5.Frailty.Gaus.ExpDecay10 <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_10 + BadEnd_ExpDecay_10
                                       +polity2 + area_1000 + lmtnest + ELF 
                                       +frailty.gaussian(ccode),data=TC4or5) 

TC4or5.Frailty.Gaus.ExpDecay10.disag <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_10 + Forceful_ExpDecay_10 + Peaceful_ExpDecay_10 + Promoted_ExpDecay_10
                                             +polity2 + area_1000 + lmtnest + ELF 
                                             +frailty.gaussian(ccode),data=TC4or5)

TC4or5.Frailty.Gaus.ExpDecay15 <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_15 + BadEnd_ExpDecay_15
                                       +polity2 + area_1000 + lmtnest + ELF 
                                       +frailty.gaussian(ccode),data=TC4or5) 

TC4or5.Frailty.Gaus.ExpDecay15.disag <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_15 + Forceful_ExpDecay_15 + Peaceful_ExpDecay_15 + Promoted_ExpDecay_15
                                             +polity2 + area_1000 + lmtnest + ELF 
                                             +frailty.gaussian(ccode),data=TC4or5)

TC4or5.Frailty.Gaus.ExpDecay20 <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_20 + BadEnd_ExpDecay_20
                                       +polity2 + area_1000 + lmtnest + ELF 
                                       +frailty.gaussian(ccode),data=TC4or5) 

TC4or5.Frailty.Gaus.ExpDecay20.disag <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_20 + Forceful_ExpDecay_20 + Peaceful_ExpDecay_20 + Promoted_ExpDecay_20
                                             +polity2 + area_1000 + lmtnest + ELF 
                                             +frailty.gaussian(ccode),data=TC4or5) 

##############################################

TC5.Frailty.Gaus.t10 <-coxph(TC5.Gap ~ Peace_wTC + Fighting_wTC   + GoodEnd_t10 + BadEnd_t10
                    +polity2 + area_1000 + lmtnest + ELF 
                    +frailty.gaussian(ccode),data=TC5) 

TC5.Frailty.Gaus.t20 <-coxph(TC5.Gap ~ Peace_wTC + Fighting_wTC   + GoodEnd_t20 + BadEnd_t20
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gaussian(ccode),data=TC5) 

TC5.Frailty.Gaus.tever <-coxph(TC5.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_end + Forceful_end + Peaceful_end + Promoted_end
                               +polity2 + area_1000 + lmtnest + ELF 
                               +frailty.gaussian(ccode),data=TC5) 

TC5.Frailty.Gaus.LinDecay <-coxph(TC5.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_LinDecay + BadEnd_LinDecay
                                  +polity2 + area_1000 + lmtnest + ELF 
                                  +frailty.gaussian(ccode),data=TC5) 

TC5.Frailty.Gaus.LinDecay.disag <-coxph(TC5.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_LinDecay + Forceful_LinDecay + Peaceful_LinDecay + Promoted_LinDecay
                                        +polity2 + area_1000 + lmtnest + ELF 
                                        +frailty.gaussian(ccode),data=TC5) 

TC5.Frailty.Gaus.ExpDecay5 <-coxph(TC5.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_5 + BadEnd_ExpDecay_5
                                   +polity2 + area_1000 + lmtnest + ELF 
                                   +frailty.gaussian(ccode),data=TC5) 

TC5.Frailty.Gaus.ExpDecay5.disag <-coxph(TC5.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_5 + Forceful_ExpDecay_5 + Peaceful_ExpDecay_5 + Promoted_ExpDecay_5
                                         +polity2 + area_1000 + lmtnest + ELF 
                                         +frailty.gaussian(ccode),data=TC5) 

TC5.Frailty.Gaus.ExpDecay10 <-coxph(TC5.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_10 + BadEnd_ExpDecay_10
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +frailty.gaussian(ccode),data=TC5) 

TC5.Frailty.Gaus.ExpDecay10.disag <-coxph(TC5.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_10 + Forceful_ExpDecay_10 + Peaceful_ExpDecay_10 + Promoted_ExpDecay_10
                                          +polity2 + area_1000 + lmtnest + ELF 
                                          +frailty.gaussian(ccode),data=TC5) 

TC5.Frailty.Gaus.ExpDecay15 <-coxph(TC5.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_15 + BadEnd_ExpDecay_15
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +frailty.gaussian(ccode),data=TC5) 

TC5.Frailty.Gaus.ExpDecay15.disag <-coxph(TC5.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_15 + Forceful_ExpDecay_15 + Peaceful_ExpDecay_15 + Promoted_ExpDecay_15
                                          +polity2 + area_1000 + lmtnest + ELF 
                                          +frailty.gaussian(ccode),data=TC5) 

TC5.Frailty.Gaus.ExpDecay20 <-coxph(TC5.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_20 + BadEnd_ExpDecay_20
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +frailty.gaussian(ccode),data=TC5) 

TC5.Frailty.Gaus.ExpDecay20.disag <-coxph(TC5.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_20 + Forceful_ExpDecay_20 + Peaceful_ExpDecay_20 + Promoted_ExpDecay_20
                                          +polity2 + area_1000 + lmtnest + ELF 
                                          +frailty.gaussian(ccode),data=TC5) 

##############################################

TC5or6.Frailty.Gaus.t10 <-coxph(TC5or6.Gap ~ Peace_wTC + Fighting_wTC   + GoodEnd_t10 + BadEnd_t10
                                +polity2 + area_1000 + lmtnest + ELF 
                                +frailty.gaussian(ccode),data=TC5or6) 

TC5or6.Frailty.Gaus.t20 <-coxph(TC5or6.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_t20 + BadEnd_t20
                                +polity2 + area_1000 + lmtnest + ELF 
                                +frailty.gaussian(ccode),data=TC5or6) 

TC5or6.Frailty.Gaus.tever <-coxph(TC5or6.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_end + Forceful_end + Peaceful_end + Promoted_end
                                  +polity2 + area_1000 + lmtnest + ELF 
                                  +frailty.gaussian(ccode),data=TC5or6) 

TC5or6.Frailty.Gaus.LinDecay <-coxph(TC5or6.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_LinDecay + BadEnd_LinDecay
                                       +polity2 + area_1000 + lmtnest + ELF 
                                     +frailty.gaussian(ccode),data=TC5or6) 

TC5or6.Frailty.Gaus.LinDecay.disag <-coxph(TC5or6.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_LinDecay + Forceful_LinDecay + Peaceful_LinDecay + Promoted_LinDecay
                                           +polity2 + area_1000 + lmtnest + ELF 
                                           +frailty.gaussian(ccode),data=TC5or6) 

TC5or6.Frailty.Gaus.ExpDecay5 <-coxph(TC5or6.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_5 + BadEnd_ExpDecay_5
                                      +polity2 + area_1000 + lmtnest + ELF 
                                      +frailty.gaussian(ccode),data=TC5or6) 

TC5or6.Frailty.Gaus.ExpDecay5.disag <-coxph(TC5or6.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_5 + Forceful_ExpDecay_5 + Peaceful_ExpDecay_5 + Promoted_ExpDecay_5
                                            +polity2 + area_1000 + lmtnest + ELF 
                                            +frailty.gaussian(ccode),data=TC5or6) 

TC5or6.Frailty.Gaus.ExpDecay10 <-coxph(TC5or6.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_10 + BadEnd_ExpDecay_10
                                       +polity2 + area_1000 + lmtnest + ELF 
                                       +frailty.gaussian(ccode),data=TC5or6)

TC5or6.Frailty.Gaus.ExpDecay10.disag <-coxph(TC5or6.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_10 + Forceful_ExpDecay_10 + Peaceful_ExpDecay_10 + Promoted_ExpDecay_10
                                             +polity2 + area_1000 + lmtnest + ELF 
                                             +frailty.gaussian(ccode),data=TC5or6)

TC5or6.Frailty.Gaus.ExpDecay15 <-coxph(TC5or6.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_15 + BadEnd_ExpDecay_15
                                       +polity2 + area_1000 + lmtnest + ELF 
                                       +frailty.gaussian(ccode),data=TC5or6)

TC5or6.Frailty.Gaus.ExpDecay15.disag <-coxph(TC5or6.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_15 + Forceful_ExpDecay_15 + Peaceful_ExpDecay_15 + Promoted_ExpDecay_15
                                             +polity2 + area_1000 + lmtnest + ELF 
                                             +frailty.gaussian(ccode),data=TC5or6)

TC5or6.Frailty.Gaus.ExpDecay20 <-coxph(TC5or6.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_20 + BadEnd_ExpDecay_20
                                       +polity2 + area_1000 + lmtnest + ELF 
                                       +frailty.gaussian(ccode),data=TC5or6)

TC5or6.Frailty.Gaus.ExpDecay20.disag <-coxph(TC5or6.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_20 + Forceful_ExpDecay_20 + Peaceful_ExpDecay_20 + Promoted_ExpDecay_20
                                             +polity2 + area_1000 + lmtnest + ELF 
                                             +frailty.gaussian(ccode),data=TC5or6)

##############################################

TC5orMore.Frailty.Gaus.t10 <-coxph(TC5orMore.Gap ~ Peace_wTC + Fighting_wTC + GoodEnd_t10 + BadEnd_t10
                                +polity2 + area_1000 + lmtnest + ELF 
                                +frailty.gaussian(ccode),data=TC5orMore) 

TC5orMore.Frailty.Gaus.t20 <-coxph(TC5orMore.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_t20 + BadEnd_t20
                                +polity2 + area_1000 + lmtnest + ELF 
                                +frailty.gaussian(ccode),data=TC5orMore) 

TC5orMore.Frailty.Gaus.tever <-coxph(TC5orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_end + Forceful_end + Peaceful_end + Promoted_end
                                  +polity2 + area_1000 + lmtnest + ELF 
                                  +frailty.gaussian(ccode),data=TC5orMore) 

TC5orMore.Frailty.Gaus.LinDecay <-coxph(TC5orMore.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_LinDecay + BadEnd_LinDecay
                                     +polity2 + area_1000 + lmtnest + ELF 
                                     +frailty.gaussian(ccode),data=TC5orMore) 

TC5orMore.Frailty.Gaus.LinDecay.disag <-coxph(TC5orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_LinDecay + Forceful_LinDecay + Peaceful_LinDecay + Promoted_LinDecay
                                        +polity2 + area_1000 + lmtnest + ELF 
                                        +frailty.gaussian(ccode),data=TC5orMore) 

TC5orMore.Frailty.Gaus.ExpDecay5 <-coxph(TC5orMore.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_5 + BadEnd_ExpDecay_5
                                      +polity2 + area_1000 + lmtnest + ELF 
                                      +frailty.gaussian(ccode),data=TC5orMore) 

TC5orMore.Frailty.Gaus.ExpDecay5.disag <-coxph(TC5orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_5 + Forceful_ExpDecay_5 + Peaceful_ExpDecay_5 + Promoted_ExpDecay_5
                                               +polity2 + area_1000 + lmtnest + ELF 
                                               +frailty.gaussian(ccode),data=TC5orMore) 

TC5orMore.Frailty.Gaus.ExpDecay10 <-coxph(TC5orMore.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_10 + BadEnd_ExpDecay_10
                                       +polity2 + area_1000 + lmtnest + ELF 
                                       +frailty.gaussian(ccode),data=TC5orMore) 

TC5orMore.Frailty.Gaus.ExpDecay10.disag <-coxph(TC5orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_10 + Forceful_ExpDecay_10 + Peaceful_ExpDecay_10 + Promoted_ExpDecay_10
                                                +polity2 + area_1000 + lmtnest + ELF 
                                                +frailty.gaussian(ccode),data=TC5orMore) 

TC5orMore.Frailty.Gaus.ExpDecay15 <-coxph(TC5orMore.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_15 + BadEnd_ExpDecay_15
                                       +polity2 + area_1000 + lmtnest + ELF 
                                       +frailty.gaussian(ccode),data=TC5orMore) 

TC5orMore.Frailty.Gaus.ExpDecay15.disag <-coxph(TC5orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_15 + Forceful_ExpDecay_15 + Peaceful_ExpDecay_15 + Promoted_ExpDecay_15
                                                +polity2 + area_1000 + lmtnest + ELF 
                                                +frailty.gaussian(ccode),data=TC5orMore) 

TC5orMore.Frailty.Gaus.ExpDecay20 <-coxph(TC5orMore.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_20 + BadEnd_ExpDecay_20
                                                  +polity2 + area_1000 + lmtnest + ELF 
                                       +frailty.gaussian(ccode),data=TC5orMore)

TC5orMore.Frailty.Gaus.ExpDecay20.disag <-coxph(TC5orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_20 + Forceful_ExpDecay_20 + Peaceful_ExpDecay_20 + Promoted_ExpDecay_20
                                                +polity2 + area_1000 + lmtnest + ELF 
                                                +frailty.gaussian(ccode),data=TC5orMore)

##############################################

TC6orMore.Frailty.Gaus.t10 <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC   + GoodEnd_t10 + BadEnd_t10
                                   +polity2 + area_1000 + lmtnest + ELF 
                                   +frailty.gaussian(ccode),data=TC6orMore) 

TC6orMore.Frailty.Gaus.t20 <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_t20 + BadEnd_t20
                                   +polity2 + area_1000 + lmtnest + ELF 
                                   +frailty.gaussian(ccode),data=TC6orMore) 

TC6orMore.Frailty.Gaus.tever <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_end + Forceful_end + Peaceful_end + Promoted_end
                                     +polity2 + area_1000 + lmtnest + ELF 
                                     +frailty.gaussian(ccode),data=TC6orMore) 

TC6orMore.Frailty.Gaus.LinDecay.disag <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_LinDecay + Forceful_LinDecay + Peaceful_LinDecay + Promoted_LinDecay
                                        +polity2 + area_1000 + lmtnest + ELF 
                                        +frailty.gaussian(ccode),data=TC6orMore) 

TC6orMore.Frailty.Gaus.LinDecay <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_LinDecay + BadEnd_LinDecay
                                              +polity2 + area_1000 + lmtnest + ELF 
                                              +frailty.gaussian(ccode),data=TC6orMore)

TC6orMore.Frailty.Gaus.ExpDecay5 <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_5 + BadEnd_ExpDecay_5
                                         +polity2 + area_1000 + lmtnest + ELF 
                                         +frailty.gaussian(ccode),data=TC6orMore) 

TC6orMore.Frailty.Gaus.ExpDecay5.disag <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_5 + Forceful_ExpDecay_5 + Peaceful_ExpDecay_5 + Promoted_ExpDecay_5
                                               +polity2 + area_1000 + lmtnest + ELF 
                                               +frailty.gaussian(ccode),data=TC6orMore) 

TC6orMore.Frailty.Gaus.ExpDecay10 <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_10 + BadEnd_ExpDecay_10
                                          +polity2 + area_1000 + lmtnest + ELF 
                                          +frailty.gaussian(ccode),data=TC6orMore) 

TC6orMore.Frailty.Gaus.ExpDecay10.disag <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_10 + Forceful_ExpDecay_10 + Peaceful_ExpDecay_10 + Promoted_ExpDecay_10
                                                +polity2 + area_1000 + lmtnest + ELF 
                                                +frailty.gaussian(ccode),data=TC6orMore) 

TC6orMore.Frailty.Gaus.ExpDecay15 <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_15 + BadEnd_ExpDecay_15
                                          +polity2 + area_1000 + lmtnest + ELF 
                                          +frailty.gaussian(ccode),data=TC6orMore) 

TC6orMore.Frailty.Gaus.ExpDecay15.disag <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_15 + Forceful_ExpDecay_15 + Peaceful_ExpDecay_15 + Promoted_ExpDecay_15
                                                +polity2 + area_1000 + lmtnest + ELF 
                                                +frailty.gaussian(ccode),data=TC6orMore) 

TC6orMore.Frailty.Gaus.ExpDecay20 <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_ExpDecay_20 + BadEnd_ExpDecay_20
                                          +polity2 + area_1000 + lmtnest + ELF 
                                          +frailty.gaussian(ccode),data=TC6orMore)

TC6orMore.Frailty.Gaus.ExpDecay20.disag <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_20 + Forceful_ExpDecay_20 + Peaceful_ExpDecay_20 + Promoted_ExpDecay_20
                                                +polity2 + area_1000 + lmtnest + ELF 
                                                +frailty.gaussian(ccode),data=TC6orMore)

##############################################

##############################################

TC7orMore.Frailty.Gaus.t10 <-coxph(TC7orMore.Gap ~ Peace_wTC + Fighting_wTC   + GoodEnd_t10 + BadEnd_t10
                    +polity2 + area_1000 + lmtnest + ELF 
                    +frailty.gaussian(ccode),data=TC7orMore) 

TC7orMore.Frailty.Gaus.t20 <-coxph(TC7orMore.Gap ~ Peace_wTC + Fighting_wTC + GoodEnd_t20 + BadEnd_t20
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gaussian(ccode),data=TC7orMore)

TC7orMore.Frailty.Gaus.tever <-coxph(TC7orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_end + Forceful_end + Peaceful_end + Promoted_end
                               +polity2 + area_1000 + lmtnest + ELF 
                               +frailty.gaussian(ccode),data=TC7orMore) 

TC7orMore.Frailty.Gaus.LinDecay.disag <-coxph(TC7orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_LinDecay + Forceful_LinDecay + Peaceful_LinDecay + Promoted_LinDecay
                                  +polity2 + area_1000 + lmtnest + ELF 
                                  +frailty.gaussian(ccode),data=TC7orMore) 

TC7orMore.Frailty.Gaus.ExpDecay5.disag <-coxph(TC7orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_5 + Forceful_ExpDecay_5 + Peaceful_ExpDecay_5 + Promoted_ExpDecay_5
                                   +polity2 + area_1000 + lmtnest + ELF 
                                   +frailty.gaussian(ccode),data=TC7orMore) 

TC7orMore.Frailty.Gaus.ExpDecay10.disag <-coxph(TC7orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_10 + Forceful_ExpDecay_10 + Peaceful_ExpDecay_10 + Promoted_ExpDecay_10
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +frailty.gaussian(ccode),data=TC7orMore) 

TC7orMore.Frailty.Gaus.ExpDecay15.disag <-coxph(TC7orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_15 + Forceful_ExpDecay_15 + Peaceful_ExpDecay_15 + Promoted_ExpDecay_15
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +frailty.gaussian(ccode),data=TC7orMore) 

TC7orMore.Frailty.Gaus.ExpDecay20.disag <-coxph(TC7orMore.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_ExpDecay_20 + Forceful_ExpDecay_20 + Peaceful_ExpDecay_20 + Promoted_ExpDecay_20
                                    +polity2 + area_1000 + lmtnest + ELF 
                                    +frailty.gaussian(ccode),data=TC7orMore) 

##############################################


summary(TC1.Frailty.Gaus)
summary(TC2.Frailty.Gaus.t10)
summary(TC3.Frailty.Gaus.t10)
summary(TC4.Frailty.Gaus.t10)
summary(TC5.Frailty.Gaus.t10)
summary(TC6.Frailty.Gaus.t10) # 10 events
summary(TC7.Frailty.Gaus.t10) # 6 events
summary(TC8.Frailty.Gaus.t10) # 4 events, goes haywire
# 9-15 would not converge. 

summary(TC2.Frailty.Gaus.t20)
summary(TC3.Frailty.Gaus.t20)
summary(TC4.Frailty.Gaus.t20)
summary(TC5.Frailty.Gaus.t20)
summary(TC6.Frailty.Gaus.t20) # NOT GOOD
summary(TC7.Frailty.Gaus.t20)
summary(TC8.Frailty.Gaus.t20)
summary(TC9.Frailty.Gaus.t20)
# 10-15 didn't converge

summary(TC2.Frailty.Gaus.tever)
summary(TC3.Frailty.Gaus.tever)
summary(TC4.Frailty.Gaus.tever)
summary(TC5.Frailty.Gaus.tever)
summary(TC6.Frailty.Gaus.tever) 
summary(TC7.Frailty.Gaus.tever)
summary(TC8.Frailty.Gaus.tever)
summary(TC9.Frailty.Gaus.tever) #No convergence

summary(TC2.Frailty.Gaus.LinDecay)
summary(TC3.Frailty.Gaus.LinDecay)
summary(TC4.Frailty.Gaus.LinDecay)
summary(TC5.Frailty.Gaus.LinDecay)
summary(TC6.Frailty.Gaus.LinDecay) 
summary(TC7.Frailty.Gaus.LinDecay)
summary(TC8.Frailty.Gaus.LinDecay)
summary(TC9.Frailty.Gaus.LinDecay)

summary(TC2.Frailty.Gaus.ExpDecay5)
summary(TC3.Frailty.Gaus.ExpDecay5)
summary(TC4.Frailty.Gaus.ExpDecay5)
summary(TC5.Frailty.Gaus.ExpDecay5)
summary(TC6.Frailty.Gaus.ExpDecay5) 
summary(TC7.Frailty.Gaus.ExpDecay5)
summary(TC8.Frailty.Gaus.ExpDecay5)
summary(TC9.Frailty.Gaus.ExpDecay5)

summary(TC2.Frailty.Gaus.ExpDecay10)
summary(TC3.Frailty.Gaus.ExpDecay10)
summary(TC4.Frailty.Gaus.ExpDecay10)
summary(TC5.Frailty.Gaus.ExpDecay10)
summary(TC6.Frailty.Gaus.ExpDecay10) 
summary(TC7.Frailty.Gaus.ExpDecay10)
summary(TC8.Frailty.Gaus.ExpDecay10)
summary(TC9.Frailty.Gaus.ExpDecay10)

summary(TC2.Frailty.Gaus.ExpDecay15)
summary(TC3.Frailty.Gaus.ExpDecay15)
summary(TC4.Frailty.Gaus.ExpDecay15)
summary(TC5.Frailty.Gaus.ExpDecay15)
summary(TC6.Frailty.Gaus.ExpDecay15) 
summary(TC7.Frailty.Gaus.ExpDecay15)
summary(TC8.Frailty.Gaus.ExpDecay15)
summary(TC9.Frailty.Gaus.ExpDecay15)

summary(TC2.Frailty.Gaus.ExpDecay20)
summary(TC3.Frailty.Gaus.ExpDecay20)
summary(TC4.Frailty.Gaus.ExpDecay20)
summary(TC5.Frailty.Gaus.ExpDecay20)
summary(TC6.Frailty.Gaus.ExpDecay20) 
summary(TC7.Frailty.Gaus.ExpDecay20)
summary(TC8.Frailty.Gaus.ExpDecay20)
summary(TC9.Frailty.Gaus.ExpDecay20)

summary(TC1.Frailty.Gaus)
summary(TC2or3.Frailty.Gaus.LinDecay)
summary(TC4or5.Frailty.Gaus.LinDecay)
summary(TC6orMore.Frailty.Gaus.LinDecay)

summary(TC1.Frailty.Gaus)
summary(TC2or3.Frailty.Gaus.ExpDecay15)
summary(TC4or5.Frailty.Gaus.ExpDecay15)
summary(TC6orMore.Frailty.Gaus.ExpDecay15)

# Stratified models for gamma frailty

TC1.Frailty.Gam <-coxph(TC1.Gap ~ polity2 + area_1000 + lmtnest + ELF 
                         +frailty.gamma(ccode),data=TC1) 

TC2.Frailty.Gam.t10 <-coxph(TC2.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gamma(ccode),data=TC2) 

TC2.Frailty.Gam.t20 <-coxph(TC2.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gamma(ccode),data=TC2) 

TC3.Frailty.Gam.t10 <-coxph(TC3.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gamma(ccode),data=TC3) 

TC3.Frailty.Gam.t20 <-coxph(TC3.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gamma(ccode),data=TC3) 

TC4.Frailty.Gam.t10 <-coxph(TC4.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gamma(ccode),data=TC4) 

TC4.Frailty.Gam.t20 <-coxph(TC4.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gamma(ccode),data=TC4) 

TC5.Frailty.Gam.t10 <-coxph(TC5.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gamma(ccode),data=TC5) 

TC5.Frailty.Gam.t20 <-coxph(TC5.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gamma(ccode),data=TC5) 

TC6.Frailty.Gam.t10 <-coxph(TC6.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gamma(ccode),data=TC6) 

TC6.Frailty.Gam.t20 <-coxph(TC6.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gamma(ccode),data=TC6) 

TC7.Frailty.Gam.t10 <-coxph(TC7.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gamma(ccode),data=TC7) 

TC7.Frailty.Gam.t20 <-coxph(TC7.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gamma(ccode),data=TC7) 

TC8.Frailty.Gam.t10 <-coxph(TC8.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gamma(ccode),data=TC8) 

TC8.Frailty.Gam.t20 <-coxph(TC8.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gamma(ccode),data=TC8) 

TC9.Frailty.Gam.t10 <-coxph(TC9.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gamma(ccode),data=TC9) 

TC9.Frailty.Gam.t20 <-coxph(TC9.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20
                             +polity2 + area_1000 + lmtnest + ELF 
                             +frailty.gamma(ccode),data=TC9) 

TC10.Frailty.Gam.t10 <-coxph(TC10.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                              +polity2 + area_1000 + lmtnest + ELF 
                              +frailty.gamma(ccode),data=TC10) 

TC10.Frailty.Gam.t20 <-coxph(TC10.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20
                              +polity2 + area_1000 + lmtnest + ELF 
                              +frailty.gamma(ccode),data=TC10) 

TC11.Frailty.Gam.t10 <-coxph(TC11.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                              +polity2 + area_1000 + lmtnest + ELF 
                              +frailty.gamma(ccode),data=TC11) 

TC11.Frailty.Gam.t20 <-coxph(TC11.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20
                              +polity2 + area_1000 + lmtnest + ELF 
                              +frailty.gamma(ccode),data=TC11)

TC12.Frailty.Gam.t10 <-coxph(TC12.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                              +polity2 + area_1000 + lmtnest + ELF 
                              +frailty.gamma(ccode),data=TC12) 

TC12.Frailty.Gam.t20 <-coxph(TC12.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20
                              +polity2 + area_1000 + lmtnest + ELF 
                              +frailty.gamma(ccode),data=TC12) 

TC13.Frailty.Gam.t10 <-coxph(TC13.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                              +polity2 + area_1000 + lmtnest + ELF 
                              +frailty.gamma(ccode),data=TC13) 

TC13.Frailty.Gam.t20 <-coxph(TC13.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20
                              +polity2 + area_1000 + lmtnest + ELF 
                              +frailty.gamma(ccode),data=TC13) 

TC14.Frailty.Gam.t10 <-coxph(TC14.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                              +polity2 + area_1000 + lmtnest + ELF 
                              +frailty.gamma(ccode),data=TC14) 

TC14.Frailty.Gam.t20 <-coxph(TC14.Gap ~ Peace_wTC + Fighting_wTC  + Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20
                              +polity2 + area_1000 + lmtnest + ELF 
                              +frailty.gamma(ccode),data=TC14) 

TC15.Frailty.Gam.t10 <-coxph(TC15.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                              +polity2 + area_1000 + lmtnest + ELF 
                              +frailty.gamma(ccode),data=TC15) 

TC15.Frailty.Gam.t20 <-coxph(TC15.Gap ~ Peace_wTC + Fighting_wTC   + Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20
                              +polity2 + area_1000 + lmtnest + ELF 
                              +frailty.gamma(ccode),data=TC15) 

summary(TC1.Frailty.Gam)
summary(TC2.Frailty.Gam.t10)
summary(TC3.Frailty.Gam.t10)
summary(TC4.Frailty.Gam.t10)
summary(TC5.Frailty.Gam.t10)
summary(TC6.Frailty.Gam.t10) 
summary(TC7.Frailty.Gam.t10) 
summary(TC8.Frailty.Gam.t10) 
 

summary(TC2.Frailty.Gam.t20)
summary(TC3.Frailty.Gam.t20)
summary(TC4.Frailty.Gam.t20)
summary(TC5.Frailty.Gam.t20)
summary(TC6.Frailty.Gam.t20) 
summary(TC7.Frailty.Gam.t20)
summary(TC8.Frailty.Gam.t20)
summary(TC9.Frailty.Gam.t20)

############################ Tables for 10-23 ######

stargazer(TCAll.PWP.LinDecay, TC.1.PWP, TC2or3.PWP.LinDecay, TC4or5.PWP.LinDecay, TC6orMore.PWP.LinDecay,
          title = "Conditional Gap Time Models with Linear Decay (20 Years)",
          covariate.labels = c("Peace w/ TC", "Fighting w/ TC", "Good End",
                               "Bad End", "Polity", "Area", "Mountains", "ELF"),
          dep.var.caption = "New TC Emergence",
          dep.var.labels.include = FALSE,
          column.labels = c("All Years","First TC", "TC 2-3", "TC 4-5", "TC 6+"), 
          df = FALSE)


stargazer(TCAll.PWP.LinDecay.disag, TC.1.PWP, TC2or3.PWP.LinDecay.disag, TC4or5.PWP.LinDecay.disag, TC6orMore.PWP.LinDecay.disag,
          title = "Conditional Gap Time Models with Linear Decay (20 Years) - Disaggregated",
          covariate.labels = c("Peace w/ TC", "Fighting w/ TC", "Absorbed" , "Forceful", "Peaceful",
                               "Promoted", "Polity", "Area", "Mountains", "ELF"),
          dep.var.caption = "New TC Emergence",
          dep.var.labels.include = FALSE,
          column.labels = c("All Years","First TC", "TC 2-3", "TC 4-5", "TC 6+"), 
          df = FALSE)


stargazer(TCAll.PWP.ExpDecay15, TC.1.PWP, TC2or3.PWP.ExpDecay15, TC4or5.PWP.ExpDecay15, TC6orMore.PWP.ExpDecay15,
          title = "Conditional Gap Time Models with Exponential Decay 15",
          covariate.labels = c("Peace w/ TC", "Fighting w/ TC", "Good End",
                               "Bad End", "Polity", "Area", "Mountains", "ELF"),
          dep.var.caption = "New TC Emergence",
          dep.var.labels.include = FALSE,
          column.labels = c("All Years","First TC", "TC 2-3", "TC 4-5", "TC 6+"), 
          df = FALSE
          )

stargazer(TCAll.PWP.ExpDecay15.disag, TC.1.PWP, TC2or3.PWP.ExpDecay15.disag, TC4or5.PWP.ExpDecay15.disag, TC6orMore.PWP.ExpDecay15.disag,
          title = "Conditional Gap Time Models with Exponential Decay 15 - Disaggregated",
          covariate.labels = c("Peace w/ TC", "Fighting w/ TC", "Absorbed" , "Forceful", "Peaceful",
                               "Promoted", "Polity", "Area", "Mountains", "ELF"),
          dep.var.caption = "New TC Emergence",
          dep.var.labels.include = FALSE,
          column.labels = c("All Years","First TC", "TC 2-3", "TC 4-5", "TC 6+"), 
          df = FALSE
)



##### Experimenting ###########################


IndoDRC <-  TCDat %>% filter(statename %in% c("Indonesia", "ZaireDRC"))
Forceful <- filter(IndoDRC, Forceful == 1)
Peaceful <- filter(IndoDRC, Peaceful == 1 | Promoted == 1)

Indo <-  TCDat %>% filter(statename %in% c("Indonesia"))
ForcefulI <- filter(Indo, Forceful == 1)
PeacefulI <- filter(Indo, Peaceful == 1 | Promoted == 1)

DRC <-  TCDat %>% filter(statename %in% c("ZaireDRC"))
ForcefulD <- filter(DRC, Forceful == 1)
PeacefulD <- filter(DRC, Peaceful == 1 | Promoted == 1)

IndoDRCPlot1 <- ggplot() +
  geom_line(aes(y = TC_tally, x = year, colour = statename), size = 1, data = IndoDRC, stat="identity") +
  geom_point(aes(y = TC_tally, x = year, colour = statename), shape = 21, fill = "red", size = 2, data = Forceful, stat="identity") +
  geom_point(aes(y = TC_tally, x = year, colour = statename), shape = 21, size = 2,fill = "blue", data = Peaceful, stat="identity") +
  scale_x_continuous(breaks=seq(1945,2010,5))+
  ggtitle("Comparison of TCs Experienced by Indonesia and DRC") + labs(x="Year", y="Total Number of TCs")+
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank())+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        legend.key=element_rect(fill="white", colour="white"))

IndoDRCPlot1

IndoDRCPlot2 <- ggplot() +
  geom_line(aes(y = num_TC, x = year, colour = statename), size = 1, data = IndoDRC, stat="identity") +
  geom_point(aes(y = num_TC, x = year, colour = statename), shape = 21, fill = "red", size = 2, data = Forceful, stat="identity") +
  geom_point(aes(y = num_TC, x = year, colour = statename), shape = 21, size = 2,fill = "blue", data = Peaceful, stat="identity") +
  scale_x_continuous(breaks=seq(1945,2010,5))+
  scale_y_continuous(breaks=seq(0,7,1))+
  ggtitle("Comparison of TCs Experienced by Indonesia and DRC") + labs(x="Year", y="Number of Active TCs")+
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank())+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        legend.key=element_rect(fill="white", colour="white"))

IndoDRCPlot2 


IndoPlot <- ggplot() +
  geom_line(aes(y = num_TC, x = year, colour = statename), size = 1, data = Indo, stat="identity") +
  geom_point(aes(y = num_TC, x = year, colour = statename), shape = 21, fill = "red", size = 2, data = ForcefulI, stat="identity") +
  geom_point(aes(y = num_TC, x = year, colour = statename), shape = 21, size = 2,fill = "blue", data = PeacefulI, stat="identity") +
  scale_x_continuous(breaks=seq(1945,2010,5))+
  scale_y_continuous(breaks=seq(0,7,1))+
  ggtitle("Comparison of TCs Experienced by Indonesia") + labs(x="Year", y="Number of Active TCs")+
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank())+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        legend.key=element_rect(fill="white", colour="white"))

IndoPlot

DRCPlot <- ggplot() +
  geom_line(aes(y = num_TC, x = year), size = 1, data = DRC, stat="identity") +
  geom_point(aes(y = num_TC, x = year), shape = 21, fill = "red", size = 2, data = ForcefulD, stat="identity") +
  geom_point(aes(y = num_TC, x = year), shape = 21, size = 2,fill = "blue", data = PeacefulD, stat="identity") +
  scale_y_continuous(limits = c(0,7), breaks=seq(0,7,1))+
  scale_x_continuous(limits = c(1945,2010), breaks=seq(1945,2010,5))+
  ggtitle("Comparison of TCs Experienced by DRC") + labs(x="Year", y="Number of Active TCs")+
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank())+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        legend.key=element_rect(fill="white", colour="white"))

DRCPlot

grid.newpage()
grid.draw(rbind(ggplotGrob(IndoPlot), ggplotGrob(DRCPlot), size = "last"))



IndoDRCPlot3 <- ggplot() +
  geom_line(aes(y = num_TC, x = year), size = 1, data = IndoDRC, stat="identity") +
  geom_point(aes(y = num_TC, x = year), shape = 21, fill = "red", size = 2, data = Forceful, stat="identity") +
  geom_point(aes(y = num_TC, x = year), shape = 21, size = 2,fill = "blue", data = Peaceful, stat="identity") +
  scale_x_continuous(breaks=seq(1945,2010,5))+
  scale_y_continuous(limits = c(0,7), breaks=seq(0,7,1))+
  facet_grid(statename ~ ., scales = "free_y")+
  ggtitle("Comparison of TCs Experienced by Indonesia and DRC") + labs(x="Year", y="Number of Active TCs")+
  theme(legend.position="bottom", legend.direction="horizontal", legend.title = element_blank())+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        legend.key=element_rect(fill="white", colour="white"))


IndoDRCPlot3 
