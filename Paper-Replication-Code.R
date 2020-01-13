# Replication code for "Assessing State Making Strategies" by Douglas Lemke and Mikaela Karstens

# Required packages for replication - Please install if needed (FIX THESE TO REMOVE UNUSED)

library(survival)
library(stargazer)
library(dplyr)
library(ggplot2)

library(sampleSelection) # Package containing heckit for heckprob substitute
library(lubridate) #Used for survival data formatting
library(reReg) #recurrent events package
library(gmodels)
library(Hmisc)
library(survminer)
library(visreg)
library(coefplot)
library(broom)
library(GGally)

# Settings

options(scipen = 50) # no scientific notation

# Loading data

TCDat <- read.csv("Survival-Data-TC.csv")
names(TCDat)

summary(TCDat$area_1000)
TCDat$area_1000_log <- log(TCDat$area_1000)
summary(TCDat$area_1000_log)

# Creating subset datasets

TC2On <- TCDat[TCDat$event_num>=2,]

TC2or3 <- filter(TCDat, event_num == 2 | event_num == 3)

TC4or5 <- filter(TCDat, event_num == 4 | event_num == 5)

TC6orMore <- filter(TCDat, event_num >= 6)

TC1 <- filter(TCDat, event_num == 1)

TC2 <- filter(TCDat, event_num == 2)

TC3 <- filter(TCDat, event_num == 3)

TC4 <- filter(TCDat, event_num == 4)

TC5 <- filter(TCDat, event_num == 5)

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

TC4.Elapsed<-Surv(TC4$alt_start,TC4$alt_stop,TC4$new_TC_dummy)
TC4.Gap<-Surv(TC4$start,TC4$stop,TC4$new_TC_dummy)

TC4or5.Elapsed<-Surv(TC4or5$alt_start,TC4or5$alt_stop,TC4or5$new_TC_dummy)
TC4or5.Gap<-Surv(TC4or5$start,TC4or5$stop,TC4or5$new_TC_dummy)

TC5.Elapsed<-Surv(TC5$alt_start,TC5$alt_stop,TC5$new_TC_dummy)
TC5.Gap<-Surv(TC5$start,TC5$stop,TC5$new_TC_dummy)

TC6.Elapsed<-Surv(TC6$alt_start,TC6$alt_stop,TC6$new_TC_dummy)
TC6.Gap<-Surv(TC6$start,TC6$stop,TC6$new_TC_dummy)

TC6orMore.Elapsed<-Surv(TC6orMore$alt_start,TC6orMore$alt_stop,TC6orMore$new_TC_dummy)
TC6orMore.Gap<-Surv(TC6orMore$start,TC6orMore$stop,TC6orMore$new_TC_dummy)

####### Linear Decay PWP Models #####

TCAll.PWP <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_LinDecay + GoodEnd_LinDecay 
                           +polity2 + area_1000_log + lmtnest + ELF 
                           +strata(event_num) + cluster(ccode),data=TCDat, method="efron")

TC.1.PWP <- coxph(TC1.Gap ~ +polity2 + area_1000_log + lmtnest + ELF 
                  +strata(event_num) + cluster(ccode),data=TC1, method="efron")

TC2or3.PWP.LinDecay <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_LinDecay + GoodEnd_LinDecay 
                            +polity2 + area_1000_log + lmtnest + ELF 
                            +strata(event_num) + cluster(ccode),data=TC2or3, method="efron") 

TC4or5.PWP.LinDecay <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_LinDecay + GoodEnd_LinDecay 
                            +polity2 + area_1000_log + lmtnest + ELF 
                            +strata(event_num) + cluster(ccode),data=TC4or5, method="efron") 

TC6orMore.PWP.LinDecay <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_LinDecay + GoodEnd_LinDecay 
                               +polity2 + area_1000_log + lmtnest + ELF 
                               +strata(event_num) + cluster(ccode),data=TC6orMore, method="efron") 

summary(TCAll.PWP)

# Main Results Table

stargazer(TCAll.PWP, TC.1.PWP, TC2or3.PWP.LinDecay, TC4or5.PWP.LinDecay, TC6orMore.PWP.LinDecay,
          type = "latex", #out = "C:/Users/mjw65/Dropbox/State Making Strategies/State-Making-Strategies/table.html", 
          title = "PWP Gap Time Model Results",
          dep.var.labels = c("All TCs", "First TC", " TC 2 or 3", "TC 4 or 5", "TC 6+"),
          covariate.labels=c("Peace w/ TC", "Fighting w/ TC", "Forceful Reintegration", "Favorable Outcome for TC", "Polity 2 Score", "State Area (logged)","Mountains", "ELF"),
          keep.stat = c("n"))

# Coefficient Plot