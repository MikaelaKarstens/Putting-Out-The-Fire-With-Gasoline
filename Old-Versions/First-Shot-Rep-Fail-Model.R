# Repeat Failure Model 

library(foreign) # download .dta files for Stata < 13
library(haven) # download .dta files for Stata 13+
library(sampleSelection) # Package containing heckit for heckprob substitute
library(dplyr) #Used for dataset laterations and survival time calculations
library(lubridate) #Used for survival data formatting
library(survival)
library(stargazer)
library(reReg) #recurrent events package

options(scipen = 50) # no scientific notation

#setwd("~/State-Making-Strategies")
data <- read_dta("Deterring_TC.dta") #Uses haven

names(data) #print variable names to make sure everything read in alright.

SDat <- read.csv("Repeat-Failure-Data.csv")
names(SDat)
summary(SDat)
table(SDat$t.start)
table(SDat$t.stop)
table(SDat$event)
table(SDat$statename)

##### Work from 10-15 and 10/16 #######

# Time to first event

TCDat <- read.csv("Survival-Data-TC.csv")
names(TCDat)

reObj <- Surv(TCDat$alt_start,TCDat$alt_stop,TCDat$new_TC_dummy)
plotEvents(reObj, TCDat)

TC1st <- TCDat[TCDat$event_num==1,]
TC.1st<-Surv(TC1st$alt_start,TC1st$alt_stop,TC1st$new_TC_dummy)

TC.Cox.1st <-coxph(TC.1st~ polity2 + area_1000 + lmtnest + ELF + Country_Age + 
                     cluster(ccode),data=TC1st,method="efron")
summary(TC.Cox.1st)

stargazer(TC.Cox.1st, single.row = T)

# AG model with time since obs started

TC.AGS<-Surv(TCDat$alt_start,TCDat$alt_stop,TCDat$new_TC_dummy) 

TC.Cox.AG<-coxph(TC.AGS~Peace_wTC + Fighting_wTC + Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20 +cluster(ccode),data=TCDat,method="efron")
summary(TC.Cox.AG)
#Why are Peace_wTC and Fighting_wTC curent year while others are past? Try without.


TC.Cox.AG2<-coxph(TC.AGS~Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20 +cluster(ccode),data=TCDat,method="efron")
summary(TC.Cox.AG2)

#Add controls

TC.Cox.AG3<-coxph(TC.AGS~Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20 +
                    polity2 + area_1000 + lmtnest + ELF + Country_Age +  
                    cluster(ccode),data=TCDat,method="efron")
summary(TC.Cox.AG3)

# Redo models using gap time

TC.AGS2<-Surv(TCDat$start,TCDat$stop,TCDat$new_TC_dummy) 

TC.Cox.AG4<-coxph(TC.AGS2~Peace_wTC + Fighting_wTC + Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20 +cluster(ccode),data=TCDat,method="efron")
summary(TC.Cox.AG4)

TC.Cox.AG5<-coxph(TC.AGS2~Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20 +cluster(ccode),data=TCDat,method="efron")
summary(TC.Cox.AG5)

#Add controls

TC.Cox.AG6<-coxph(TC.AGS~Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20 +
                    polity2 + area_1000 + lmtnest + ELF + Country_Age + 
                    cluster(ccode),data=TCDat,method="efron")
summary(TC.Cox.AG6)
stargazer(TC.Cox.AG6, single.row = T)

##### Now moving on the the conditional elapsed time models

TC.PWPES<-Surv(TCDat$alt_start,TCDat$alt_stop,TCDat$new_TC_dummy) 

TC.Cox.PWPE<-coxph(TC.PWPES~ Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20
                     +polity2 + area_1000 + lmtnest + ELF + Country_Age 
                   +strata(event_num)+cluster(ccode),data=TCDat, method="efron")

summary(TC.Cox.PWPE)

# Now with gap time

TC.PWPES2<-Surv(TCDat$start,TCDat$stop,TCDat$new_TC_dummy) 

TC.Cox.PWPE2<-coxph(TC.PWPES2 ~ Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20
                   +polity2 + area_1000 + lmtnest + ELF + Country_Age
                   +strata(event_num)+cluster(ccode),data=TCDat, method="efron")

summary(TC.Cox.PWPE2)

stargazer(TC.Cox.PWPE2, single.row = T)

########## Cox gap time with gaussian frailty
s <- Surv(TCDat$start,TCDat$stop,TCDat$new_TC_dummy) 
cox.F<-coxph(s~ Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20
             +polity2 + area_1000 + lmtnest + ELF + Country_Age
             +frailty.gaussian(ccode),data=TCDat) 
summary(cox.F) 
stargazer(cox.F)

########## Cox gap time with gamma frailty
s <- Surv(TCDat$start,TCDat$stop,TCDat$new_TC_dummy) 
cox.F2<-coxph(s~ Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20
             +polity2 + area_1000 + lmtnest + ELF + Country_Age
             +frailty.gamma(ccode),data=TCDat) 
summary(cox.F2) 
stargazer(cox.F2)

###### Only set of countries with more than 1 TC

TC2On <- TCDat[TCDat$event_num>=2,]
TC.2On.Elapsed<-Surv(TC2On$alt_start,TC2On$alt_stop,TC2On$new_TC_dummy)
TC.2On.Gap<-Surv(TC2On$start,TC2On$stop,TC2On$new_TC_dummy)

TC.Cox.PWPE.2On<-coxph(TC.2On.Gap ~ Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20
                    +polity2 + area_1000 + lmtnest + ELF + Country_Age
                    +strata(event_num)+cluster(ccode),data=TC2On, method="efron")

summary(TC.Cox.PWPE.2On)


cox.F.2On<-coxph(TC.2On.Gap ~ Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20
              +polity2 + area_1000 + lmtnest + ELF + Country_Age
              +frailty.gaussian(ccode),data=TC2On) 

summary(cox.F.2On)

cox.F.2On.t10<-coxph(TC.2On.Gap ~ Fighting_wTC + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                 +polity2 + area_1000 + lmtnest + ELF
                 +frailty.gaussian(ccode),data=TC2On) 

summary(cox.F.2On.t10)

cox.F.2On.t10.persist<-coxph(TC.2On.Gap ~ TC_persists_alt + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                     +polity2 + area_1000 + lmtnest + ELF
                     +frailty.gaussian(ccode),data=TC2On) 

summary(cox.F.2On.t10.persist)
# DOES NOT WORK!!!!!!!!!!!!!!!



cox.F.2On2<-coxph(TC.2On.Gap ~ Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20
                 +polity2 + area_1000 + lmtnest + ELF + Country_Age
                 +frailty.gamma(ccode),data=TC2On) 

summary(cox.F.2On2)

###### Full Strata ########

TC1 <- TCDat[TCDat$event_num==1,]
TC1.Elapsed<-Surv(TC1$alt_start,TC1$alt_stop,TC1$new_TC_dummy)
TC1.Gap<-Surv(TC1$start,TC1$stop,TC1$new_TC_dummy)

TC2 <- TCDat[TCDat$event_num==2,]
TC2.Elapsed<-Surv(TC2$alt_start,TC2$alt_stop,TC2$new_TC_dummy)
TC2.Gap<-Surv(TC2$start,TC2$stop,TC2$new_TC_dummy)

TC3 <- TCDat[TCDat$event_num==3,]
TC3.Elapsed<-Surv(TC3$alt_start,TC3$alt_stop,TC3$new_TC_dummy)
TC3.Gap<-Surv(TC3$start,TC3$stop,TC3$new_TC_dummy)

TC4 <- TCDat[TCDat$event_num==4,]
TC4.Elapsed<-Surv(TC4$alt_start,TC4$alt_stop,TC4$new_TC_dummy)
TC4.Gap<-Surv(TC4$start,TC4$stop,TC4$new_TC_dummy)

TC5 <- TCDat[TCDat$event_num==5,]
TC5.Elapsed<-Surv(TC5$alt_start,TC5$alt_stop,TC5$new_TC_dummy)
TC5.Gap<-Surv(TC5$start,TC5$stop,TC5$new_TC_dummy)

TC6 <- TCDat[TCDat$event_num==6,]
TC6.Elapsed<-Surv(TC6$alt_start,TC6$alt_stop,TC6$new_TC_dummy)
TC6.Gap<-Surv(TC6$start,TC6$stop,TC6$new_TC_dummy)

TC7 <- TCDat[TCDat$event_num==7,]
TC7.Elapsed<-Surv(TC7$alt_start,TC7$alt_stop,TC7$new_TC_dummy)
TC7.Gap<-Surv(TC7$start,TC7$stop,TC7$new_TC_dummy)

TC8 <- TCDat[TCDat$event_num==8,]
TC8.Elapsed<-Surv(TC8$alt_start,TC8$alt_stop,TC8$new_TC_dummy)
TC8.Gap<-Surv(TC8$start,TC8$stop,TC8$new_TC_dummy)

TC9 <- TCDat[TCDat$event_num==9,]
TC9.Elapsed<-Surv(TC9$alt_start,TC9$alt_stop,TC9$new_TC_dummy)
TC9.Gap<-Surv(TC9$start,TC9$stop,TC9$new_TC_dummy)

TC10 <- TCDat[TCDat$event_num==10,]
TC10.Elapsed<-Surv(TC10$alt_start,TC10$alt_stop,TC10$new_TC_dummy)
TC10.Gap<-Surv(TC10$start,TC10$stop,TC10$new_TC_dummy)

TC11 <- TCDat[TCDat$event_num==11,]
TC11.Elapsed<-Surv(TC11$alt_start,TC11$alt_stop,TC11$new_TC_dummy)
TC11.Gap<-Surv(TC11$start,TC11$stop,TC11$new_TC_dummy)

TC12 <- TCDat[TCDat$event_num==12,]
TC12.Elapsed<-Surv(TC12$alt_start,TC12$alt_stop,TC12$new_TC_dummy)
TC12.Gap<-Surv(TC12$start,TC12$stop,TC12$new_TC_dummy)

TC13 <- TCDat[TCDat$event_num==13,]
TC13.Elapsed<-Surv(TC13$alt_start,TC13$alt_stop,TC13$new_TC_dummy)
TC13.Gap<-Surv(TC13$start,TC13$stop,TC13$new_TC_dummy)

TC14 <- TCDat[TCDat$event_num==14,]
TC14.Elapsed<-Surv(TC14$alt_start,TC14$alt_stop,TC14$new_TC_dummy)
TC14.Gap<-Surv(TC14$start,TC14$stop,TC14$new_TC_dummy)

TC15 <- TCDat[TCDat$event_num==15,]
TC15.Elapsed<-Surv(TC15$alt_start,TC15$alt_stop,TC15$new_TC_dummy)
TC15.Gap<-Surv(TC15$start,TC15$stop,TC15$new_TC_dummy)

TC1.Frailty <-coxph(TC1.Gap ~ TC_persists_alt + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                             +polity2 + area_1000 + lmtnest + ELF + Country_Age
                             +frailty.gaussian(ccode),data=TC1) 

TC2.Frailty <-coxph(TC2.Gap ~ TC_persists_alt + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                    +polity2 + area_1000 + lmtnest + ELF + Country_Age
                    +frailty.gaussian(ccode),data=TC2) 

TC3.Frailty <-coxph(TC3.Gap ~ TC_persists_alt + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                    +polity2 + area_1000 + lmtnest + ELF + Country_Age
                    +frailty.gaussian(ccode),data=TC3) 

TC4.Frailty <-coxph(TC4.Gap ~ TC_persists_alt + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                    +polity2 + area_1000 + lmtnest + ELF + Country_Age
                    +frailty.gaussian(ccode),data=TC4) 

TC5.Frailty <-coxph(TC5.Gap ~ TC_persists_alt + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                    +polity2 + area_1000 + lmtnest + ELF + Country_Age
                    +frailty.gaussian(ccode),data=TC5) 

TC6.Frailty <-coxph(TC6.Gap ~ TC_persists_alt + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                    +polity2 + area_1000 + lmtnest + ELF + Country_Age
                    +frailty.gaussian(ccode),data=TC6) 

TC7.Frailty <-coxph(TC7.Gap ~ TC_persists_alt + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                    +polity2 + area_1000 + lmtnest + ELF + Country_Age
                    +frailty.gaussian(ccode),data=TC7) 

TC8.Frailty <-coxph(TC8.Gap ~ TC_persists_alt + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                    +polity2 + area_1000 + lmtnest + ELF + Country_Age
                    +frailty.gaussian(ccode),data=TC8) 

#These models don't work past here. At TC 12, no longer enough df. 
TC9.Frailty <-coxph(TC9.Gap ~ TC_persists_alt + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                    +polity2 + area_1000 + lmtnest + ELF + Country_Age
                    +frailty.gaussian(ccode),data=TC9) 

TC10.Frailty <-coxph(TC10.Gap ~ TC_persists_alt + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                    +polity2 + area_1000 + lmtnest + ELF + Country_Age
                    +frailty.gaussian(ccode),data=TC10) 

TC11.Frailty <-coxph(TC11.Gap ~ TC_persists_alt + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                    +polity2 + area_1000 + lmtnest + ELF + Country_Age
                    +frailty.gaussian(ccode),data=TC11) 

TC12.Frailty <-coxph(TC12.Gap ~ TC_persists_alt + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                    +polity2 + area_1000 + lmtnest + ELF + Country_Age
                    +frailty.gaussian(ccode),data=TC12) 

TC13.Frailty <-coxph(TC13.Gap ~ TC_persists_alt + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                    +polity2 + area_1000 + lmtnest + ELF + Country_Age
                    +frailty.gaussian(ccode),data=TC13) 

TC14.Frailty <-coxph(TC14.Gap ~ TC_persists_alt + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                    +polity2 + area_1000 + lmtnest + ELF + Country_Age
                    +frailty.gaussian(ccode),data=TC14) 

TC15.Frailty <-coxph(TC15.Gap ~ TC_persists_alt + Absorbed_t10 + Forceful_t10 + Peaceful_t10 + Promoted_t10
                    +polity2 + area_1000 + lmtnest + ELF + Country_Age
                    +frailty.gaussian(ccode),data=TC15) 

summary(TC1.Frailty)
summary(TC2.Frailty)
summary(TC3.Frailty)
summary(TC4.Frailty)
summary(TC5.Frailty)
summary(TC6.Frailty)
summary(TC7.Frailty)
summary(TC8.Frailty)
summary(TC9.Frailty)
summary(TC10.Frailty)
summary(TC11.Frailty)
summary(TC12.Frailty)
summary(TC13.Frailty)
summary(TC14.Frailty)
summary(TC15.Frailty)


###########

reObj <- with(SDat, reSurv(t.start, t.stop, statename, event))
plotEvents(reObj~SDat, SDat)

write.csv(data, "DougData.csv")

summary(reObj)

function (reObj ~ ccode, data=SDat, B = 200, method = c("cox.LWYY", "cox.GL", 
                                             "cox.HW", "am.GL", "am.XCHWY",
                                             "sc.XCYH"), se = c("NULL",
                                              "bootstrap", "resampling"), 
                                          contrasts = NULL, control = list()) 


############################################################################################################# Testing ##

#Replicating vignette for reReg
install.packages("frailtypack")
library(frailtypack)
library(reReg)

data(readmission, package = "frailtypack") #getting data


# Failed Experimenting with package from Therneau 2019 at Mayo Clinic

install.packages("coxme")
install.packages("kinship2")
library(kinship2)
library(coxme)
data(minnbreast)

makefig <- function(file) {
  pdf(paste(file, "pdf", sep='.'), width=7, height=5)
  par(mar=c(5.1, 4.1, .1, .1))
                          }

names(minnbreast)

with(minnbreast, table(sex, cancer, exclude=NULL))

mped <- with(minnbreast, pedigree(id, fatherid, motherid, sex,
              affected=cancer, famid=famid, status=proband))

plot(mped["8"]) #figure 1. This is the pedigree for a family in this case. Black is breast cancer. 


