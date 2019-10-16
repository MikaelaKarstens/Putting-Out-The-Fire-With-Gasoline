# Repeat Failure Model

library(foreign) # download .dta files for Stata < 13
library(haven) # download .dta files for Stata 13+
library(sampleSelection) # Package containing heckit for heckprob substitute
library(dplyr) #Used for dataset laterations and survival time calculations
library(lubridate) #Used for survival data formatting
library(survival)
library(stargazer)
library(reReg) #recurrent events package

options(scipen = 999) # no scientific notation

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


