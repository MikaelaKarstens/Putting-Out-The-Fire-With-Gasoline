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

options(scipen = 999) # no scientific notation

# setwd("C:/Users/mjw65/Dropbox/State Making Strategies/State-Making-Strategies") #209 Pond Lab

TCDat <- read.csv("Survival-Data-TC.csv")
names(TCDat)

# Summary Statistics & Crosstabs

CrossTable(TCDat$new_TC_dummy, TCDat$Absorbed_t20)

CrossTable(TCDat$new_TC_dummy, TCDat$Forceful_t20)

CrossTable(TCDat$new_TC_dummy, TCDat$Peaceful_t20)

CrossTable(TCDat$new_TC_dummy, TCDat$Promoted_t20)

CrossTable(TCDat$new_TC_dummy, TCDat$TC_persists_alt)

CrossTable(TCDat$new_TC_dummy, TCDat$Peace_wTC)

CrossTable(TCDat$new_TC_dummy, TCDat$Fighting_wTC)

CrossTable(TCDat$new_TC_dummy, TCDat$polity2)

plot(TCDat$new_TC_dummy, TCDat$area_1000)

plot(TCDat$new_TC_dummy, TCDat$lmtnest)

plot(TCDat$new_TC_dummy, TCDat$ELF)

plot(TCDat$new_TC_dummy, TCDat$Country_Age)

# Creating subset datasets

TC2On <- TCDat[TCDat$event_num>=2,]

TC1 <- TCDat[TCDat$event_num==1,]

TC2 <- TCDat[TCDat$event_num==2,]

TC3 <- TCDat[TCDat$event_num==3,]

TC4 <- TCDat[TCDat$event_num==4,]

TC5 <- TCDat[TCDat$event_num==5,]

TC6 <- TCDat[TCDat$event_num==6,]

TC7 <- TCDat[TCDat$event_num==7,]

TC8 <- TCDat[TCDat$event_num==8,]

TC9 <- TCDat[TCDat$event_num==9,]

TC10 <- TCDat[TCDat$event_num==10,]

TC11 <- TCDat[TCDat$event_num==11,]

TC12 <- TCDat[TCDat$event_num==12,]

TC13 <- TCDat[TCDat$event_num==13,]

TC14 <- TCDat[TCDat$event_num==14,]

TC15 <- TCDat[TCDat$event_num==15,]


# Creating all survival objects

TC2On.Elapsed<-Surv(TC2On$alt_start,TC2On$alt_stop,TC2On$new_TC_dummy)
TC2On.Gap<-Surv(TC2On$start,TC2On$stop,TC2On$new_TC_dummy)

TC1.Elapsed<-Surv(TC1$alt_start,TC1$alt_stop,TC1$new_TC_dummy)
TC1.Gap<-Surv(TC1$start,TC1$stop,TC1$new_TC_dummy)

TC2.Elapsed<-Surv(TC2$alt_start,TC2$alt_stop,TC2$new_TC_dummy)
TC2.Gap<-Surv(TC2$start,TC2$stop,TC2$new_TC_dummy)

TC3.Elapsed<-Surv(TC3$alt_start,TC3$alt_stop,TC3$new_TC_dummy)
TC3.Gap<-Surv(TC3$start,TC3$stop,TC3$new_TC_dummy)

TC4.Elapsed<-Surv(TC4$alt_start,TC4$alt_stop,TC4$new_TC_dummy)
TC4.Gap<-Surv(TC4$start,TC4$stop,TC4$new_TC_dummy)

TC5.Elapsed<-Surv(TC5$alt_start,TC5$alt_stop,TC5$new_TC_dummy)
TC5.Gap<-Surv(TC5$start,TC5$stop,TC5$new_TC_dummy)

TC6.Elapsed<-Surv(TC6$alt_start,TC6$alt_stop,TC6$new_TC_dummy)
TC6.Gap<-Surv(TC6$start,TC6$stop,TC6$new_TC_dummy)

TC7.Elapsed<-Surv(TC7$alt_start,TC7$alt_stop,TC7$new_TC_dummy)
TC7.Gap<-Surv(TC7$start,TC7$stop,TC7$new_TC_dummy)

TC8.Elapsed<-Surv(TC8$alt_start,TC8$alt_stop,TC8$new_TC_dummy)
TC8.Gap<-Surv(TC8$start,TC8$stop,TC8$new_TC_dummy)

TC9.Elapsed<-Surv(TC9$alt_start,TC9$alt_stop,TC9$new_TC_dummy)
TC9.Gap<-Surv(TC9$start,TC9$stop,TC9$new_TC_dummy)

TC10.Elapsed<-Surv(TC10$alt_start,TC10$alt_stop,TC10$new_TC_dummy)
TC10.Gap<-Surv(TC10$start,TC10$stop,TC10$new_TC_dummy)

TC11.Elapsed<-Surv(TC11$alt_start,TC11$alt_stop,TC11$new_TC_dummy)
TC11.Gap<-Surv(TC11$start,TC11$stop,TC11$new_TC_dummy)

TC12.Elapsed<-Surv(TC12$alt_start,TC12$alt_stop,TC12$new_TC_dummy)
TC12.Gap<-Surv(TC12$start,TC12$stop,TC12$new_TC_dummy)

TC13.Elapsed<-Surv(TC13$alt_start,TC13$alt_stop,TC13$new_TC_dummy)
TC13.Gap<-Surv(TC13$start,TC13$stop,TC13$new_TC_dummy)

TC14.Elapsed<-Surv(TC14$alt_start,TC14$alt_stop,TC14$new_TC_dummy)
TC14.Gap<-Surv(TC14$start,TC14$stop,TC14$new_TC_dummy)

TC15.Elapsed<-Surv(TC15$alt_start,TC15$alt_stop,TC15$new_TC_dummy)
TC15.Gap<-Surv(TC15$start,TC15$stop,TC15$new_TC_dummy)

# Anderson-Gill Model - Variance Correction

TC.AG.2On.1 <- coxph(TC2On.Gap ~ Peace_wTC + Fighting_wTC + Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20 
                 +cluster(ccode),data=TC2On, method="efron")
summary(TC.AG.2On.1)

TC.AG.2On.2 <- coxph(TC2On.Gap ~ Peace_wTC + Fighting_wTC + Absorbed_t20 + Forceful_t20 + Peaceful_t20 + Promoted_t20 
                     +cluster(ccode),data=TC2On, method="efron")
summary(TC.AG.2On.2)




