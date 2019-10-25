##### Replication File #####


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
library(coefplot)
library(broom)
library(ggplot2)
library(GGally)



options(scipen = 50) # no scientific notation

# setwd("C:/Users/mjw65/Dropbox/State Making Strategies/State-Making-Strategies") #209 Pond Lab

TCDat <- read.csv("Survival-Data-TC.csv")
names(TCDat)


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

####### Linear Decay PWP Models #####

TCDat$PrevTC <- TCDat$event_num - 1

TCAll.PWP.LinDecay <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_LinDecay + BadEnd_LinDecay
                           +polity2 + area_1000 + lmtnest + ELF 
                           +strata(event_num) + cluster(ccode),data=TCDat, method="efron")

TC.1.PWP <- coxph(TC1.Gap ~ +polity2 + area_1000 + lmtnest + ELF 
                  +strata(event_num) + cluster(ccode),data=TC1, method="efron")

TC2or3.PWP.LinDecay <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_LinDecay + BadEnd_LinDecay
                            +polity2 + area_1000 + lmtnest + ELF 
                            +strata(event_num) + cluster(ccode),data=TC2or3, method="efron") 

TC4or5.PWP.LinDecay <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_LinDecay + BadEnd_LinDecay
                            +polity2 + area_1000 + lmtnest + ELF 
                            +strata(event_num) + cluster(ccode),data=TC4or5, method="efron") 

TC6orMore.PWP.LinDecay <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + GoodEnd_LinDecay + BadEnd_LinDecay
                               +polity2 + area_1000 + lmtnest + ELF 
                               +strata(event_num) + cluster(ccode),data=TC6orMore, method="efron") 

stargazer(TCAll.PWP.LinDecay)

##################################

model1 <- data.frame(Variable = rownames(TCAll.PWP.LinDecay$coefficients),
                     Coefficient = TCAll.PWP.LinDecay$coefficients,
                     SE = TCAll.PWP.LinDecay$)

coefs1 <- tidy(TCAll.PWP.LinDecay, conf.int = TRUE, exponentiate = T)
coefs1$HazardRatio <- coefs1$estimate
coefs1$Variable <- c("Peace", "Fighting", "Favorable for TC", "Favorable for SS", "Polity", "Area", "Terrain", "ELF")

coefs2 <- tidy(TC.1.PWP, conf.int = TRUE, exponentiate = T)
coefs2$HazardRatio <- coefs2$estimate
coefs2$Variable <- c( "Polity", "Area", "Terrain", "ELF")

coefs3 <- tidy(TC2or3.PWP.LinDecay, conf.int = TRUE, exponentiate = T)
coefs3$HazardRatio <- coefs3$estimate
coefs3$Variable <- c("Peace", "Fighting", "Favorable for TC", "Favorable for SS", "Polity", "Area", "Terrain", "ELF")

coefs4 <- tidy(TC4or5.PWP.LinDecay, conf.int = TRUE, exponentiate = T)
coefs4$HazardRatio <- coefs4$estimate
coefs4$Variable <- c("Peace", "Fighting", "Favorable for TC", "Favorable for SS", "Polity", "Area", "Terrain", "ELF")

coefs5 <- tidy(TC6orMore.PWP.LinDecay, conf.int = TRUE, exponentiate = T)
coefs5$HazardRatio <- coefs5$estimate
coefs5$Variable <- c("Peace", "Fighting", "Favorable for TC", "Favorable for SS", "Polity", "Area", "Terrain", "ELF")


coefs$HazardRatio <- coefs$estimate
ggcoef(coefs, mapping = aes_string(y = "Variable", x = "HazardRatio"), color = "orange", vline_intercept = 1)

dat <- filter(TCDat, last_year == 1)
table(dat$TC_tally)
hist(dat$TC_tally)

combined <- data.frame(Variable = coefs1$Variable,
                   HazardRatio = coefs1$HazardRatio,
                   SE = summary(TCAll.PWP.LinDecay)$coef[,3],
                   TC = "Combined")
combined$Variable <- factor(combined$Variable, levels=unique(as.character(combined$Variable)) )

first <- data.frame(Variable = coefs2$Variable,
                       HazardRatio = coefs2$HazardRatio,
                       SE = summary(TC.1.PWP)$coef[,3],
                       TC = "First TC")

mod2 <- data.frame(Variable = coefs3$Variable,
                       HazardRatio = coefs3$HazardRatio,
                       SE = summary(TC2or3.PWP.LinDecay)$coef[,3],
                       TC = "2-3")

mod3 <- data.frame(Variable = coefs4$Variable,
                       HazardRatio = coefs4$HazardRatio,
                       SE = summary(TC4or5.PWP.LinDecay)$coef[,3],
                       TC = "4-5")

mod4 <- data.frame(Variable = coefs5$Variable,
                       HazardRatio = coefs5$HazardRatio,
                       SE = summary(TC6orMore.PWP.LinDecay)$coef[,3],
                       TC = "6+")



allmod <- data.frame(rbind(mod4, mod3, mod2, first))
allmod$Variable <- factor(allmod$Variable, levels=unique(as.character(allmod$Variable)) )


interval1 <- -qnorm((1-0.9)/2)  # 90% multiplier
interval2 <- -qnorm((1-0.95)/2)  # 95% multiplier

zp2 <- ggplot(combined)
zp2 <- zp2 + geom_hline(yintercept = 1, colour = gray(1/2), lty = 2)
zp2 <- zp2 + geom_linerange(aes(x = Variable, ymin = HazardRatio - SE*interval1,
                                ymax = HazardRatio + SE*interval1),
                            lwd = 1, position = position_dodge(width = 1/2))
zp2 <- zp2 + geom_pointrange(aes(x = Variable, y = HazardRatio, ymin = HazardRatio - SE*interval2,
                                 ymax = HazardRatio + SE*interval2),
                             lwd = 1/2, position = position_dodge(width = 1/2),
                             shape = 21, fill = "WHITE")
zp2 <- zp2 + xlim(rev(levels(combined$Variable))) + coord_flip() + theme_bw()
#zp1 <- zp1 + ggtitle("Comparing several models")
print(zp2)  # The trick to these is position_dodge().

zp1 <- ggplot(allmod, aes(colour = TC))
zp1 <- zp1 + geom_hline(yintercept = 1, colour = gray(1/2), lty = 2)
zp1 <- zp1 + geom_linerange(aes(x = Variable, ymin = HazardRatio - SE*interval1,
                                ymax = HazardRatio + SE*interval1),
                            lwd = 1, position = position_dodge(width = 1/2))
zp1 <- zp1 + geom_pointrange(aes(x = Variable, y = HazardRatio, ymin = HazardRatio - SE*interval2,
                                 ymax = HazardRatio + SE*interval2),
                             lwd = 1/2, position = position_dodge(width = 1/2),
                             shape = 21, fill = "WHITE")
zp1 <- zp1 + xlim(rev(levels(allmod$Variable))) + coord_flip() + theme_bw()
#zp1 <- zp1 + ggtitle("Comparing several models")
print(zp1)  # The trick to these is position_dodge().
