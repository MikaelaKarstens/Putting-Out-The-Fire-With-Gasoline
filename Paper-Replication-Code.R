# Replication code for "Assessing State Making Strategies" by Douglas Lemke and Mikaela Karstens
#TEST
# Required packages for replication - Please install if needed (FIX THESE TO REMOVE UNUSED)

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

################ Experimenting. 

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

# Main Results Table

stargazer(TCAll.PWP, TC.1.PWP, TC2or3.PWP.LinDecay, TC4or5.PWP.LinDecay, TC6orMore.PWP.LinDecay,
          type = "latex", #out = "C:/Users/mjw65/Dropbox/State Making Strategies/State-Making-Strategies/table.html", 
          title = "PWP Gap Time Model Results",
          dep.var.labels = c("All TCs", "First TC", " TC 2 or 3", "TC 4 or 5", "TC 6+"),
          covariate.labels=c("Peace w/ TC", "Fighting w/ TC", "Forceful Reintegration", "Favorable Outcome for TC", "Polity 2 Score", "State Area (logged)","Mountains", "ELF"),
          keep.stat = c("n"))

# Additional models

# Absorbed by other TC as a control

TCAll.AbsorbedControl <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_LinDecay + GoodEnd_LinDecay + Absorbed_LinDecay
                  +polity2 + area_1000_log + lmtnest + ELF 
                  +strata(event_num) + cluster(ccode),data=TCDat, method="efron")

TC.1.AbsorbedControl <- coxph(TC1.Gap ~ +polity2 + area_1000_log + lmtnest + ELF 
                  +strata(event_num) + cluster(ccode),data=TC1, method="efron")

TC2or3.AbsorbedControl <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_LinDecay + GoodEnd_LinDecay + Absorbed_LinDecay
                            +polity2 + area_1000_log + lmtnest + ELF 
                            +strata(event_num) + cluster(ccode),data=TC2or3, method="efron") 

TC4or5.AbsorbedControl <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_LinDecay + GoodEnd_LinDecay + Absorbed_LinDecay
                            +polity2 + area_1000_log + lmtnest + ELF 
                            +strata(event_num) + cluster(ccode),data=TC4or5, method="efron") 

TC6orMore.AbsorbedControl <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_LinDecay + GoodEnd_LinDecay + Absorbed_LinDecay
                               +polity2 + area_1000_log + lmtnest + ELF 
                               +strata(event_num) + cluster(ccode),data=TC6orMore, method="efron") 

summary(TCAll.AbsorbedControl)

stargazer(TCAll.AbsorbedControl, TC.1.AbsorbedControl, TC2or3.AbsorbedControl, TC4or5.AbsorbedControl, TC6orMore.AbsorbedControl,
          type = "latex", #out = "C:/Users/mjw65/Dropbox/State Making Strategies/State-Making-Strategies/table.html", 
          title = "PWP Gap Time Model Results - Absorbed as a Control",
          dep.var.labels = c("All TCs", "First TC", " TC 2 or 3", "TC 4 or 5", "TC 6+"),
          covariate.labels=c("Peace w/ TC", "Fighting w/ TC", "Forceful Reintegration", "Favorable Outcome for TC", "Absorbed", "Polity 2 Score", "State Area (logged)","Mountains", "ELF"),
          keep.stat = c("n"))

# All End Types Disaggregated

TCAll.Disag <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_LinDecay + Promoted_LinDecay + Peaceful_LinDecay + Absorbed_LinDecay
                              +polity2 + area_1000_log + lmtnest + ELF 
                              +strata(event_num) + cluster(ccode),data=TCDat, method="efron")

TC.1.Disag <- coxph(TC1.Gap ~ +polity2 + area_1000_log + lmtnest + ELF 
                              +strata(event_num) + cluster(ccode),data=TC1, method="efron")

TC2or3.Disag <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_LinDecay + Promoted_LinDecay + Peaceful_LinDecay + Absorbed_LinDecay
                               +polity2 + area_1000_log + lmtnest + ELF 
                               +strata(event_num) + cluster(ccode),data=TC2or3, method="efron") 

TC4or5.Disag <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_LinDecay + Promoted_LinDecay + Peaceful_LinDecay + Absorbed_LinDecay
                               +polity2 + area_1000_log + lmtnest + ELF 
                               +strata(event_num) + cluster(ccode),data=TC4or5, method="efron") 

TC6orMore.Disag <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_LinDecay + Promoted_LinDecay + Peaceful_LinDecay + Absorbed_LinDecay
                                  +polity2 + area_1000_log + lmtnest + ELF 
                                  +strata(event_num) + cluster(ccode),data=TC6orMore, method="efron") 

summary(TCAll.Disag)

stargazer(TCAll.Disag, TC.1.Disag, TC2or3.Disag, TC4or5.Disag, TC6orMore.Disag,
          type = "latex", #out = "C:/Users/mjw65/Dropbox/State Making Strategies/State-Making-Strategies/table.html", 
          title = "PWP Gap Time Model Results - Disaggregated",
          dep.var.labels = c("All TCs", "First TC", " TC 2 or 3", "TC 4 or 5", "TC 6+"),
          covariate.labels=c("Peace w/ TC", "Fighting w/ TC", "Forceful Reintegration", "Promoted", "Peaceful Reintegration", "Absorbed", "Polity 2 Score", "State Area (logged)","Mountains", "ELF"),
          keep.stat = c("n"))

# Exponential Decay 20 year half life

TCAll.Exponential20 <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_ExpDecay_20 + GoodEnd_ExpDecay_20
                  +polity2 + area_1000_log + lmtnest + ELF 
                  +strata(event_num) + cluster(ccode),data=TCDat, method="efron")

TC.1.Exponential20 <- coxph(TC1.Gap ~ +polity2 + area_1000_log + lmtnest + ELF 
                  +strata(event_num) + cluster(ccode),data=TC1, method="efron")

TC2or3.Exponential20 <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_ExpDecay_20 + GoodEnd_ExpDecay_20 
                            +polity2 + area_1000_log + lmtnest + ELF 
                            +strata(event_num) + cluster(ccode),data=TC2or3, method="efron") 

TC4or5.Exponential20 <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_ExpDecay_20 + GoodEnd_ExpDecay_20 
                            +polity2 + area_1000_log + lmtnest + ELF 
                            +strata(event_num) + cluster(ccode),data=TC4or5, method="efron") 

TC6orMore.Exponential20 <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_ExpDecay_20 + GoodEnd_ExpDecay_20 
                               +polity2 + area_1000_log + lmtnest + ELF 
                               +strata(event_num) + cluster(ccode),data=TC6orMore, method="efron") 

summary(TCAll.Exponential20)


stargazer(TCAll.Exponential20, TC.1.Exponential20, TC2or3.Exponential20, TC4or5.Exponential20, TC6orMore.Exponential20,
          type = "latex", #out = "C:/Users/mjw65/Dropbox/State Making Strategies/State-Making-Strategies/table.html", 
          title = "PWP Gap Time Model Results - Exponential Decay with 20 year Half Life",
          dep.var.labels = c("All TCs", "First TC", " TC 2 or 3", "TC 4 or 5", "TC 6+"),
          covariate.labels=c("Peace w/ TC", "Fighting w/ TC", "Forceful Reintegration", "Favorable Outcome for TC", "Polity 2 Score", "State Area (logged)","Mountains", "ELF"),
          keep.stat = c("n"))

# Exponential Decay 10 year half life

TCAll.Exponential10 <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_ExpDecay_10 + GoodEnd_ExpDecay_10
                            +polity2 + area_1000_log + lmtnest + ELF 
                            +strata(event_num) + cluster(ccode),data=TCDat, method="efron")

TC.1.Exponential10 <- coxph(TC1.Gap ~ +polity2 + area_1000_log + lmtnest + ELF 
                            +strata(event_num) + cluster(ccode),data=TC1, method="efron")

TC2or3.Exponential10 <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_ExpDecay_10 + GoodEnd_ExpDecay_10 
                             +polity2 + area_1000_log + lmtnest + ELF 
                             +strata(event_num) + cluster(ccode),data=TC2or3, method="efron") 

TC4or5.Exponential10 <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_ExpDecay_10 + GoodEnd_ExpDecay_10 
                             +polity2 + area_1000_log + lmtnest + ELF 
                             +strata(event_num) + cluster(ccode),data=TC4or5, method="efron") 

TC6orMore.Exponential10 <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_ExpDecay_10 + GoodEnd_ExpDecay_10 
                                +polity2 + area_1000_log + lmtnest + ELF 
                                +strata(event_num) + cluster(ccode),data=TC6orMore, method="efron") 

summary(TCAll.Exponential10)


stargazer(TCAll.Exponential10, TC.1.Exponential10, TC2or3.Exponential10, TC4or5.Exponential10, TC6orMore.Exponential10,
          type = "latex", #out = "C:/Users/mjw65/Dropbox/State Making Strategies/State-Making-Strategies/table.html", 
          title = "PWP Gap Time Model Results - Exponential Decay with 10 year Half Life",
          dep.var.labels = c("All TCs", "First TC", " TC 2 or 3", "TC 4 or 5", "TC 6+"),
          covariate.labels=c("Peace w/ TC", "Fighting w/ TC", "Forceful Reintegration", "Favorable Outcome for TC", "Polity 2 Score", "State Area (logged)","Mountains", "ELF"),
          keep.stat = c("n"))

# 20 Year Binary

TCAll.t20 <-coxph(TCDat.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_t20 + GoodEnd_t20
                            +polity2 + area_1000_log + lmtnest + ELF 
                            +strata(event_num) + cluster(ccode),data=TCDat, method="efron")

TC.1.t20 <- coxph(TC1.Gap ~ +polity2 + area_1000_log + lmtnest + ELF 
                            +strata(event_num) + cluster(ccode),data=TC1, method="efron")

TC2or3.t20 <-coxph(TC2or3.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_t20 + GoodEnd_t20 
                             +polity2 + area_1000_log + lmtnest + ELF 
                             +strata(event_num) + cluster(ccode),data=TC2or3, method="efron") 

TC4or5.t20 <-coxph(TC4or5.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_t20 + GoodEnd_t20 
                             +polity2 + area_1000_log + lmtnest + ELF 
                             +strata(event_num) + cluster(ccode),data=TC4or5, method="efron") 

TC6orMore.t20 <-coxph(TC6orMore.Gap ~ Peace_wTC + Fighting_wTC  + Forceful_t20 + GoodEnd_t20 
                                +polity2 + area_1000_log + lmtnest + ELF 
                                +strata(event_num) + cluster(ccode),data=TC6orMore, method="efron") 

summary(TCAll.t20)


stargazer(TCAll.t20, TC.1.t20, TC2or3.t20, TC4or5.t20, TC6orMore.t20,
          type = "latex", #out = "C:/Users/mjw65/Dropbox/State Making Strategies/State-Making-Strategies/table.html", 
          title = "PWP Gap Time Model Results - 20 Year Binary Indicator",
          dep.var.labels = c("All TCs", "First TC", " TC 2 or 3", "TC 4 or 5", "TC 6+"),
          covariate.labels=c("Peace w/ TC", "Fighting w/ TC", "Forceful Reintegration", "Favorable Outcome for TC", "Polity 2 Score", "State Area (logged)","Mountains", "ELF"),
          keep.stat = c("n"))

# Hurdle Model - This seems so wrong. Its just lumping everything together. 

hurdlemod <- hurdle(TC_tally ~ Peace_wTC + Fighting_wTC  + Forceful_LinDecay + GoodEnd_LinDecay 
                    +polity2 + area_1000_log + lmtnest + ELF | polity2 + area_1000_log + lmtnest + ELF , data = TCDat)
summary(hurdlemod)

# Censored Probit
