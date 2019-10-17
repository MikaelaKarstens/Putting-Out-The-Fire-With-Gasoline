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

options(scipen = 999) # no scientific notation

# setwd("C:/Users/mjw65/Dropbox/State Making Strategies/State-Making-Strategies") #209 Pond Lab

TCDat <- read.csv("Survival-Data-TC.csv")
names(TCDat)
