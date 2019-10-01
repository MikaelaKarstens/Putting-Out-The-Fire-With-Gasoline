# Repeat Failure Model

library(foreign) # download .dta files for Stata < 13
library(haven) # download .dta files for Stata 13+
library(sampleSelection) # Package containing heckit for heckprob substitute

data <- read_dta("Deterring_TC.dta") #Uses haven

names(data) #print variable names to make sure everything read in alright.

