# Replication of Doug's initial code

library(foreign) # download .dta files for Stata < 13
library(haven) # download .dta files for Stata 13+

data <- read_dta("Deterring_TC.dta") #Uses haven

names(data) #print variable names to make sure everything read in alright.

