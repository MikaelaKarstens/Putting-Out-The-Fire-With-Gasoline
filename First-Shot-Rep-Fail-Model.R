# Repeat Failure Model

library(foreign) # download .dta files for Stata < 13
library(haven) # download .dta files for Stata 13+
library(sampleSelection) # Package containing heckit for heckprob substitute
library(dplyr) #Used for dataset laterations and survival time calculations
library(lubridate) #Used for survival data formatting
library(survival)
library(reReg) #recurrent events package

data <- read_dta("Deterring_TC.dta") #Uses haven

names(data) #print variable names to make sure everything read in alright.

SurvDat <- read.csv("Repeat-Failure-Data-10-3.csv")

reObj <- reSurv(SurvDat$t.start, SurvDat$t.stop, SurvDat$ccode, SurvDat$event)
plotEvents(reObj)

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


