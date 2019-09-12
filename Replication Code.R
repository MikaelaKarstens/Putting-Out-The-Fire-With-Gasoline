# Replication of Doug's initial code

library(foreign) # download .dta files for Stata < 13
library(haven) # download .dta files for Stata 13+
library(sampleSelection) # Package containing heckit for heckprob substitute

data <- read_dta("Deterring_TC.dta") #Uses haven

names(data) #print variable names to make sure everything read in alright.

# Replicating heckprob results

# heckprob new_dfs_dummy Peace_WDFS Fighting_WDFS_alt Absorbed_t20 Forceful_t20 Peaceful_t20
# Promoted_t20 DFS_tally newdfs_spline newdfs_spline_sq newdfs_spline_cube 
# sel(first_dfs_alt = Devel polity2 area_1000 lmtnest ELF Country_Age 
#    Country_Age_Sq Country_Age_Cube) vce(cluster code)

Table2 <- heckit(selection = factor(first_TC_alt) ~ Devel + polity2 + area_1000 + lmtnest + ELF + Country_Age +
          Country_Age_Sq + Country_Age_Cube,
          outcome = factor(new_TC_dummy) ~ Peace_wTC + Fighting_wTC + Absorbed_t20 + Forceful_t20 +
          Peaceful_t20 + Promoted_t20 + TC_tally + new_TC_spline + newTC_spline_sq + newTC_spline_cube,
          data = data)

Table2ML <- selection(factor(first_TC_alt) ~ Devel + polity2 + area_1000 + lmtnest + ELF + Country_Age +
                     Country_Age_Sq + Country_Age_Cube,
                   factor(new_TC_dummy) ~ Peace_wTC + Fighting_wTC + Absorbed_t20 + Forceful_t20 +
                     Peaceful_t20 + Promoted_t20 + TC_tally + new_TC_spline + newTC_spline_sq + newTC_spline_cube,
                   data = data)

summary(Table2) 
summary(Table2ML)

# NONE OF THIS CONVERGES