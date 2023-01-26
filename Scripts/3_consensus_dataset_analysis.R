##########################
#### user: Lisbeth Hordley
#### date: June 2022
#### info: Rates of establishment over time (consensus dataset)

# Packages
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(AER)
library(msm)
library(segmented)
library(ggeffects)
library(DHARMa)
options(scipen=999)

# Load data
cdata <- read.csv("Coloniser data/moth_coloniser_nspp_consensus.csv", header=TRUE)

#########################################################
################# 1. Entire time series ################# 
#########################################################

# Run GLMs to see how trend of all species and species groups change over time
# Try this for whole time period and compare to segmented regression - which fit better?
# This will tell us whether the trend differs over time, or stays the same
# And whether the trend is significant

### 1a. All species

# Poisson GLM
all_spp_m <- glm(all_nspp ~ date_num, data=cdata, family="poisson")
summary(all_spp_m) 

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = all_spp_m)
plot(simulationOutput) 
testDispersion(simulationOutput)
# no assumptions violated

# Segmented GLM
all_spp_seg <- segmented(all_spp_m,
                         seg.Z = ~date_num)

summary(all_spp_seg) 
# Davies' test - tests for the hypothesis that leftSlope=Rightslop (no difference in the two slopes - breakpoint does not exist)
davies.test(all_spp_m, seg.Z= ~date_num) # non-significant - accept null that breakpoint does not exist
slope(all_spp_seg)

# Compare to model with no break - use BIC following Macgregor moth biomass paper 
BIC(all_spp_m) 
BIC(all_spp_seg) 
# model with no break fits best

# What is the incidence rate over time?
output <- coef(summary(all_spp_m))
all_spp_r <- data.frame(Startdate = 1905, Enddate=2015, # start date of model
                        rate = exp(output["date_num","Estimate"]), # exponential of the adventive slope
                        rate_SE = 
                          deltamethod(~exp(x1), 
                                      output["date_num","Estimate"], 
                                      output["date_num","Std. Error"]^2, ses=TRUE)) # standard error around exponential of immigrant slope

# Produce confidence intervals for the rate estimates
setDT(all_spp_r)
all_spp_r[, rate_lower := rate - 1.96*rate_SE]
all_spp_r[, rate_upper := rate + 1.96*rate_SE]
all_spp_r

#

### 1b. Immigrant species
# Poisson GLM
imm_spp_m <- glm(imm_nspp ~ date_num, data=cdata, family="poisson")
summary(imm_spp_m) 

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = imm_spp_m)
plot(simulationOutput)
testDispersion(simulationOutput)

# Segmented GLM
imm_spp_seg <- segmented(imm_spp_m,
                         seg.Z = ~date_num)

summary(imm_spp_seg) # break is at 4 (i.e. 1935)
# Davies' test - tests for the hypothesis that leftSlope=Rightslop (no difference in the two slopes - breakpoint does not exist)
davies.test(imm_spp_m, seg.Z= ~date_num) # non-significant - accept null that breakpoint does not exist
slope(imm_spp_seg)

# Compare to model with no break
BIC(imm_spp_m) 
BIC(imm_spp_seg) 
# within 2 points so take the simplest model (non-segmented)

# Plot trend over time with raw data
pred <- ggpredict(imm_spp_m, terms="date_num")
pred$mid_date <- cdata$mid_date

imm_spp_p <- ggplot(pred, aes(x=mid_date, y=predicted))+
  geom_line()+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.1)+
  geom_point(data=cdata, aes(x=mid_date, y=imm_nspp))+
  scale_x_continuous(breaks=seq(1905,2015,by=20))+
  labs(x="Year", y="Number of species")+
  scale_y_continuous(breaks=seq(0,30,by=5), limits=c(0,30))+
  theme_classic()
imm_spp_p

# What is the incidence rate over time? (from non-segmented model)
output <- coef(summary(imm_spp_m))
imm_spp_r <- data.frame(Startdate = 1905, Enddate=2015, # start date of model
                        rate = exp(output["date_num","Estimate"]), # exponential of the adventive slope
                        rate_SE = 
                          deltamethod(~exp(x1), 
                                      output["date_num","Estimate"], 
                                      output["date_num","Std. Error"]^2, ses=TRUE)) # standard error around exponential of immigrant slope

# Produce confidence intervals for the rate estimates
setDT(imm_spp_r)
imm_spp_r[, rate_lower := rate - 1.96*rate_SE]
imm_spp_r[, rate_upper := rate + 1.96*rate_SE]
imm_spp_r

#

### 1c. Adventive species
# Poisson GLM
adv_spp_m <- glm(adv_nspp ~ date_num, data=cdata, family="poisson")
summary(adv_spp_m) # date is significant

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = adv_spp_m)
plot(simulationOutput) 
plotResiduals(simulationOutput, cdata$date_num, quantreg = T) 
res = recalculateResiduals(simulationOutput, group = cdata$date_num) # recalculate residuals to aggregate residuals per time step
testTemporalAutocorrelation(res, time = unique(cdata$date_num)) # non-sig
testDispersion(simulationOutput)
# no assumptions violated

# Segmented GLM
adv_spp_seg <- segmented(adv_spp_m,
                         seg.Z = ~date_num)

summary(adv_spp_seg) # break is at 8.6 (i.e. ~1975)
# Davies' test - tests for the hypothesis that leftSlope=Rightslop (no difference in the two slopes - breakpoint does not exist)
davies.test(adv_spp_m, seg.Z= ~date_num) # non-significant - accept null that breakpoint does not exist
slope(adv_spp_seg)

# Compare to model with no break
BIC(adv_spp_m) 
BIC(adv_spp_seg) 
# non-segmented model has lower BIC

# Plot trend over time with raw data
pred <- ggpredict(adv_spp_m, terms="date_num")
pred$mid_date <- cdata$mid_date

adv_spp_p <- ggplot(pred, aes(x=mid_date, y=predicted))+
  geom_line()+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.1)+
  geom_point(data=cdata, aes(x=mid_date, y=adv_nspp))+
  scale_x_continuous(breaks=seq(1905,2015,by=20))+
  labs(x="Year", y="Number of species")+
  scale_y_continuous(breaks=seq(0,30,by=5), limits=c(0,30))+
  theme_classic()
adv_spp_p

# What is the incidence rate over time?
output <- coef(summary(adv_spp_m))
adv_spp_r <- data.frame(Startdate = 1905, Enddate=2015, # start date of model
                        rate = exp(output["date_num","Estimate"]), # exponential of the adventive slope
                        rate_SE = 
                          deltamethod(~exp(x1), 
                                      output["date_num","Estimate"], 
                                      output["date_num","Std. Error"]^2, ses=TRUE)) # standard error around exponential of immigrant slope

# Produce confidence intervals for the rate estimates
setDT(adv_spp_r)
adv_spp_r[, rate_lower := rate - 1.96*rate_SE]
adv_spp_r[, rate_upper := rate + 1.96*rate_SE]
adv_spp_r

#

### 1d. Immigrant + native host plants
# Poisson GLM
imm_nat_spp_m <- glm(imm_nat_nspp ~ date_num, data=cdata, family="poisson")
summary(imm_nat_spp_m) # date is significant

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = imm_nat_spp_m)
plot(simulationOutput)
plotResiduals(simulationOutput, cdata$date_num, quantreg = T) 
res = recalculateResiduals(simulationOutput, group = cdata$date_num) # recalculate residuals to aggregate residuals per time step
testTemporalAutocorrelation(res, time = unique(cdata$date_num)) # non-sig
testDispersion(simulationOutput)
# no assumptions violated

# Segmented GLM
imm_nat_spp_seg <- segmented(imm_nat_spp_m,
                             seg.Z = ~date_num)

summary(imm_nat_spp_seg) # coefficients are far too high - segmented model is struggling with too many zeros
# just stick with non-segmented model
BIC(imm_nat_spp_m) # 59.65

# Plot trend over time with raw data
pred <- ggpredict(imm_nat_spp_m, terms="date_num")
pred$mid_date <- cdata$mid_date

imm_nat_spp_p <- ggplot(pred, aes(x=mid_date, y=predicted))+
  geom_line()+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.1)+
  geom_point(data=cdata, aes(x=mid_date, y=imm_nat_nspp))+
  scale_x_continuous(breaks=seq(1905,2015,by=20))+
  labs(x="Year", y="Number of species")+
  theme_classic()
imm_nat_spp_p

# What is the incidence rate over time?
output <- coef(summary(imm_nat_spp_m))
imm_nat_spp_r <- data.frame(Startdate = 1905, Enddate=2015, # start date of model
                            rate = exp(output["date_num","Estimate"]), # exponential of the adventive slope
                            rate_SE = 
                              deltamethod(~exp(x1), 
                                          output["date_num","Estimate"], 
                                          output["date_num","Std. Error"]^2, ses=TRUE)) # standard error around exponential of immigrant slope

# Produce confidence intervals for the rate estimates
setDT(imm_nat_spp_r)
imm_nat_spp_r[, rate_lower := rate - 1.96*rate_SE]
imm_nat_spp_r[, rate_upper := rate + 1.96*rate_SE]
imm_nat_spp_r
# ~20% significant increase in the number of species colonising per decade on average

#

### 1e. Immigrant + non-native host plants
# Poisson GLM
imm_nn_spp_m <- glm(imm_nn_nspp ~ date_num, data=cdata, family="poisson")
summary(imm_nn_spp_m) # date is non-significant

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = imm_nn_spp_m)
plot(simulationOutput) 
res = recalculateResiduals(simulationOutput, group = cdata$date_num) # recalculate residuals to aggregate residuals per time step
testTemporalAutocorrelation(res, time = unique(cdata$date_num)) # non-sig 
testDispersion(simulationOutput)

# Segmented GLM
imm_nn_spp_seg <- segmented(imm_nn_spp_m,
                            seg.Z = ~date_num)

summary(imm_nn_spp_seg) # break is at 4.8 (i.e. ~1945)
slope(imm_nn_spp_seg) # both slopes are non-significant - no significant trend over time
davies.test(imm_nn_spp_m,
            seg.Z = ~date_num)

BIC(imm_nn_spp_m)
BIC(imm_nn_spp_seg) 
# within 2 points so take the simplest model (non-segmented)

# Plot trend over time with raw data - dashed line and no error bars to show non-significance
pred <- ggpredict(imm_nn_spp_m, terms="date_num")
pred$mid_date <- imm_nn_spp$mid_date

imm_nn_spp_p <- ggplot(pred, aes(x=mid_date, y=predicted))+
  geom_line(linetype=2)+
  #geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.1)+
  geom_point(data=imm_nn_spp, aes(x=mid_date, y=nspp))+
  scale_x_continuous(breaks=seq(1905,2015,by=20))+
  labs(x="Year", y="Number of species")+
  theme_classic()
imm_nn_spp_p

# What is the incidence rate over time?
output <- coef(summary(imm_nn_spp_m))
imm_nn_spp_r <- data.frame(Startdate = 1905, Enddate=2015, # start date of model
                           rate = exp(output["date_num","Estimate"]), # exponential of the adventive slope
                           rate_SE = 
                             deltamethod(~exp(x1), 
                                         output["date_num","Estimate"], 
                                         output["date_num","Std. Error"]^2, ses=TRUE)) # standard error around exponential of immigrant slope

# Produce confidence intervals for the rate estimates
setDT(imm_nn_spp_r)
imm_nn_spp_r[, rate_lower := rate - 1.96*rate_SE]
imm_nn_spp_r[, rate_upper := rate + 1.96*rate_SE]
imm_nn_spp_r

#

### 1f. Adventives + native host plants
# Poisson GLM
adv_nat_spp_m <- glm(adv_nat_nspp ~ date_num, data=cdata, family="poisson")
summary(adv_nat_spp_m) # date is significant

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = adv_nat_spp_m)
plot(simulationOutput)
res = recalculateResiduals(simulationOutput, group = cdata$date_num) # recalculate residuals to aggregate residuals per time step
testTemporalAutocorrelation(res, time = unique(cdata$date_num)) # non-sig
testDispersion(simulationOutput)
# no assumptions violated

# Segmented GLM
adv_nat_spp_seg <- segmented(adv_nat_spp_m,
                             seg.Z = ~date_num)

summary(adv_nat_spp_seg) # break is at 8 (i.e. ~1975)
# Davies' test - tests for the hypothesis that leftSlope=Rightslop (no difference in the two slopes - breakpoint does not exist)
davies.test(adv_nat_spp_m, seg.Z= ~date_num) # non-significant - accept null that breakpoint does not exist
slope(adv_nat_spp_seg)

# Compare to model with no break
BIC(adv_nat_spp_m) 
BIC(adv_nat_spp_seg) 
# within 2 points so take the simplest model (non-segmented)

# Plot trend over time with raw data
pred <- ggpredict(adv_nat_spp_m, terms="date_num")
pred$mid_date <- adv_nat_spp$mid_date

adv_nat_spp_p <- ggplot(pred, aes(x=mid_date, y=predicted))+
  geom_line()+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.1)+
  geom_point(data=adv_nat_spp, aes(x=mid_date, y=nspp))+
  scale_x_continuous(breaks=seq(1905,2015,by=20))+
  labs(x="Year", y="Number of species")+
  theme_classic()
adv_nat_spp_p

# What is the incidence rate over time?
output <- coef(summary(adv_nat_spp_m))
adv_nat_spp_r <- data.frame(Startdate = 1905, Enddate=2015, # start date of model
                            rate = exp(output["date_num","Estimate"]), # exponential of the adventive slope
                            rate_SE = 
                              deltamethod(~exp(x1), 
                                          output["date_num","Estimate"], 
                                          output["date_num","Std. Error"]^2, ses=TRUE)) # standard error around exponential of immigrant slope

# Produce confidence intervals for the rate estimates
setDT(adv_nat_spp_r)
adv_nat_spp_r[, rate_lower := rate - 1.96*rate_SE]
adv_nat_spp_r[, rate_upper := rate + 1.96*rate_SE]
adv_nat_spp_r
# ~24% significant increase in the number of species colonising per decade on average

#

### 1g. Adventive + non-native host plants
# Poisson GLM
adv_nn_spp_m <- glm(adv_nn_nspp ~ date_num, family="poisson", data=cdata)
summary(adv_nn_spp_m) # date is significant

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = adv_nn_spp_m)
plot(simulationOutput)
res = recalculateResiduals(simulationOutput, group = cdata$date_num) # recalculate residuals to aggregate residuals per time step
testTemporalAutocorrelation(res, time = unique(cdata$date_num)) # non-sig
testDispersion(simulationOutput)

# Segmented GLM
adv_nn_spp_seg <- segmented(adv_nn_spp_m,
                            seg.Z = ~date_num)

summary(adv_nn_spp_seg) # coefficients are far too high - segmented model is struggling with too many zeros
# just stick with non-segmented model
BIC(adv_nn_spp_m) 

# Plot trend over time with raw data
pred <- ggpredict(adv_nn_spp_m, terms="date_num")
pred$mid_date <- cdata$mid_date

adv_nn_spp_p <- ggplot(pred, aes(x=mid_date, y=predicted))+
  geom_line()+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.1)+
  geom_point(data=cdata, aes(x=mid_date, y=adv_nn_nspp))+
  scale_x_continuous(breaks=seq(1905,2015,by=20))+
  labs(x="Year", y="Number of species")+
  theme_classic()
adv_nn_spp_p

# What is the incidence rate over time?
output <- coef(summary(adv_nn_spp_m))
adv_nn_spp_r <- data.frame(Startdate = 1905, Enddate=2015, # start date of model
                           rate = exp(output["date_num","Estimate"]), # exponential of the adventive slope
                           rate_SE = 
                             deltamethod(~exp(x1), 
                                         output["date_num","Estimate"], 
                                         output["date_num","Std. Error"]^2, ses=TRUE)) # standard error around exponential of immigrant slope

# Produce confidence intervals for the rate estimates
setDT(adv_nn_spp_r)
adv_nn_spp_r[, rate_lower := rate - 1.96*rate_SE]
adv_nn_spp_r[, rate_upper := rate + 1.96*rate_SE]
adv_nn_spp_r



######################################################################
################# 2. Last decade (2010-2019) removed ################# 
######################################################################

# Remove last decade
cdata <- cdata[!cdata$mid_date=="2015",]

# Run GLMs to see how trend of all species and species groups change over time
# Try this for whole time period and compare to segmented regression - which fit better?
# This will tell us whether the trend differs over time, or stays the same
# And whether the trend is significant

### 2a. All species

# Poisson GLM
all_spp_m <- glm(all_nspp ~ date_num, data=cdata, family="poisson")
summary(all_spp_m) # date is significant

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = all_spp_m)
plot(simulationOutput) 
testDispersion(simulationOutput)
# no assumptions violated

# Segmented GLM
all_spp_seg <- segmented(all_spp_m,
                         seg.Z = ~date_num)

summary(all_spp_seg) 
# Davies' test - tests for the hypothesis that leftSlope=Rightslop (no difference in the two slopes - breakpoint does not exist)
davies.test(all_spp_m, seg.Z= ~date_num) # non-significant - accept null that breakpoint does not exist
slope(all_spp_seg)

# Compare to model with no break - use BIC following Macgregor moth biomass paper 
BIC(all_spp_m) 
BIC(all_spp_seg) 
# model with no break fits best

# What is the incidence rate over time?
output <- coef(summary(all_spp_m))
all_spp_r <- data.frame(Startdate = 1905, Enddate=2015, # start date of model
                        rate = exp(output["date_num","Estimate"]), # exponential of the adventive slope
                        rate_SE = 
                          deltamethod(~exp(x1), 
                                      output["date_num","Estimate"], 
                                      output["date_num","Std. Error"]^2, ses=TRUE)) # standard error around exponential of immigrant slope

# Produce confidence intervals for the rate estimates
setDT(all_spp_r)
all_spp_r[, rate_lower := rate - 1.96*rate_SE]
all_spp_r[, rate_upper := rate + 1.96*rate_SE]
all_spp_r

#

### 2b. Immigrant species
# Poisson GLM
imm_spp_m <- glm(imm_nspp ~ date_num, data=cdata, family="poisson")
summary(imm_spp_m) # date is significant

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = imm_spp_m)
plot(simulationOutput) 
testDispersion(simulationOutput)
# no assumptions violated

# Segmented GLM
imm_spp_seg <- segmented(imm_spp_m,
                         seg.Z = ~date_num)

summary(imm_spp_seg) 
# Davies' test - tests for the hypothesis that leftSlope=Rightslop (no difference in the two slopes - breakpoint does not exist)
davies.test(imm_spp_m, seg.Z= ~date_num) # non-significant - accept null that breakpoint does not exist
slope(imm_spp_seg)

# Compare to model with no break
BIC(imm_spp_m) 
BIC(imm_spp_seg) 
# within 2 points
# take the simplest model - non-segmented

# What is the incidence rate over time? (from non-segmented model)
output <- coef(summary(imm_spp_m))
imm_spp_r <- data.frame(Startdate = 1905, Enddate=2015, # start date of model
                        rate = exp(output["date_num","Estimate"]), # exponential of the adventive slope
                        rate_SE = 
                          deltamethod(~exp(x1), 
                                      output["date_num","Estimate"], 
                                      output["date_num","Std. Error"]^2, ses=TRUE)) # standard error around exponential of immigrant slope

# Produce confidence intervals for the rate estimates
setDT(imm_spp_r)
imm_spp_r[, rate_lower := rate - 1.96*rate_SE]
imm_spp_r[, rate_upper := rate + 1.96*rate_SE]
imm_spp_r

#

### 2c. Adventive species
# Poisson GLM
adv_spp_m <- glm(adv_nspp ~ date_num, data=cdata, family="poisson")
summary(adv_spp_m) # date is significant

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = adv_spp_m)
plot(simulationOutput) 
plotResiduals(simulationOutput, cdata$date_num, quantreg = T) 
res = recalculateResiduals(simulationOutput, group = cdata$date_num) # recalculate residuals to aggregate residuals per time step
testTemporalAutocorrelation(res, time = unique(cdata$date_num)) # non-sig
testDispersion(simulationOutput)
# no assumptions violated

# Segmented GLM
adv_spp_seg <- segmented(adv_spp_m,
                         seg.Z = ~date_num)

summary(adv_spp_seg) 
# Davies' test - tests for the hypothesis that leftSlope=Rightslop (no difference in the two slopes - breakpoint does not exist)
davies.test(adv_spp_m, seg.Z= ~date_num) # non-significant - accept null that breakpoint does not exist
slope(adv_spp_seg)

# Compare to model with no break
BIC(adv_spp_m) 
BIC(adv_spp_seg) 
# segmented model has lower BIC

# Note: graph of segmented model result plotted at end of script

#

### 2d. Immigrant + native host plants
# Poisson GLM
imm_nat_spp_m <- glm(imm_nat_nspp ~ date_num, data=cdata, family="poisson")
summary(imm_nat_spp_m) # date is significant

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = imm_nat_spp_m)
plot(simulationOutput)
res = recalculateResiduals(simulationOutput, group = cdata$date_num) # recalculate residuals to aggregate residuals per time step
testTemporalAutocorrelation(res, time = unique(cdata$date_num)) # non-sig
testDispersion(simulationOutput)
# no assumptions violated

# Segmented GLM
imm_nat_spp_seg <- segmented(imm_nat_spp_m,
                             seg.Z = ~date_num)

summary(imm_nat_spp_seg)# coefficients are far too high - segmented model is struggling with too many zeros
# just stick with non-segmented model
BIC(imm_nat_spp_m) 

# What is the incidence rate over time?
output <- coef(summary(imm_nat_spp_m))
imm_nat_spp_r <- data.frame(Startdate = 1905, Enddate=2015, # start date of model
                            rate = exp(output["date_num","Estimate"]), # exponential of the adventive slope
                            rate_SE = 
                              deltamethod(~exp(x1), 
                                          output["date_num","Estimate"], 
                                          output["date_num","Std. Error"]^2, ses=TRUE)) # standard error around exponential of immigrant slope

# Produce confidence intervals for the rate estimates
setDT(imm_nat_spp_r)
imm_nat_spp_r[, rate_lower := rate - 1.96*rate_SE]
imm_nat_spp_r[, rate_upper := rate + 1.96*rate_SE]
imm_nat_spp_r

#

### 2e. Immigrant + non-native host plants
# Poisson GLM
imm_nn_spp_m <- glm(imm_nn_nspp ~ date_num, data=cdata, family="poisson")
summary(imm_nn_spp_m) # date is non-significant

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = imm_nn_spp_m)
plot(simulationOutput) 
testDispersion(simulationOutput)

# Segmented GLM
imm_nn_spp_seg <- segmented(imm_nn_spp_m,
                            seg.Z = ~date_num)

summary(imm_nn_spp_seg) 
slope(imm_nn_spp_seg) 
davies.test(imm_nn_spp_m,
            seg.Z = ~date_num)

BIC(imm_nn_spp_m)
BIC(imm_nn_spp_seg)
# within 2 points so take the simplest model (non-segmented)

# What is the incidence rate over time?
output <- coef(summary(imm_nn_spp_m))
imm_nn_spp_r <- data.frame(Startdate = 1905, Enddate=2015, # start date of model
                           rate = exp(output["date_num","Estimate"]), # exponential of the adventive slope
                           rate_SE = 
                             deltamethod(~exp(x1), 
                                         output["date_num","Estimate"], 
                                         output["date_num","Std. Error"]^2, ses=TRUE)) # standard error around exponential of immigrant slope

# Produce confidence intervals for the rate estimates
setDT(imm_nn_spp_r)
imm_nn_spp_r[, rate_lower := rate - 1.96*rate_SE]
imm_nn_spp_r[, rate_upper := rate + 1.96*rate_SE]
imm_nn_spp_r

#

### 2f. Adventives + native host plants
# Poisson GLM
adv_nat_spp_m <- glm(adv_nat_nspp ~ date_num, data=cdata, family="poisson")
summary(adv_nat_spp_m) 

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = adv_nat_spp_m)
plot(simulationOutput)
res = recalculateResiduals(simulationOutput, group = cdata$date_num) # recalculate residuals to aggregate residuals per time step
testTemporalAutocorrelation(res, time = unique(cdata$date_num)) # non-sig
testDispersion(simulationOutput)
# no assumptions violated

# Segmented GLM
adv_nat_spp_seg <- segmented(adv_nat_spp_m,
                             seg.Z = ~date_num)

summary(adv_nat_spp_seg) 
# Davies' test - tests for the hypothesis that leftSlope=Rightslop (no difference in the two slopes - breakpoint does not exist)
davies.test(adv_nat_spp_m, seg.Z= ~date_num) 
slope(adv_nat_spp_seg)

# Compare to model with no break
BIC(adv_nat_spp_m) 
BIC(adv_nat_spp_seg)
# segmented model has lower BIC

# Note: graph of segmented model result plotted at end of script

#

### 2g. Adventive + non-native host plants
# Poisson GLM
adv_nn_spp_m <- glm(adv_nn_nspp ~ date_num, family="poisson", data=cdata)
summary(adv_nn_spp_m) 

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = adv_nn_spp_m)
plot(simulationOutput)
res = recalculateResiduals(simulationOutput, group = cdata$date_num) # recalculate residuals to aggregate residuals per time step
testTemporalAutocorrelation(res, time = unique(cdata$date_num)) # non-sig
testDispersion(simulationOutput)

# Segmented GLM
adv_nn_spp_seg <- segmented(adv_nn_spp_m,
                            seg.Z = ~date_num)

summary(adv_nn_spp_seg)# coefficients are far too high - segmented model is struggling with too many zeros
# just stick with non-segmented model
BIC(adv_nn_spp_m) 

# What is the incidence rate over time?
output <- coef(summary(adv_nn_spp_m))
adv_nn_spp_r <- data.frame(Startdate = 1905, Enddate=2015, # start date of model
                           rate = exp(output["date_num","Estimate"]), # exponential of the adventive slope
                           rate_SE = 
                             deltamethod(~exp(x1), 
                                         output["date_num","Estimate"], 
                                         output["date_num","Std. Error"]^2, ses=TRUE)) # standard error around exponential of immigrant slope

# Produce confidence intervals for the rate estimates
setDT(adv_nn_spp_r)
adv_nn_spp_r[, rate_lower := rate - 1.96*rate_SE]
adv_nn_spp_r[, rate_upper := rate + 1.96*rate_SE]
adv_nn_spp_r





####################################
############## GRAPHS ############## 
####################################


# Graph of segmented model: adventives
newdataesw <- expand.grid(date_num = cdata$date_num, nspp = 0)
dat <- predict(adv_spp_seg, newdata = newdataesw, type="response", se.fit=TRUE)
dat <- as.data.frame(dat)
dat$mid_date <- c(1905,1915,1925,1935,1945,1955,1965,1975,1985,1995,2005)
setDT(dat)
dat[, rate_lower := fit - 1.96*se.fit]
dat[, rate_upper := fit + 1.96*se.fit]

adv_spp_seg_p <- ggplot(dat, aes(x=mid_date, y=fit))+
  geom_line()+
  geom_ribbon(aes(ymin=rate_lower, ymax=rate_upper), alpha=0.1)+
  geom_point(data=cdata, aes(x=mid_date, y=adv_nspp))+
  scale_x_continuous(breaks=seq(1905,2015,by=20))+
  labs(x="Year", y="Number of species")+
  theme_classic()
adv_spp_seg_p # not convinced this is right but it's close to the plot below
ggsave(adv_spp_seg_p, file="Graphs/FigureS1a.png")

# Graph of segmented model: adventives + native host plants
newdataesw <- expand.grid(date_num = cdata$date_num, nspp = 0)
dat2 <- predict(adv_nat_spp_seg, newdata = newdataesw, type="response", se.fit=TRUE)
dat2 <- as.data.frame(dat2)
dat2$mid_date <- c(1905,1915,1925,1935,1945,1955,1965,1975,1985,1995,2005)
setDT(dat2)
dat2[, rate_lower := fit - 1.96*se.fit]
dat2[, rate_upper := fit + 1.96*se.fit]

adv_nat_spp_seg_p <- ggplot(dat2, aes(x=mid_date, y=fit))+
  geom_line()+
  geom_ribbon(aes(ymin=rate_lower, ymax=rate_upper), alpha=0.1)+
  geom_point(data=cdata, aes(x=mid_date, y=adv_nat_nspp))+
  scale_x_continuous(breaks=seq(1905,2015,by=20))+
  labs(x="Year", y="Number of species")+
  theme_classic()
adv_nat_spp_seg_p
ggsave(adv_nat_spp_seg_p, file="Graphs/FigureS1b.png")













