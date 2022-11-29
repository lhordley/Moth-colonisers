##########################
#### user: Lisbeth Hordley
#### date: June 2022
#### info: Analysing changing drivers of moth colonisation

# Packages
library(data.table)
library(ggplot2)
library(DHARMa)
library(dplyr)
library(tidyverse)
library(AER)
library(msm)
library(segmented)
library(ggeffects)
#library(glmmTMB)
options(scipen=999)

setwd("~/Moth colonisers")

# Load data
cdata <- read.csv("../Moth colonisers/moth_coloniser_dataMP.csv", header=TRUE)

cdata$mid_date <- factor(cdata$mid_date, levels=c("1905","1915","1925","1935","1945","1955","1965","1975","1985","1995","2005", "2015"))
cdata$coloniser_mode <- factor(cdata$coloniser_mode, levels=c("A", "I"))
cdata$host_plant_status <- factor(cdata$host_plant_status, levels=c("Native", "Non-native"))

# create new datasets 
cdata <- cdata[!is.na(cdata$coloniser_mode),] # 144 species 
# need to remove two species now considered overlooked residents (Carpatolechia alburnella and Dendrolimus pini)
cdata <- cdata[!cdata$scientific_name=="Carpatolechia alburnella",]
cdata <- cdata[!cdata$scientific_name=="Dendrolimus pini",]
# 142 species now 

# # calc percentage of micros
# table(cdata$micro_macro) # 42 macros, 100 micros
# (101/142)*100 # 71%

cdata1 <- cdata %>% group_by(mid_date, .drop=FALSE) %>% summarise(all_nspp=n())
cdata2 <- cdata %>% filter(coloniser_mode=="I") %>% group_by(mid_date, .drop=FALSE) %>% summarise(imm_nspp=n())
cdata3 <- cdata %>% filter(coloniser_mode=="A") %>% group_by(mid_date, .drop=FALSE) %>% summarise(adv_nspp=n())
cdata4 <- cdata %>% filter(coloniser_mode=="I", host_plant_status=="Native") %>% group_by(mid_date, .drop=FALSE) %>% 
  summarise(imm_nat_nspp=n())
cdata5 <- cdata %>% filter(coloniser_mode=="I", host_plant_status=="Non-native") %>% group_by(mid_date, .drop=FALSE) %>% 
  summarise(imm_nn_nspp=n())
cdata6 <- cdata %>% filter(coloniser_mode=="A", host_plant_status=="Native") %>% group_by(mid_date, .drop=FALSE) %>% 
  summarise(adv_nat_nspp=n())
cdata7 <- cdata %>% filter(coloniser_mode=="A", host_plant_status=="Non-native") %>% group_by(mid_date, .drop=FALSE) %>% 
  summarise(adv_nn_nspp=n())

cdata_mp <- list(cdata1, cdata2, cdata3, cdata4, cdata5, cdata6, cdata7) %>% reduce(full_join, by = "mid_date")
cdata_mp <- as.data.frame(cdata_mp)
cdata_mp$date_num <- case_when(
  cdata_mp$mid_date==1905 ~ 1,
  cdata_mp$mid_date==1915 ~ 2,
  cdata_mp$mid_date==1925 ~ 3,
  cdata_mp$mid_date==1935 ~ 4,
  cdata_mp$mid_date==1945 ~ 5,
  cdata_mp$mid_date==1955 ~ 6,
  cdata_mp$mid_date==1965 ~ 7,
  cdata_mp$mid_date==1975 ~ 8,
  cdata_mp$mid_date==1985 ~ 9,
  cdata_mp$mid_date==1995 ~ 10,
  cdata_mp$mid_date==2005 ~ 11,
  TRUE ~ 12
)
write.csv(cdata_mp, file="Data/moth_coloniser_nspp_MP.csv", row.names=FALSE)

# Load data
cdata <- read.csv("../Moth colonisers/moth_coloniser_dataCERT.csv", header=TRUE)

cdata$mid_date <- factor(cdata$mid_date, levels=c("1905","1915","1925","1935","1945","1955","1965","1975","1985","1995","2005", "2015"))
cdata$coloniser_mode <- factor(cdata$coloniser_mode, levels=c("A", "I"))
cdata$host_plant_status <- factor(cdata$host_plant_status, levels=c("Native", "Non-native"))

# create new datasets 
cdata <- cdata[!is.na(cdata$coloniser_mode),] # 107 species 
cdata1 <- cdata %>% group_by(mid_date, .drop=FALSE) %>% summarise(all_nspp=n())
cdata2 <- cdata %>% filter(coloniser_mode=="I") %>% group_by(mid_date, .drop=FALSE) %>% summarise(imm_nspp=n())
cdata3 <- cdata %>% filter(coloniser_mode=="A") %>% group_by(mid_date, .drop=FALSE) %>% summarise(adv_nspp=n())
cdata4 <- cdata %>% filter(coloniser_mode=="I", host_plant_status=="Native") %>% group_by(mid_date, .drop=FALSE) %>% 
  summarise(imm_nat_nspp=n())
cdata5 <- cdata %>% filter(coloniser_mode=="I", host_plant_status=="Non-native") %>% group_by(mid_date, .drop=FALSE) %>% 
  summarise(imm_nn_nspp=n())
cdata6 <- cdata %>% filter(coloniser_mode=="A", host_plant_status=="Native") %>% group_by(mid_date, .drop=FALSE) %>% 
  summarise(adv_nat_nspp=n())
cdata7 <- cdata %>% filter(coloniser_mode=="A", host_plant_status=="Non-native") %>% group_by(mid_date, .drop=FALSE) %>% 
  summarise(adv_nn_nspp=n())

cdata_cert <- list(cdata1, cdata2, cdata3, cdata4, cdata5, cdata6, cdata7) %>% reduce(full_join, by = "mid_date")
cdata_cert <- as.data.frame(cdata_cert)
cdata_cert$date_num <- case_when(
  cdata_cert$mid_date==1905 ~ 1,
  cdata_cert$mid_date==1915 ~ 2,
  cdata_cert$mid_date==1925 ~ 3,
  cdata_cert$mid_date==1935 ~ 4,
  cdata_cert$mid_date==1945 ~ 5,
  cdata_cert$mid_date==1955 ~ 6,
  cdata_cert$mid_date==1965 ~ 7,
  cdata_cert$mid_date==1975 ~ 8,
  cdata_cert$mid_date==1985 ~ 9,
  cdata_cert$mid_date==1995 ~ 10,
  cdata_cert$mid_date==2005 ~ 11,
  TRUE ~ 12
)
write.csv(cdata_cert, file="Data/moth_coloniser_nspp_CERT.csv", row.names=FALSE)

##### Rates of establishment over time ##### 
rm(list = ls())
cdata <- read.csv("Data/moth_coloniser_nspp_MP.csv", header=TRUE)

cdata <- cdata[!cdata$mid_date=="2015",]

cdata <- cdata %>%
  mutate(nspp_lag = all_nspp - lag(all_nspp))
cdata$rate_decade <- (cdata$nspp_lag/cdata$all_nspp)*100
mean(na.omit(cdata$rate_decade))
sd(na.omit(cdata$rate_decade))

# Run GLMs to see how trend of all species and species groups change over time
# Try this for whole time period and compare to segmented regression - which fit better?
# This will tell us whether the trend differs over time, or stays the same
# And whether the trend is significant

### All species

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

summary(all_spp_seg) # break is at 3.9 (i.e. around 1935)
# Davies' test - tests for the hypothesis that leftSlope=Rightslop (no difference in the two slopes - breakpoint does not exist)
davies.test(all_spp_m, seg.Z= ~date_num) # non-significant - accept null that breakpoint does not exist
slope(all_spp_seg)

# Compare to model with no break - use BIC following Macgregor moth biomass paper 
BIC(all_spp_m) # 67.88
BIC(all_spp_seg) # 71.01
# model with no break fits best
anova(all_spp_m, all_spp_seg, test="LRT")

# Plot trend over time with raw data
pred <- ggpredict(all_spp_m, terms="date_num")
pred$mid_date <- cdata$mid_date

all_spp_p <- ggplot(pred, aes(x=mid_date, y=predicted))+
  geom_line()+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.1)+
  geom_point(data=cdata, aes(x=mid_date, y=all_nspp))+
  scale_x_continuous(breaks=seq(1905,2015,by=20))+
  labs(x="Year", y="Number of species")+
  scale_y_continuous(breaks=seq(0,30,by=5))+
  theme_classic()
all_spp_p
ggsave(all_spp_p, file="Graphs/all_spp_mp_nonseg.png")

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
# ~19% significant increase in the number of species colonising per decade on average


#

### Immigrant species
# Poisson GLM
imm_spp_m <- glm(imm_nspp ~ date_num, data=cdata, family="poisson")
summary(imm_spp_m) # date is significant

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = imm_spp_m)
plot(simulationOutput) # quantile deviations
plotResiduals(simulationOutput, imm_spp$date_num, quantreg = T) 
res = recalculateResiduals(simulationOutput, group = imm_spp$date_num) # recalculate residuals to aggregate residuals per time step
testTemporalAutocorrelation(res, time = unique(imm_spp$date_num)) # non-sig - model is fine?
testDispersion(simulationOutput)
# no assumptions violated

# Segmented GLM
imm_spp_seg <- segmented(imm_spp_m,
                         seg.Z = ~date_num)

summary(imm_spp_seg) # break is at 4 (i.e. 1935)
# Davies' test - tests for the hypothesis that leftSlope=Rightslop (no difference in the two slopes - breakpoint does not exist)
davies.test(imm_spp_m, seg.Z= ~date_num) # non-significant - accept null that breakpoint does not exist

plot(imm_spp_seg)

# Compare to model with no break
BIC(imm_spp_m) # 65.84
BIC(imm_spp_seg) # 67.45
# within 2 points
anova(imm_spp_m, imm_spp_seg, test="LRT") # non-significant
# take the simplest model - non-segmented
slope(imm_spp_seg)


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
ggsave(imm_spp_p, file="Graphs/imm_spp_mp_nonseg.png")

# Also plot segmented result
newdataess <- expand.grid(date_num = c(1,3.217,4,5,6,7,8,9,10,11,12), nspp = 0)
dat <- data.frame(predict(imm_spp_seg, newdata=newdataess, type="response", se.fit=TRUE))
dat$mid_date <- c(1905,1925,1935,1945,1955,1965,1975,1985,1995,2005,2015)

p <- predict(imm_spp_seg, se.fit=TRUE)$se.fit

# calc confidence intervals
setDT(dat)
dat[, rate_lower := fit - 1.96*se.fit]
dat[, rate_upper := fit + 1.96*se.fit]

imm_spp_p2 <- ggplot(dat, aes(x=mid_date, y=fit))+
  geom_line()+
  geom_ribbon(aes(ymin=rate_lower, ymax=rate_upper), alpha=0.1)+
  geom_point(data=imm_spp, aes(x=mid_date, y=nspp))+
  scale_x_continuous(breaks=seq(1905,2015,by=20))+
  labs(x="Year", y="Number of species")+
  theme_classic()
imm_spp_p2 # not convinced this is right but it's close to the plot below
ggsave(imm_spp_p2, file="Graphs/imm_spp_cert_seg.png")

plot(imm_spp_seg, conf.level = 0.95, shade = TRUE) # faster rate of increase between 1905 and 1925/1935 (depending on dataset) (but also more error)
# the y values are not right though - they should be the exponential because it's a poisson model?

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
# ~16% significant increase in the number of species colonising per decade on average

# What is the incidence rate over time? (from segmented model)
output2 <- data.frame(slope(imm_spp_seg))
imm_spp_seg_r <- data.frame(Startdate = 1905, Enddate=2015, # start date of model
                        slope1_rate = exp(output2["slope1","date_num.Est."]), # exponential of the adventive slope
                        slope1_rate_SE = 
                          deltamethod(~exp(x1), 
                                      output2["slope1","date_num.Est."], 
                                      output2["slope1","date_num.St.Err."]^2, ses=TRUE),
                        slope2_rate = exp(output2["slope2","date_num.Est."]), # exponential of the adventive slope
                        slope2_rate_SE = 
                          deltamethod(~exp(x1), 
                                      output2["slope2","date_num.Est."], 
                                      output2["slope2","date_num.St.Err."]^2, ses=TRUE) ) # standard error around exponential of immigrant slope

# Produce confidence intervals for the rate estimates
setDT(imm_spp_seg_r)
imm_spp_seg_r[, slope1_rate_lower := slope1_rate - 1.96*slope1_rate_SE]
imm_spp_seg_r[, slope1_rate_upper := slope1_rate + 1.96*slope1_rate_SE]
imm_spp_seg_r[, slope2_rate_lower := slope2_rate - 1.96*slope2_rate_SE]
imm_spp_seg_r[, slope2_rate_upper := slope2_rate + 1.96*slope2_rate_SE]
imm_spp_seg_r
# not right

#

### Adventive species
# Poisson GLM
adv_spp_m <- glm(adv_nspp ~ date_num, data=cdata, family="poisson")
summary(adv_spp_m) # date is significant

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = adv_spp_m)
plot(simulationOutput) # quantile deviations
plotResiduals(simulationOutput, adv_spp$date_num, quantreg = T) 
res = recalculateResiduals(simulationOutput, group = adv_spp$date_num) # recalculate residuals to aggregate residuals per time step
testTemporalAutocorrelation(res, time = unique(adv_spp$date_num)) # non-sig - model is fine?
testDispersion(simulationOutput)
# no assumptions violated

# Segmented GLM
adv_spp_seg <- segmented(adv_spp_m,
                         seg.Z = ~date_num)

summary(adv_spp_seg) # break is at 8.6 (i.e. ~1975)
# Davies' test - tests for the hypothesis that leftSlope=Rightslop (no difference in the two slopes - breakpoint does not exist)
davies.test(adv_spp_m, seg.Z= ~date_num) # non-significant - accept null that breakpoint does not exist

# Compare to model with no break
BIC(adv_spp_m) # 51.15
BIC(adv_spp_seg) # 55.35
# non-segmented model has lower BIC

# NOTE: V. SIMILAR MODELS FOR TD CATEGORIES - SEGMENTED MODEL SLIGHTLY LOWER BIC (53.75 VS 53.22)
anova(adv_spp_m, adv_spp_seg, test="LRT") # just non-sig (0.06)
slope(adv_spp_seg)

plot(adv_spp_seg, conf.level = 0.95, shade=TRUE) # decline between 1905 and 1920s (large error - not many points), then increase
# still think non-segmented is better..

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
ggsave(adv_spp_p, file="Graphs/adv_spp_mp_nonseg.png")

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
# ~25% significant increase in the number of species colonising per decade on average

## Try and plot segmented model results (for certain + 2015 removed) with ggplot

plot(cdata$date_num, cdata$adv_nspp, pch=19, xaxt="n", xlab="", ylab="", bty = "n", xaxt = "n", yaxt = "n", ylim=c(-1,20))
box("plot", bty = "l", lwd = 1.5)
axis(1, at=1:11, labels=c(1905,1915,1925,1935,1945,1955,1965,1975,1985,1995,2005), lwd=0, lwd.ticks = 1.5, tck=-0.01)
axis(side = 2, lwd = 0, lwd.ticks = 1.5, las = 1, tck=-0.01)
plot(adv_spp_seg, conf.level = 0.95, shade=TRUE, link=FALSE, add=TRUE, lwd=1)
title(ylab="Number of species", line=2.3, cex.lab=1.2, family="Calibri Light")
title(xlab="Year", line=2.3, cex.lab=1.2, family="Calibri Light")


newdataess <- expand.grid(date_num = c(1,2,3,4,5,6,7,8.882,9,11), nspp = 0)
dat <- data.frame(predict(adv_spp_seg, newdata=newdataess, type="response", se.fit=TRUE))
dat$mid_date <- c(1905,1915,1925,1935,1945,1955,1965,1975,1985,2005)
dat$sig <- c("Non-significant", "Non-significant", "Non-significant","Non-significant","Non-significant","Non-significant","Non-significant","Non-significant",
             "Significant","Significant")
p <- predict(adv_spp_seg, se.fit=TRUE)$se.fit

# calc confidence intervals
newdataesw <- expand.grid(date_num = cdata$date_num, nspp = 0)
dat <- predict(adv_spp_seg, newdata = newdataesw, type="response", se.fit=TRUE)
dat <- as.data.frame(dat)
dat$mid_date <- c(1905,1915,1925,1935,1945,1955,1965,1975,1985,1995,2005)
setDT(dat)
dat[, rate_lower := fit - 1.96*se.fit]
dat[, rate_upper := fit + 1.96*se.fit]

adv_spp_p2 <- ggplot(dat, aes(x=mid_date, y=fit))+
  geom_line()+
  geom_ribbon(aes(ymin=rate_lower, ymax=rate_upper), alpha=0.1)+
  geom_point(data=cdata, aes(x=mid_date, y=adv_nspp))+
  scale_x_continuous(breaks=seq(1905,2015,by=20))+
  labs(x="Year", y="Number of species")+
  theme_classic()
adv_spp_p2 # not convinced this is right but it's close to the plot below
ggsave(adv_spp_p2, file="Graphs/adv_spp_cert_seg.png")


dat <- predict(adv_spp_seg, newdata = newdataesw, type="response", se.fit=TRUE)
dat <- as.data.frame(dat)
dat$mid_date <- c(1905,1915,1925,1935,1945,1955,1965,1975,1985,1995,2005)




## SO FAR: NO CHANGE IN THE RATE OF ESTABLISHMENT OVER TIME - RATE IS CONSISTENT BETWEEN 1905 AND 2015 FOR IMMIGRANTS AND ADVENTIVES

# NOTE: SUBSET AND TD CATEGORIES - IMMIGRANTS MIGHT BE BETTER WITH SEGMENTED MODEL (FASTER RATE FIRST, THEN SLOWER RATE OF INCREASE)
# ADVENTIVES ALSO POSSIBLE BETTER FIT FOR SEGMENTED (BUT LESS CONVINCED)

### Immigrant + native host plants
# Poisson GLM
imm_nat_spp_m <- glm(imm_nat_nspp ~ date_num, data=cdata, family="poisson")
summary(imm_nat_spp_m) # date is significant

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = imm_nat_spp_m)
plot(simulationOutput)
testDispersion(simulationOutput)
# no assumptions violated

# Segmented GLM
imm_nat_spp_seg <- segmented(imm_nat_spp_m,
                         seg.Z = ~date_num)

summary(imm_nat_spp_seg) # break is at 4 (i.e. ~1935)
# Davies' test - tests for the hypothesis that leftSlope=Rightslop (no difference in the two slopes - breakpoint does not exist)
davies.test(imm_nat_spp_m, seg.Z= ~date_num) # non-significant - accept null that breakpoint does not exist

# Compare to model with no break
BIC(imm_nat_spp_m) # 59.65
BIC(imm_nat_spp_seg) # 63.97
# non-segmented model has lower BIC
anova(imm_nat_spp_m, imm_nat_spp_seg, test="LRT")
slope(imm_nat_spp_seg)

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
ggsave(imm_nat_spp_p, file="Graphs/imm_nat_spp_mp_nonseg.png")

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

### Immigrant + non-native host plants
# Poisson GLM
imm_nn_spp_m <- glm(imm_nn_nspp ~ date_num, data=cdata, family="poisson")
summary(imm_nn_spp_m) # date is non-significant

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = imm_nn_spp_m)
plot(simulationOutput) # quantile deviations
res = recalculateResiduals(simulationOutput, group = imm_nn_spp$date_num) # recalculate residuals to aggregate residuals per time step
testTemporalAutocorrelation(res, time = unique(imm_nn_spp$date_num)) # non-sig - model is fine?
testDispersion(simulationOutput)

# Segmented GLM
imm_nn_spp_seg <- segmented(imm_nn_spp_m,
                             seg.Z = ~date_num)

summary(imm_nn_spp_seg) # break is at 4.8 (i.e. ~1945)
slope(imm_nn_spp_seg) # both slopes are non-significant - no significant trend over time
davies.test(imm_nn_spp_m,
            seg.Z = ~date_num)

BIC(imm_nn_spp_m)
BIC(imm_nn_spp_seg) # segmented is better fit but no trend in time
anova(imm_nn_spp_m, imm_nn_spp_seg, test="LRT")
slope(imm_nn_spp_seg)

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
ggsave(imm_nn_spp_p, file="Graphs/imm_nn_spp_mp_nonseg.png")

## SUMMARY: NATIVE HOST PLANT IMMIGRANT SPECIES SIGNIFICANTLY INCREASE OVER TIME AT THE SAME RATE, 
## NON-NATIVE HOST PLANT IMMIGRANT SPECIES SHOW NO SIGNIFICANT CHANGE OVER TIME

# 2015 removed date is significant
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
#

### Adventives + native host plants
# Poisson GLM
adv_nat_spp_m <- glm(adv_nat_nspp ~ date_num, data=cdata, family="poisson")
summary(adv_nat_spp_m) # date is significant

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = adv_nat_spp_m)
plot(simulationOutput)
testDispersion(simulationOutput)
# no assumptions violated

# Segmented GLM
adv_nat_spp_seg <- segmented(adv_nat_spp_m,
                             seg.Z = ~date_num)

summary(adv_nat_spp_seg) # break is at 8 (i.e. ~1975)
# Davies' test - tests for the hypothesis that leftSlope=Rightslop (no difference in the two slopes - breakpoint does not exist)
davies.test(adv_nat_spp_m, seg.Z= ~date_num) # non-significant - accept null that breakpoint does not exist
#NOTE: SIGNIFICANT FOR TD CATEGORIES
slope(adv_nat_spp_seg)

# Compare to model with no break
BIC(adv_nat_spp_m) # 43.02
BIC(adv_nat_spp_seg) # 46.14
# non-segmented model has lower BIC

# NOTE: V. SIMILAR BIC (40.55 VS 40.99) WITH CERTAIN SPECIES
# SEGMENTED HAS LOWER BIC FOR TD COMMENTS

anova(adv_nat_spp_m, adv_nat_spp_seg, test="LRT") # non-significant with certain species, significant with TD

plot(adv_nat_spp_seg, conf.level=0.95, shade=TRUE) # DECLINE BETWEEN 1905 AND ~1965, THEN INCREASE TO 2015
# POSSIBLY DRIVEN BY STRING OF ZEROS BETWEEN 1955 AND 1985
slope(adv_nat_spp_seg)

# Also plot segmented result
mylocations <- c("topleft", "topright")
par(mfrow=c(1,2))

plot(cdata$date_num, cdata$adv_nspp, pch=19, xaxt="n", xlab="", ylab="", bty = "n", xaxt = "n", yaxt = "n", ylim=c(-1,20))
box("plot", bty = "l", lwd = 1.5)
axis(1, at=1:11, labels=c(1905,1915,1925,1935,1945,1955,1965,1975,1985,1995,2005), lwd=0, lwd.ticks = 1.5, tck=-0.01)
axis(side = 2, lwd = 0, lwd.ticks = 1.5, las = 1, tck=-0.01)
plot(adv_spp_seg, conf.level = 0.95, shade=TRUE, link=FALSE, add=TRUE, lwd=1)
title(ylab="Number of species", line=2.3, cex.lab=1.2, family="Calibri Light")
title(xlab="Year", line=2.3, cex.lab=1.2, family="Calibri Light")


plot(cdata$date_num, cdata$adv_nat_nspp, pch=19, xaxt="n", xlab="", ylab="", bty = "n", xaxt = "n", yaxt = "n", ylim=c(-1,12))
box("plot", bty = "l", lwd = 1.5)
axis(1, at=1:11, labels=c(1905,1915,1925,1935,1945,1955,1965,1975,1985,1995,2005), lwd=0, lwd.ticks = 1.5, tck=-0.01)
axis(side = 2, lwd = 0, lwd.ticks = 1.5, las = 1, tck=-0.01)
plot(adv_nat_spp_seg, conf.level = 0.95, shade=TRUE, link=FALSE, add=TRUE, lwd=1)
title(ylab="Number of species", line=2.3, cex.lab=1.2)
title(xlab="Year", line=2.3, cex.lab=1.2)

newdataesw <- expand.grid(date_num = cdata$date_num, nspp = 0)
dat2 <- predict(adv_nat_spp_seg, newdata = newdataesw, type="response", se.fit=TRUE)
dat2 <- as.data.frame(dat2)
dat2$mid_date <- c(1905,1915,1925,1935,1945,1955,1965,1975,1985,1995,2005)
setDT(dat2)
dat2[, rate_lower := fit - 1.96*se.fit]
dat2[, rate_upper := fit + 1.96*se.fit]

adv_nat_spp_p2 <- ggplot(dat2, aes(x=mid_date, y=fit))+
  geom_line()+
  geom_ribbon(aes(ymin=rate_lower, ymax=rate_upper), alpha=0.1)+
  geom_point(data=cdata, aes(x=mid_date, y=adv_nat_nspp))+
  scale_x_continuous(breaks=seq(1905,2015,by=20))+
  labs(x="Year", y="Number of species")+
  theme_classic()
adv_nat_spp_p2 # not convinced this is right but it's close to the plot below
ggsave(adv_nat_spp_p2, file="Graphs/adv_nat_spp_cert_seg.png")

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
ggsave(adv_nat_spp_p, file="Graphs/adv_nat_spp_mp_nonseg.png")

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

### Adventive + non-native host plants
# Poisson GLM
adv_nn_spp_m <- glm(adv_nn_nspp ~ date_num, family="poisson", data=cdata)
summary(adv_nn_spp_m) # date is significant

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = adv_nn_spp_m)
plot(simulationOutput)
testDispersion(simulationOutput)
# possible underdispersion (cert + 2015 removed)

# Segmented GLM
adv_nn_spp_seg <- segmented(adv_nn_spp_m,
                             seg.Z = ~date_num)

summary(adv_nn_spp_seg) # break is at 2 (i.e. ~1915)
# Davies' test - tests for the hypothesis that leftSlope=Rightslop (no difference in the two slopes - breakpoint does not exist)
davies.test(adv_nn_spp_m, seg.Z= ~date_num) # non-significant - accept null that breakpoint does not exist
slope(adv_nn_spp_seg)
plot(adv_nn_spp_seg)

df <- cdata[,c(7,9)]
dput(df)
df2 <- structure(list(nspp = c(1L, 1L, 0L, 1L, 2L, 0L, 2L, 0L, 
                                      1L, 3L, 8L, 3L), date_num = 1:12), class = "data.frame", row.names = c(NA, 
                                                                                                             -12L))
df <- structure(list(adv_nn_nspp = c(0L, 0L, 2L, 1L, 0L, 2L, 3L, 2L, 
                                     1L, 2L, 6L, 3L), date_num = 1:12), class = "data.frame", row.names = c(NA, 
                                                                                                            -12L))

# Compare to model with no break
BIC(adv_nn_spp_m) # 39.39
BIC(adv_nn_spp_seg) # 41.62
# non-segmented model has lower BIC

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
ggsave(adv_nn_spp_p, file="Graphs/adv_nn_spp_mp_nonseg.png")

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
# ~22% significant increase in the number of species colonising per decade on average

## SUMMARY: BOTH ADVENTIVE GROUPS RATE INCREASE DOESN'T CHANGE OVER TIME, BOTH INCREASING ~20% PER DECADE ON AVERAGE
## ALL RESULTS V. SIMILAR (APART FROM IMMIGRANT NATIVE SEGMENTED MODEL BETTER) WITH CERTAIN SPECIES (N=104)

## TD CATEGORIES DIFFERENCE: IMMIGRANTS BETTER WITH SEGMENTED MODEL, ADVENTIVES TOO (BUT LESS CERTAIN), PLUS ADVENTIVES ON NATIVE
## HOST-PLANTS ALSO BETTER WITH SEGMENTED MODEL

###

## put all plots together?

# 1. All species
# 2. Immigrant species
# 3. Adventive species
# 4. Immigrant + native host plants
# 5. Immigrant + non-native host plants
# 6. Adventive + native host plants
# 7. Adventive + non-native host plants

# Plot total, immigrants and adventives on one plot
pred7 <- ggpredict(all_spp_m, terms="date_num")
pred8 <- ggpredict(imm_spp_m, terms="date_num")
pred9 <- ggpredict(adv_spp_m, terms="date_num")
pred7$group <- "all_spp"
pred8$group <- "immigrants"
pred9$group <- "adventives"

pred_final3 <- rbind(pred7, pred8, pred9)
pred_final3$group <- factor(pred_final3$group, levels=c("all_spp", "immigrants", "adventives"))

all_final <- ggplot(pred_final3, aes(x=x, y=predicted, fill=group))+
  geom_line(data=pred_final3, aes(colour=group))+
  geom_ribbon(data=pred_final3, aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.1)+
  #geom_point(data=raw_dat3, aes(x=date_num, y=nspp, colour=coloniser_mode))+
  #scale_x_continuous(breaks=seq(1905,2015,by=20))+
  labs(x="Year", y="Number of species")+
  theme_classic()
all_final

ggsave(all_final, file="Graphs/all_imm_adv_nonseg.png")

# Try plotting 3 immigrant lines on one plot
pred1 <- ggpredict(imm_spp_m, terms="date_num")
pred2 <- ggpredict(imm_nat_spp_m, terms="date_num")
pred3 <- ggpredict(imm_nn_spp_m, terms="date_num")
# pred3$conf.low <- NA
# pred3$conf.high <- NA

pred1$group <- "immigrant"
pred1$significance <- "significant"
pred1$mid_date <- cdata$mid_date
pred2$group <- "immigrant_native"
pred2$significance <- "significant"
pred2$mid_date <- cdata$mid_date
pred3$group <- "immigrant_nonnative"
pred3$significance <- "non-significant"
pred3$mid_date <- cdata$mid_date

imm <- cdata[,c(1,3,5:6)]
library(data.table)
imm2 <- melt(setDT(imm), id.vars = "mid_date", variable.name = "group")
imm2$group <- as.character(imm2$group)
imm2$group[imm2$group == "imm_nspp"] <- "immigrant"
imm2$group[imm2$group == "imm_nat_nspp"] <- "immigrant_native"
imm2$group[imm2$group == "imm_nn_nspp"] <- "immigrant_nonnative"


pred_final <- rbind(pred1, pred2, pred3)
pred_final <- merge(pred_final, imm2, by=c("group", "mid_date"))
legend_title <- ""
pred_final$significance <- factor(pred_final$significance, levels=c("significant", "non-significant"))
imm_final <- ggplot(data=pred_final)+
  geom_line(aes(x=mid_date, y=predicted, colour=group, linetype=significance), lwd=2)+
  geom_ribbon(aes(x=mid_date, y=predicted, ymin=conf.low, ymax=conf.high, fill=group), alpha=0.2)+
  geom_point(aes(x=mid_date, y=value, colour=group, shape=group), size=3)+
  scale_fill_manual(legend_title,labels=c("Immigrants","Immigrants on native \nhost plants", "Immigrants on non-native \nhost plants"),
                    values = c("#1B9E77", "#D95F02", "#7570B3"))+
  scale_colour_manual(legend_title,labels=c("Immigrants","Immigrants on native \nhost plants", "Immigrants on non-native \nhost plants"),
                    values = c("#1B9E77", "#D95F02", "#7570B3"))+
  scale_shape_manual(legend_title,labels=c("Immigrants","Immigrants on native \nhost plants", "Immigrants on non-native \nhost plants"),
                     values = c(19,17,15))+
  guides(linetype = "none", fill = guide_legend(byrow = TRUE))+
  scale_x_continuous(breaks=seq(1905,2015,by=20))+
  scale_y_continuous(breaks=seq(0,30,by=5), limits=c(0,25))+
  labs(x="Year", y="Number of species")+
  theme_classic()+
  theme(legend.position = "bottom", text=element_text(size=24))
imm_final
ggsave(imm_final, file="Graphs/all_imm_nonseg.png")


# Try plotting 3 adventive lines on one plot
pred4 <- ggpredict(adv_spp_m, terms="date_num")
pred5 <- ggpredict(adv_nat_spp_m, terms="date_num")
pred6 <- ggpredict(adv_nn_spp_m, terms="date_num")
pred4$group <- "adventive"
pred4$mid_date <- cdata$mid_date
pred5$group <- "adventive_native"
pred5$mid_date <- cdata$mid_date
pred6$group <- "adventive_nonnative"
pred6$mid_date <- cdata$mid_date

pred_final2 <- rbind(pred4, pred5, pred6)
adv <- cdata[,c(1,4,7:8)]
library(data.table)
adv2 <- melt(setDT(adv), id.vars = "mid_date", variable.name = "group")
adv2$group <- as.character(adv2$group)
adv2$group[adv2$group == "adv_nspp"] <- "adventive"
adv2$group[adv2$group == "adv_nat_nspp"] <- "adventive_native"
adv2$group[adv2$group == "adv_nn_nspp"] <- "adventive_nonnative"
pred_final2 <- merge(pred_final2, adv2, by=c("group", "mid_date"))

adv_final <- ggplot(pred_final2)+
  geom_line(aes(x=mid_date, y=predicted, colour=group), lwd=2)+
  geom_ribbon(aes(x=mid_date, y=predicted, ymin=conf.low, ymax=conf.high, fill=group), alpha=0.2)+
  geom_point(aes(x=mid_date, y=value, colour=group, shape=group), size=3)+
  scale_fill_manual(legend_title,labels=c("Adventives","Adventives on native \nhost plants", "Adventives on non-native \nhost plants"),
                    values = c("#1B9E77", "#D95F02", "#7570B3"))+
  scale_colour_manual(legend_title,labels=c("Adventives","Adventives on native \nhost plants", "Adventives on non-native \nhost plants"),
                      values = c("#1B9E77", "#D95F02", "#7570B3"))+
  scale_shape_manual(legend_title,labels=c("Adventives","Adventives on native \nhost plants", "Adventives on non-native \nhost plants"),
                     values = c(19,17,15))+
  scale_x_continuous(breaks=seq(1905,2015,by=20))+
  scale_y_continuous(breaks=seq(0,30,by=5), limits=c(0,25))+
  labs(x="Year", y="Number of species")+
  guides(fill = guide_legend(byrow = TRUE))+
  theme_classic()+
  theme(legend.position = "bottom", text=element_text(size=24))
adv_final
ggsave(adv_final, file="Graphs/all_adv_nonseg.png")

## adventive plots with segmented regressions
dat$group <- "adventive"
dat2$group <- "adventive_native"
colnames(dat) <- c("predicted", "std.error","residual.scale","mid_date","conf.low", "conf.high","group")
colnames(dat2) <- c("predicted", "std.error","residual.scale","mid_date","conf.low", "conf.high","group")
dat$residual.scale <- NULL
dat2$residual.scale <- NULL

pred6 <- ggpredict(adv_nn_spp_m, terms="date_num")
pred6$group <- "adventive_nonnative"
pred6$mid_date <- cdata$mid_date
pred6$x <- NULL

pred_final2 <- rbind(dat, dat2, pred6)

adv_final <- ggplot(pred_final2, aes(x=mid_date, y=predicted, fill=group))+
  geom_line(aes(colour=group), lwd=2)+
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.2)+
  #geom_point(data=adv_nn_spp, aes(x=mid_date, y=nspp))+
  scale_fill_manual(legend_title,labels=c("Adventives","Adventives on native \nhost plants", "Adventives on non-native \nhost plants"),
                    values = c("#1B9E77", "#D95F02", "#7570B3"))+
  scale_colour_manual(legend_title,labels=c("Adventives","Adventives on native \nhost plants", "Adventives on non-native \nhost plants"),
                      values = c("#1B9E77", "#D95F02", "#7570B3"))+
  scale_x_continuous(breaks=seq(1905,2005,by=20))+
  scale_y_continuous(breaks=seq(0,20,by=5), limits=c(-1,20))+
  labs(x="Year", y="Number of species")+
  guides(fill = guide_legend(byrow = TRUE))+
  theme_classic()+
  theme(legend.position = "bottom", text=element_text(size=24))
adv_final
ggsave(adv_final, file="Graphs/all_adv_seg.png", width=10, height=10, units="in")


# put all plots together?
library(ggpubr)
final <- ggarrange(imm_final, adv_final, 
          labels = c("(a)", "(b)"),
          ncol = 2, nrow = 1, font.label=list(color="black",size=24))
ggsave("Graphs/final_imm_adv.png", final, width = 20, height = 10, units = "in")

# put adventives and immigrants on one plot for supplementary material
pred_final3 <- pred_final[pred_final$group=="immigrant",]
pred_final3$significance <- NULL
adv <- pred_final2[pred_final2$group=="adventive",]
pred_final3 <- rbind(pred_final3, adv)

imm_adv <- ggplot(data=pred_final3)+
  geom_line(aes(x=mid_date, y=predicted, colour=group), lwd=1)+
  geom_ribbon(aes(x=mid_date, y=predicted, ymin=conf.low, ymax=conf.high, fill=group), alpha=0.2)+
  geom_point(aes(x=mid_date, y=value, colour=group, shape=group), size=2)+
  scale_fill_manual(legend_title,labels=c("Adventives","Immigrants"),
                    values = c("#1B9E77", "#7570B3"))+
  scale_colour_manual(legend_title,labels=c("Adventives","Immigrants"),
                      values = c("#1B9E77", "#7570B3"))+
  scale_shape_manual(legend_title,labels=c("Adventives","Immigrants"),
                     values = c(19,15))+
  guides(linetype = "none", fill = guide_legend(byrow = TRUE))+
  scale_x_continuous(breaks=seq(1905,2015,by=20))+
  scale_y_continuous(breaks=seq(0,30,by=5), limits=c(0,25))+
  labs(x="Year", y="Number of species")+
  theme_classic()+
  theme(legend.position = "bottom", text=element_text(size=18))
imm_adv
ggsave(imm_adv, file="Graphs/imm_adv_nonseg.png")



########################################################################################
########################################################################################
########################################################################################
########################################################################################

## Group/time interaction
cdata <- read.csv("../Moth colonisers/moth_coloniser_dataMP.csv", header=TRUE)

cdata$mid_date <- factor(cdata$mid_date, levels=c("1905","1915","1925","1935","1945","1955","1965","1975","1985","1995","2005", "2015"))
cdata$coloniser_mode <- factor(cdata$coloniser_mode, levels=c("A", "I"))
cdata$host_plant_status <- factor(cdata$host_plant_status, levels=c("Native", "Non-native"))

# create new datasets 
cdata <- cdata[!is.na(cdata$coloniser_mode),] # 144 species 
# need to remove two species now considered overlooked residents (Carpatolechia alburnella and Dendrolimus pini)
cdata <- cdata[!cdata$scientific_name=="Carpatolechia alburnella",]
cdata <- cdata[!cdata$scientific_name=="Dendrolimus pini",]
# 142 species now 
# cdata <- cdata[!cdata$mid_date=="2015",]

# Immigrants vs Adventives
cdata$coloniser_mode <- factor(cdata$coloniser_mode, levels=c("A", "I"))
imm_ad_spp <- cdata %>% group_by(mid_date, coloniser_mode, .drop=FALSE) %>% summarise(nspp=n())
imm_ad_spp <- as.data.frame(imm_ad_spp)
imm_ad_spp$date_num <- case_when(
  imm_ad_spp$mid_date==1905 ~ 1,
  imm_ad_spp$mid_date==1915 ~ 2,
  imm_ad_spp$mid_date==1925 ~ 3,
  imm_ad_spp$mid_date==1935 ~ 4,
  imm_ad_spp$mid_date==1945 ~ 5,
  imm_ad_spp$mid_date==1955 ~ 6,
  imm_ad_spp$mid_date==1965 ~ 7,
  imm_ad_spp$mid_date==1975 ~ 8,
  imm_ad_spp$mid_date==1985 ~ 9,
  imm_ad_spp$mid_date==1995 ~ 10,
  imm_ad_spp$mid_date==1995 ~ 12,
  TRUE ~ 12
)
mod <- glm(nspp ~ date_num*coloniser_mode, data=imm_ad_spp, family="poisson")
summary(mod) # non-significant for all 4 datasets

# Immigrant native vs non-native
cdata2 <- cdata[!is.na(cdata$host_plant_status),] # 107 species 
cdata_test2$mid_date <- factor(cdata_test2$mid_date, levels=c("1905","1915","1925","1935","1945","1955","1965","1975","1985","1995","2005"))
cdata2$mid_date <- factor(cdata2$mid_date, levels=c("1905","1915","1925","1935","1945","1955","1965","1975","1985","1995","2005","2015"))
cdata2$host_plant_status <- factor(cdata2$host_plant_status, levels=c("Native", "Non-native"))
cdata_test2$host_plant_status <- factor(cdata_test2$host_plant_status, levels=c("Native", "Non-native"))

imm_natnn_spp <- cdata2 %>% filter(coloniser_mode=="I") %>% group_by(mid_date, host_plant_status, .drop=FALSE) %>% summarise(nspp=n())
imm_natnn_spp$mid_date <- as.numeric(levels(imm_natnn_spp$mid_date))[imm_natnn_spp$mid_date]
imm_natnn_spp <- as.data.frame(imm_natnn_spp)
imm_natnn_spp$date_num <- case_when(
  imm_natnn_spp$mid_date==1905 ~ 1,
  imm_natnn_spp$mid_date==1915 ~ 2,
  imm_natnn_spp$mid_date==1925 ~ 3,
  imm_natnn_spp$mid_date==1935 ~ 4,
  imm_natnn_spp$mid_date==1945 ~ 5,
  imm_natnn_spp$mid_date==1955 ~ 6,
  imm_natnn_spp$mid_date==1965 ~ 7,
  imm_natnn_spp$mid_date==1975 ~ 8,
  imm_natnn_spp$mid_date==1985 ~ 9,
  imm_natnn_spp$mid_date==1995 ~ 10,
  TRUE ~ 11
)

mod2 <- glm(nspp ~ date_num*host_plant_status, data=imm_natnn_spp, family="poisson")
summary(mod2) # non-significant for all 4 datasets

# Adventive native vs non-native
ad_natnn_spp <- cdata2 %>% filter(coloniser_mode=="A") %>% group_by(mid_date, host_plant_status, .drop=FALSE) %>% summarise(nspp=n())
ad_natnn_spp$mid_date <- as.numeric(levels(ad_natnn_spp$mid_date))[ad_natnn_spp$mid_date]
ad_natnn_spp <- as.data.frame(ad_natnn_spp)
ad_natnn_spp$date_num <- case_when(
  ad_natnn_spp$mid_date==1905 ~ 1,
  ad_natnn_spp$mid_date==1915 ~ 2,
  ad_natnn_spp$mid_date==1925 ~ 3,
  ad_natnn_spp$mid_date==1935 ~ 4,
  ad_natnn_spp$mid_date==1945 ~ 5,
  ad_natnn_spp$mid_date==1955 ~ 6,
  ad_natnn_spp$mid_date==1965 ~ 7,
  ad_natnn_spp$mid_date==1975 ~ 8,
  ad_natnn_spp$mid_date==1985 ~ 9,
  ad_natnn_spp$mid_date==1995 ~ 10,
  TRUE ~ 11
)
mod3 <- glm(nspp ~ date_num*host_plant_status, data=ad_natnn_spp, family="poisson")
summary(mod3) # non-significant for all 4 datasets



