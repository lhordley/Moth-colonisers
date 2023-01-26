##########################
#### user: Lisbeth Hordley
#### date: June 2022
#### info: Rates of establishment over time (main dataset)

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
cdata <- read.csv("Coloniser data/moth_coloniser_nspp.csv", header=TRUE)

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
ggsave(all_spp_p, file="Graphs/Figure1.png")

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

### 1b. Immigrant species
# Poisson GLM
imm_spp_m <- glm(imm_nspp ~ date_num, data=cdata, family="poisson")
summary(imm_spp_m) 

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = imm_spp_m)
plot(simulationOutput) # quantile deviations
plotResiduals(simulationOutput, cdata$date_num, quantreg = T) 
res = recalculateResiduals(simulationOutput, group = cdata$date_num) # recalculate residuals to aggregate residuals per time step
testTemporalAutocorrelation(res, time = unique(cdata$date_num)) # non-sig
testDispersion(simulationOutput)
# no assumptions violated

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
# within 2 points
# take the simplest model - non-segmented

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
# ~16% significant increase in the number of species colonising per decade on average

#

### 1c. Adventive species
# Poisson GLM
adv_spp_m <- glm(adv_nspp ~ date_num, data=cdata, family="poisson")
summary(adv_spp_m) # date is significant

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = adv_spp_m)
plot(simulationOutput) # quantile deviations
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
# ~25% significant increase in the number of species colonising per decade on average

#

### 1d. Immigrant + native host plants
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
slope(imm_nat_spp_seg)

# Compare to model with no break
BIC(imm_nat_spp_m) 
BIC(imm_nat_spp_seg) 
# non-segmented model has lower BIC

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
# non-seg model better fit

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
# non-segmented model has lower BIC

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
# ~22% significant increase in the number of species colonising per decade on average


####################################
############## GRAPHS ############## 
####################################

## Plot all 3 immigrant trends (all immigrants, immigrants + native and immigrants + non-native) 

pred1 <- ggpredict(imm_spp_m, terms="date_num")
pred2 <- ggpredict(imm_nat_spp_m, terms="date_num")
pred3 <- ggpredict(imm_nn_spp_m, terms="date_num")

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


# Plot all 3 adventive trends 
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

# Put all plots together
library(ggpubr)
final <- ggarrange(imm_final, adv_final, 
                   labels = c("(a)", "(b)"),
                   ncol = 2, nrow = 1, font.label=list(color="black",size=24))
ggsave("Graphs/Figure2.png", final, width = 20, height = 10, units = "in")

# put adventives and immigrants on one plot for supplementary material
imm <- pred_final[pred_final$group=="immigrant",]
imm$significance <- NULL
adv <- pred_final2[pred_final2$group=="adventive",]
pred_final3 <- rbind(imm, adv)

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
ggsave(imm_adv, file="Graphs/FigureS2.png")


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
plot(simulationOutput) # quantile deviations
res = recalculateResiduals(simulationOutput, group = cdata$date_num) # recalculate residuals to aggregate residuals per time step
testTemporalAutocorrelation(res, time = unique(cdata$date_num)) # non-sig
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
anova(all_spp_m, all_spp_seg, test="LRT")

# What is the incidence rate over time?
output <- coef(summary(all_spp_m))
all_spp_r <- data.frame(Startdate = 1905, Enddate=2005, # start date of model
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
summary(imm_spp_m) 

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = imm_spp_m)
plot(simulationOutput) 
testDispersion(simulationOutput)

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
# within 2 points - take simplest model (non-segmented)

# What is the incidence rate over time? (from non-segmented model)
output <- coef(summary(imm_spp_m))
imm_spp_r <- data.frame(Startdate = 1905, Enddate=2005, # start date of model
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
# within 2 points so take simplest model (non-segmented)

# What is the incidence rate over time?
output <- coef(summary(adv_spp_m))
adv_spp_r <- data.frame(Startdate = 1905, Enddate=2005, # start date of model
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

### 2d. Immigrant + native host plants
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

summary(imm_nat_spp_seg) 
# Davies' test - tests for the hypothesis that leftSlope=Rightslop (no difference in the two slopes - breakpoint does not exist)
davies.test(imm_nat_spp_m, seg.Z= ~date_num) # non-significant - accept null that breakpoint does not exist
slope(imm_nat_spp_seg)

# Compare to model with no break
BIC(imm_nat_spp_m) 
BIC(imm_nat_spp_seg) 
# non-segmented model has lower BIC

# What is the incidence rate over time?
output <- coef(summary(imm_nat_spp_m))
imm_nat_spp_r <- data.frame(Startdate = 1905, Enddate=2005, # start date of model
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
res = recalculateResiduals(simulationOutput, group = cdata$date_num) # recalculate residuals to aggregate residuals per time step
testTemporalAutocorrelation(res, time = unique(cdata$date_num)) # non-sig
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
# within 2 points so take simplest model (non-segmented)

# What is the incidence rate over time?
output <- coef(summary(imm_nn_spp_m))
imm_nn_spp_r <- data.frame(Startdate = 1905, Enddate=2005, # start date of model
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
# within 2 points so take simplest model (non-segmented)

# What is the incidence rate over time?
output <- coef(summary(adv_nat_spp_m))
adv_nat_spp_r <- data.frame(Startdate = 1905, Enddate=2005, # start date of model
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

#

### 2g. Adventive + non-native host plants
# Poisson GLM
adv_nn_spp_m <- glm(adv_nn_nspp ~ date_num, family="poisson", data=cdata)
summary(adv_nn_spp_m) 

# Check model fit
simulationOutput <- simulateResiduals(fittedModel = adv_nn_spp_m)
plot(simulationOutput)
testDispersion(simulationOutput)

# Segmented GLM
adv_nn_spp_seg <- segmented(adv_nn_spp_m,
                            seg.Z = ~date_num)

summary(adv_nn_spp_seg) # coefficients are far too high - segmented model is struggling with too many zeros
# just stick with non-segmented model
BIC(adv_nn_spp_m) 

# What is the incidence rate over time?
output <- coef(summary(adv_nn_spp_m))
adv_nn_spp_r <- data.frame(Startdate = 1905, Enddate=2005, # start date of model
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

