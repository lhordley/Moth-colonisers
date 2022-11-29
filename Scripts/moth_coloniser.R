##########################
#### user: Lisbeth Hordley
#### date: May 2022
#### info: Analysing changing drivers of moth colonisation
options(scipen=999)

# Packages
library(data.table)
library(ggplot2)
library(DHARMa)
library(dplyr)
library(AER)
library(msm)
library(ggeffects)

# Read in data
cdata <- read.csv("../Moth colonisers/moth_coloniser_data2.csv", header=TRUE)

# Natural colonisers
nat_col <- cdata[!cdata$coloniser_mode=="A",]
str(nat_col)
nat_col$host_plant_status <- factor(nat_col$host_plant_status, levels=c("Native", "Non-native"))

nat_col2 <- nat_col %>% group_by(mid_date, host_plant_status, .drop=FALSE) %>% summarise(nspp=n())

# Plot raw data and smoothed curves
ggplot(nat_col2, aes(mid_date, nspp, col=host_plant_status))+
  geom_line()+ geom_point()+
  geom_smooth(se=FALSE)+
  theme(legend.position="bottom")
# native-feeding natural colonisers increase more recently

# Adventive colonisers
ad_col <- cdata[!cdata$coloniser_mode=="I",]
str(ad_col)
ad_col$host_plant_status <- factor(ad_col$host_plant_status, levels=c("Native", "Non-native"))

ad_col2 <- ad_col %>% group_by(mid_date, host_plant_status, .drop=FALSE) %>% summarise(nspp=n())

# Plot raw data and smoothed curves
ggplot(ad_col2, aes(mid_date, nspp, col=host_plant_status))+
  geom_line()+ geom_point()+
  geom_smooth(se=FALSE)+
  theme(legend.position="bottom")
## much fewer differences here

# Adventive vs Natural colonisers (regardless of host plant status)
nat_ad_col <- cdata
nat_ad_col$host_plant_status <- NULL

nat_ad_col2 <- nat_ad_col %>% group_by(mid_date, coloniser_mode, .drop=FALSE) %>% summarise(nspp=n())

# Plot raw data and smoothed curves
ggplot(nat_ad_col2, aes(mid_date, nspp, col=coloniser_mode))+
  geom_line()+ geom_point()+
  geom_smooth(se=FALSE)+
  theme(legend.position="bottom")
## both groups increased at similar rates, but fewer adventives

## First thoughts:
  # Initial increase in native colonisers overall driven by native and non-native feeders
  # Recent increase in native colonisers being driven by native feeders only = climate change
  # Initial increase in non-native feeding natural colonisers (probably human introducing ornamental plants), but recent decline/stable
  # Steady and more recent rapid increase in native feeding natural colonisers - due to climate change
  # No differences between native/non-native feeding adventives - likely driven by same mechanism?
  # Adventives overall also show rapid recent increase - maybe due to poor biosecurity?


#### TOTAL NUMBER OF NATURAL VS ADVENTIVE COLONISTS #####

# Q: Is there a difference in trend over time between natural and adventive colonists? 
# Do this for shifting start date periods as this might change over time - i.e. more similar at the end
nat_ad_col2 <- as.data.frame(nat_ad_col2)
nat_ad_col2$date_num <- case_when(
  nat_ad_col2$mid_date==1905 ~ 1,
  nat_ad_col2$mid_date==1915 ~ 2,
  nat_ad_col2$mid_date==1925 ~ 3,
  nat_ad_col2$mid_date==1935 ~ 4,
  nat_ad_col2$mid_date==1945 ~ 5,
  nat_ad_col2$mid_date==1955 ~ 6,
  nat_ad_col2$mid_date==1965 ~ 7,
  nat_ad_col2$mid_date==1975 ~ 8,
  nat_ad_col2$mid_date==1985 ~ 9,
  nat_ad_col2$mid_date==1995 ~ 10,
  nat_ad_col2$mid_date==2005 ~ 11,
  TRUE ~ 12
)

# early: 1905 - 1955
# late: 1955 - 2015
# long: 1905 - 2015
results <- rates <- coefs <- NULL
for(startdate in unique(c(1905, 1955))[1]){print(startdate)}

for(enddate in unique(c(1955,2015))[2]){print(enddate)}
  if(startdate==enddate) next 
# Subset by start and end decade
cdatat <- nat_ad_col2[nat_ad_col2$mid_date >= startdate & nat_ad_col2$mid_date <= enddate,]

  # Model fits without and with the interaction term
library(glmmTMB)
  fit1 <- glm(nspp ~ date_num+coloniser_mode, data = cdatat, family = "poisson")
  fit2 <- glm(nspp ~ coloniser_mode*date_num, data = cdatat, family="poisson")
  
  simulationOutput <- simulateResiduals(fittedModel = fit2)
  plot(simulationOutput)
  plotResiduals(simulationOutput, cdatat$date_num, quantreg = T)
  
  # Quantile deviations for: 1905, 1915, 1945
  # odd that this is worse for earlier years with more data
  # but testQuantiles p-value non-sig for 1905 and 1915
  # something weird going on with 1945
  
  # Significant dispersion test and quantile devations for 1905 - 1955 (underdispersion)
  # Dispersion fixed using generalized poisson, but still quantile deviations
  
  testDispersion(simulationOutput)
  summary(fit1)
  # intercept is the predicted nspp when mid_date=0 in the reference group of mode adventives 
  # mid_date - the simple slope of mid_date for reference group adventives (i.e. slope over time for adventives)
  # coloniser_modeI - the simple effect of mode OR the difference in nspp between immigrants and adventives at mid_date = 0
  
  # keeping coloniser mode constant (i.e. adventives), one unit increase in mid_year is associated with 0.01 increase in number of species
  # keeping the value of mid_year constant (at zero), the average number of species is 0.74 higher for immigrants than adventives
  
  summary(fit2)
  # same as above, but
  # coloniser_modeI:mid_date - the difference in the simple slopes of mid_date for immigrants vs adventives
  # slope for immigrants = slope of adventives + interaction slope
  # significant interaction = significant difference in mid_date slope between immigrants and adventives - NOT that each simple slope is different from zero
  
  # Use anova to check whether the interaction term is significant
  anova12 <- anova(fit1, fit2, test = "Chisq")
  
  AIC(fit1)
  AIC(fit2)
  
  # no difference in models for startdate: 1905, 1915, 1925, 1935, 1945 and 1955
  # interaction model better for startdate: none
  
  # Save the coefficients from the model with interaction term
  fit2output <- coef(summary(fit2))
  
  # Main results
  results <- rbind(results, data.frame(Startdate = startdate, # start date of model
                                       n_decades = uniqueN(cdatat$date_num), # number of decades in model
                                       anova_p = anova12$`Pr(>Chi)`[2], # p-value from anova - are the two models (with/without interaction) significantly different?
                                       diff_estimate = coef(summary(fit2))[4,"Estimate"], # estimate from interaction term i.e the slope - doesn't mean much on its own 
                                       diff_SE = coef(summary(fit2))[4,"Std. Error"], # standard error of interaction term
                                       diff_pval = coef(summary(fit2))[4,"Pr(>|z|)"], # p-value of interaction term
                                       od2 = dispersiontest(fit2)$estimate, # dispersion test estimate
                                       od2_pval = dispersiontest(fit2)$p.value)) # dispersion test p-value - is there any sign of overdispersion?
  
  dat <- ggpredict(fit2, terms = c("date_num", "coloniser_mode"))
  plot(dat) # very similar predicted rates over time from 1905
  
  # Save the rate estimates for native and non-native
  rates <- rbind(rates, data.frame(Startdate = startdate, # start date of model
                                   rate_adventive = exp(fit2output["date_num","Estimate"]), # exponential of the adventive slope
                                   rate_adventive_SE = 
                                     deltamethod(~exp(x1), 
                                                 fit2output["date_num","Estimate"], 
                                                 fit2output["date_num","Std. Error"]^2, ses=TRUE), # standard error around exponential of adventive slope
                                   rate_immigrant = exp(fit2output["date_num","Estimate"] + 
                                                          fit2output["coloniser_modeI:date_num","Estimate"]), # exponential of the immigrant slope (=adventive slope + interaction slope)
                                   rate_immigrant_SE = deltamethod(~exp(x1+x2), fit2output[c("date_num","coloniser_modeI:date_num"),"Estimate"],
                                                                   vcov(fit2)[3:4, 3:4], ses=TRUE))) # standard error around exponential of immigrant slope
 
  # Also look at estimates for native/non-native on date
  coefs <- rbind(coefs, data.frame(Startdate = startdate,
                                   n_decades = uniqueN(cdatat$date_num),
                                   adventive_est = fit2output["date_num", "Estimate"], # slope for reference level (adventives)
                                   adventive_SE = fit2output["date_num", "Std. Error"], # standard error for adventive slope
                                   immigrant_est = sum(fit2output[3:4, "Estimate"]), # slope for immigrant (=adventive slope + interaction slope)
                                   immigrant_SE = sqrt(sum(fit2output[3:4, "Std. Error"]^2)+ # standard error for immigrant slope
                                                         2*vcov(fit2)[3, 4])))
  
  
}

# Produce confidence intervals for the rate estimates
setDT(rates)
setDT(results)
setDT(coefs)
rates[, rate_adventive_lower := rate_adventive - 1.96*rate_adventive_SE]
rates[, rate_adventive_upper := rate_adventive + 1.96*rate_adventive_SE]
rates[, rate_immigrant_lower := rate_immigrant - 1.96*rate_immigrant_SE]
rates[, rate_immigrant_upper := rate_immigrant + 1.96*rate_immigrant_SE]
# non-significant if confidence intervals contain 1


#### NATURAL NATIVE-FEEDERS VS NON-NATIVE FEEDERS #####

# Q: Is there a difference in trend over time between native and non-native feeder natural colonisers? 
# Do this for shifting start date periods as this might change over time - i.e. more similar at the beginning

nat_col2 <- as.data.frame(nat_col2)
nat_col2$date_num <- case_when(
  nat_col2$mid_date==1905 ~ 1,
  nat_col2$mid_date==1915 ~ 2,
  nat_col2$mid_date==1925 ~ 3,
  nat_col2$mid_date==1935 ~ 4,
  nat_col2$mid_date==1945 ~ 5,
  nat_col2$mid_date==1955 ~ 6,
  nat_col2$mid_date==1965 ~ 7,
  nat_col2$mid_date==1975 ~ 8,
  nat_col2$mid_date==1985 ~ 9,
  nat_col2$mid_date==1995 ~ 10,
  nat_col2$mid_date==2005 ~ 11,
  TRUE ~ 12
)

results <- rates <- coefs <- NULL
for(startdate in seq(1905, 1965, by=10)[7]){print(startdate)}

for(enddate in unique(c(1955,2015))[2]){print(enddate)}
if(startdate==enddate) next 

# Subset by start and end decade
cdatat <- nat_col2[nat_col2$mid_date >= startdate & nat_col2$mid_date <= enddate,]
  
  # Model fits without and with the interaction term
  fit1 <- glm(nspp ~ date_num+host_plant_status, data = cdatat, family = "poisson")
  fit2 <- glm(nspp ~ host_plant_status*date_num, data = cdatat, family="poisson")
  
  simulationOutput <- simulateResiduals(fittedModel = fit1)
  plot(simulationOutput)
  
  # Quantile deviations for: 1945 (something weird going on here) - fit1 is fine though
  # but fit2 has much lower AIC (nearly 6)

  summary(fit1)
  # intercept is the predicted nspp when mid_date=0 in the reference group of mode adventives 
  # mid_date - the simple slope of mid_date for reference group adventives (i.e. slope over time for adventives)
  # coloniser_modeI - the simple effect of mode OR the difference in nspp between immigrants and adventives at mid_date = 0
  
  # keeping coloniser mode constant (i.e. adventives), one unit increase in mid_year is associated with 0.01 increase in number of species
  # keeping the value of mid_year constant (at zero), the average number of species is 0.74 higher for immigrants than adventives
  
  summary(fit2)
  # same as above, but
  # coloniser_modeI:mid_date - the difference in the simple slopes of mid_date for immigrants vs adventives
  # slope for immigrants = slope of adventives + interaction slope
  # significant interaction = significant difference in mid_date slope between immigrants and adventives - NOT that each simple slope is different from zero
  
  # Use anova to check whether the interaction term is significant
  anova12 <- anova(fit1, fit2, test = "Chisq")
  
  AIC(fit1)
  AIC(fit2)
  
  # no difference in models for startdate: 1905, 1915, 
  # interaction model better for startdate: 1925, 1935
  
  # Save the coefficients from the model with interaction term
  fit2output <- coef(summary(fit2))
  
  # Main results
  results <- rbind(results, data.frame(Startdate = startdate, # start date of model
                                       n_decades = uniqueN(cdatat$date_num), # number of decades in model
                                       anova_p = anova12$`Pr(>Chi)`[2], # p-value from anova - are the two models (with/without interaction) significantly different?
                                       diff_estimate = coef(summary(fit2))[4,"Estimate"], # estimate from interaction term i.e the slope - doesn't mean much on its own 
                                       diff_SE = coef(summary(fit2))[4,"Std. Error"], # standard error of interaction term
                                       diff_pval = coef(summary(fit2))[4,"Pr(>|z|)"], # p-value of interaction term
                                       od2 = dispersiontest(fit2)$estimate, # dispersion test estimate
                                       od2_pval = dispersiontest(fit2)$p.value)) # dispersion test p-value - is there any sign of overdispersion?
  
  dat <- ggpredict(fit2, terms = c("date_num", "host_plant_status"))
  plot(dat) # very similar predicted rates over time from 1905
  
  # Save the rate estimates for native and non-native
  rates <- rbind(rates, data.frame(Startdate = startdate, # start date of model
                                   rate_native = exp(fit2output["date_num","Estimate"]), # exponential of the native slope
                                   rate_native_SE = 
                                     deltamethod(~exp(x1), 
                                                 fit2output["date_num","Estimate"], 
                                                 fit2output["date_num","Std. Error"]^2, ses=TRUE), # standard error around exponential of native slope
                                   rate_nonnative = exp(fit2output["date_num","Estimate"] + 
                                                          fit2output["host_plant_statusNon-native:date_num","Estimate"]), # exponential of the non-native slope (=native slope + interaction slope)
                                   rate_nonnative_SE = deltamethod(~exp(x1+x2), fit2output[c("date_num","host_plant_statusNon-native:date_num"),"Estimate"],
                                                                   vcov(fit2)[3:4, 3:4], ses=TRUE))) # standard error around exponential of non-native slope
  
  # Also look at estimates for native/non-native on date
  coefs <- rbind(coefs, data.frame(Startdate = startdate,
                                   n_decades = uniqueN(cdatat$date_num),
                                   native_est = fit2output["date_num", "Estimate"], # slope for reference level (natives)
                                   native_SE = fit2output["date_num", "Std. Error"], # standard error for native slope
                                   nonnative_est = sum(fit2output[3:4, "Estimate"]), # slope for non-native (=native slope + interaction slope)
                                   nonnative_SE = sqrt(sum(fit2output[3:4, "Std. Error"]^2)+ # standard error for non-native slope
                                                         2*vcov(fit2)[3, 4])))
  
  
}

# Produce confidence intervals for the rate estimates
setDT(rates)
setDT(results)
setDT(coefs)
rates[, rate_native_lower := rate_native - 1.96*rate_native_SE]
rates[, rate_native_upper := rate_native + 1.96*rate_native_SE]
rates[, rate_nonnative_lower := rate_nonnative - 1.96*rate_nonnative_SE]
rates[, rate_nonnative_upper := rate_nonnative + 1.96*rate_nonnative_SE]
# non-significant if confidence intervals contain 1


#### ADVENTIVE NATIVE-FEEDERS VS NON-NATIVE FEEDERS #####

# Q: Is there a difference in trend over time between native and non-native feeder adventive colonisers? 
# Do this for shifting start date periods as this might change over time - i.e. more similar at the beginning

ad_col2 <- as.data.frame(ad_col2)
ad_col2$date_num <- case_when(
  ad_col2$mid_date==1905 ~ 1,
  ad_col2$mid_date==1915 ~ 2,
  ad_col2$mid_date==1925 ~ 3,
  ad_col2$mid_date==1935 ~ 4,
  ad_col2$mid_date==1945 ~ 5,
  ad_col2$mid_date==1955 ~ 6,
  ad_col2$mid_date==1965 ~ 7,
  ad_col2$mid_date==1975 ~ 8,
  ad_col2$mid_date==1985 ~ 9,
  ad_col2$mid_date==1995 ~ 10,
  ad_col2$mid_date==2005 ~ 11,
  TRUE ~ 12
)

results <- rates <- coefs <- NULL
for(startdate in seq(1905, 1955, by=10)[6]){print(startdate)}
# Subset by start decade
cdatat <- ad_col2[ad_col2$mid_date >= startdate,]

# Model fits without and with the interaction term
fit1 <- glm(nspp ~ date_num+host_plant_status, data = cdatat, family = "poisson")
fit2 <- glm(nspp ~ host_plant_status*date_num, data = cdatat, family="poisson")

simulationOutput <- simulateResiduals(fittedModel = fit1)
plot(simulationOutput)

# No assumptions violated for all models 1905-1955 

summary(fit1)
# intercept is the predicted nspp when mid_date=0 in the reference group of mode adventives 
# mid_date - the simple slope of mid_date for reference group adventives (i.e. slope over time for adventives)
# coloniser_modeI - the simple effect of mode OR the difference in nspp between immigrants and adventives at mid_date = 0

# keeping coloniser mode constant (i.e. adventives), one unit increase in mid_year is associated with 0.01 increase in number of species
# keeping the value of mid_year constant (at zero), the average number of species is 0.74 higher for immigrants than adventives

summary(fit2)
# same as above, but
# coloniser_modeI:mid_date - the difference in the simple slopes of mid_date for immigrants vs adventives
# slope for immigrants = slope of adventives + interaction slope
# significant interaction = significant difference in mid_date slope between immigrants and adventives - NOT that each simple slope is different from zero

# Use anova to check whether the interaction term is significant
anova12 <- anova(fit1, fit2, test = "Chisq")

AIC(fit1)
AIC(fit2)

# fit1 has lower AIC for: 1905, 1915, 1925, 1935, 1945

# Save the coefficients from the model with interaction term
fit2output <- coef(summary(fit2))

# Main results
results <- rbind(results, data.frame(Startdate = startdate, # start date of model
                                     n_decades = uniqueN(cdatat$date_num), # number of decades in model
                                     anova_p = anova12$`Pr(>Chi)`[2], # p-value from anova - are the two models (with/without interaction) significantly different?
                                     diff_estimate = coef(summary(fit2))[4,"Estimate"], # estimate from interaction term i.e the slope - doesn't mean much on its own 
                                     diff_SE = coef(summary(fit2))[4,"Std. Error"], # standard error of interaction term
                                     diff_pval = coef(summary(fit2))[4,"Pr(>|z|)"], # p-value of interaction term
                                     od2 = dispersiontest(fit2)$estimate, # dispersion test estimate
                                     od2_pval = dispersiontest(fit2)$p.value)) # dispersion test p-value - is there any sign of overdispersion?

dat <- ggpredict(fit1, terms = "date_num")
plot(dat) # very similar predicted rates over time from 1905

# Save the rate estimates for native and non-native
rates <- rbind(rates, data.frame(Startdate = startdate, # start date of model
                                 rate_native = exp(fit2output["date_num","Estimate"]), # exponential of the native slope
                                 rate_native_SE = 
                                   deltamethod(~exp(x1), 
                                               fit2output["date_num","Estimate"], 
                                               fit2output["date_num","Std. Error"]^2, ses=TRUE), # standard error around exponential of native slope
                                 rate_nonnative = exp(fit2output["date_num","Estimate"] + 
                                                        fit2output["host_plant_statusNon-native:date_num","Estimate"]), # exponential of the non-native slope (=native slope + interaction slope)
                                 rate_nonnative_SE = deltamethod(~exp(x1+x2), fit2output[c("date_num","host_plant_statusNon-native:date_num"),"Estimate"],
                                                                 vcov(fit2)[3:4, 3:4], ses=TRUE))) # standard error around exponential of non-native slope

# Also look at estimates for native/non-native on date
coefs <- rbind(coefs, data.frame(Startdate = startdate,
                                 n_decades = uniqueN(cdatat$date_num),
                                 native_est = fit2output["date_num", "Estimate"], # slope for reference level (natives)
                                 native_SE = fit2output["date_num", "Std. Error"], # standard error for native slope
                                 nonnative_est = sum(fit2output[3:4, "Estimate"]), # slope for non-native (=native slope + interaction slope)
                                 nonnative_SE = sqrt(sum(fit2output[3:4, "Std. Error"]^2)+ # standard error for non-native slope
                                                       2*vcov(fit2)[3, 4])))


}

# Produce confidence intervals for the rate estimates
setDT(rates)
setDT(results)
setDT(coefs)
rates[, rate_native_lower := rate_native - 1.96*rate_native_SE]
rates[, rate_native_upper := rate_native + 1.96*rate_native_SE]
rates[, rate_nonnative_lower := rate_nonnative - 1.96*rate_nonnative_SE]
rates[, rate_nonnative_upper := rate_nonnative + 1.96*rate_nonnative_SE]
# non-significant if confidence intervals contain 1


set.seed(12)
xx <- 1:100
zz <- runif(100)
yy <- 2 + 1.5*pmax(xx - 35, 0) - 1.5*pmax(xx - 70, 0) + 15*pmax(zz - .5, 0) + 
  rnorm(100,0,2)
dati <- data.frame(x = xx, y = yy, z = zz)
out.lm <- lm(y ~ x, data = dati)
o <- segmented(out.lm, seg.Z = ~x, psi = list(x = c(30,60)),
               control = seg.control(display = FALSE)
)
dat2 = data.frame(x = xx, y = broken.line(o)$fit)

library(ggplot2)
ggplot(dati, aes(x = x, y = y)) +
  geom_point() +
  geom_line(data = dat2, color = 'blue')
