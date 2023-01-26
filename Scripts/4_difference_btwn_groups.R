##########################
#### user: Lisbeth Hordley
#### date: June 2022
#### info: Difference in trends between coloniser modes + host plant status

# Packages
library(data.table)
library(ggplot2)
library(DHARMa)
library(dplyr)
library(tidyverse)

## Group/time interaction
cdata <- read.csv("Coloniser data/moth_coloniser_data_spp.csv", header=TRUE)

cdata$mid_date <- factor(cdata$mid_date, levels=c("1905","1915","1925","1935","1945","1955","1965","1975","1985","1995","2005", "2015"))
cdata$coloniser_mode <- factor(cdata$coloniser_mode, levels=c("A", "I"))
cdata$host_plant_status <- factor(cdata$host_plant_status, levels=c("Native", "Non-native"))

# Immigrants vs Adventives - calculate number of species with each colonisation mode
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

# GLM
mod <- glm(nspp ~ date_num*coloniser_mode, data=imm_ad_spp, family="poisson")
summary(mod) # non-significant 

# Immigrant native vs non-native
cdata2 <- cdata[!is.na(cdata$host_plant_status),] # 139 species
cdata2$mid_date <- factor(cdata2$mid_date, levels=c("1905","1915","1925","1935","1945","1955","1965","1975","1985","1995","2005","2015"))
cdata2$host_plant_status <- factor(cdata2$host_plant_status, levels=c("Native", "Non-native"))

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

# GLM
mod2 <- glm(nspp ~ date_num*host_plant_status, data=imm_natnn_spp, family="poisson")
summary(mod2) # non-significant 

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

# GLM
mod3 <- glm(nspp ~ date_num*host_plant_status, data=ad_natnn_spp, family="poisson")
summary(mod3) # non-significant 
