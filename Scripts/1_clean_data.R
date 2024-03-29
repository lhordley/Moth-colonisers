##########################
#### user: Lisbeth Hordley
#### date: June 2022
#### info: Clean + sort moth data

# Packages
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)
options(scipen=999)

# Load data
cdata <- read.csv("Coloniser data/moth_coloniser_data_spp_new.csv", header=TRUE) # Categories based on Parsons 2003, 2010, 2020

# Change columns to factors
cdata$mid_date <- factor(cdata$mid_date, levels=c("1905","1915","1925","1935","1945","1955","1965","1975","1985","1995","2005", "2015"))
cdata$coloniser_mode <- factor(cdata$coloniser_mode, levels=c("A", "I"))
cdata$host_plant_status <- factor(cdata$host_plant_status, levels=c("Native", "Non-native"))

# Calculate number of species in each group for each decade
cdata_allspp <- cdata %>% group_by(mid_date, .drop=FALSE) %>% summarise(all_nspp=n())
cdata_immigrants <- cdata %>% filter(coloniser_mode=="I") %>% group_by(mid_date, .drop=FALSE) %>% summarise(imm_nspp=n())
cdata_adventives <- cdata %>% filter(coloniser_mode=="A") %>% group_by(mid_date, .drop=FALSE) %>% summarise(adv_nspp=n())
cdata_imm_native <- cdata %>% filter(coloniser_mode=="I", host_plant_status=="Native") %>% group_by(mid_date, .drop=FALSE) %>% 
  summarise(imm_nat_nspp=n())
cdata_imm_nonnative <- cdata %>% filter(coloniser_mode=="I", host_plant_status=="Non-native") %>% group_by(mid_date, .drop=FALSE) %>% 
  summarise(imm_nn_nspp=n())
cdata_adv_native <- cdata %>% filter(coloniser_mode=="A", host_plant_status=="Native") %>% group_by(mid_date, .drop=FALSE) %>% 
  summarise(adv_nat_nspp=n())
cdata_adv_nonnative <- cdata %>% filter(coloniser_mode=="A", host_plant_status=="Non-native") %>% group_by(mid_date, .drop=FALSE) %>% 
  summarise(adv_nn_nspp=n())

cdata_final <- list(cdata_allspp, cdata_immigrants, cdata_adventives, cdata_imm_native, cdata_imm_nonnative, cdata_adv_native, cdata_adv_nonnative) %>% reduce(full_join, by = "mid_date")
cdata_final <- as.data.frame(cdata_final)
cdata_final$date_num <- case_when(
  cdata_final$mid_date==1905 ~ 1,
  cdata_final$mid_date==1915 ~ 2,
  cdata_final$mid_date==1925 ~ 3,
  cdata_final$mid_date==1935 ~ 4,
  cdata_final$mid_date==1945 ~ 5,
  cdata_final$mid_date==1955 ~ 6,
  cdata_final$mid_date==1965 ~ 7,
  cdata_final$mid_date==1975 ~ 8,
  cdata_final$mid_date==1985 ~ 9,
  cdata_final$mid_date==1995 ~ 10,
  cdata_final$mid_date==2005 ~ 11,
  TRUE ~ 12
)
write.csv(cdata_final, file="Coloniser data/moth_coloniser_nspp.csv", row.names=FALSE)


