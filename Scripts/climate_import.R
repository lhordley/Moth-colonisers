##########################
#### user: Lisbeth Hordley
#### date: June 2022
#### info: Decadal temperature trend

# Packages
library(ggplot2)
library(dplyr)

setwd("~/Moth colonisers")

# Load data
temp <- read.csv("HadCET_annual_temperature.csv", header=TRUE)

temp$decade <- rep(c("1905","1915","1925","1935","1945","1955","1965","1975","1985","1995","2005","2015"),each=10)

# calculate mean temperature per decade
decade_temp <- temp %>% group_by(decade) %>% summarise(mean_temp=mean(temperature))
decade_temp$decade <- as.integer(decade_temp$decade)
temp <- ggplot(decade_temp, aes(x=decade, y=mean_temp))+
  geom_line(col="black", lwd=1)+
  geom_line(data=temp, aes(x=ï..year, y=temperature), col="grey")+
  labs(x="Year", y="Mean decadal temperature")+
  scale_x_continuous(breaks=seq(1905,2015,by=20))+
  scale_y_continuous(breaks=seq(8,11,by=0.2))+
  theme_classic()
temp
ggsave(temp, file="decadal_temperature.png")

imm_temp <- merge(imm_spp, decade_temp, by.x="mid_date", by.y="decade")

cor.test(imm_temp$nspp, imm_temp$mean_temp, method="spearman")
# significant positive correlation
plot(nspp ~ mean_temp, data=imm_temp)

adv_temp <- merge(adv_spp, decade_temp, by.x="mid_date", by.y="decade")
cor.test(adv_temp$nspp, adv_temp$mean_temp, method="pearson")
# significant positive correlation
plot(nspp ~ mean_temp, data=adv_temp)



## 

# Trade import values

import <- read.csv("CorrelatesOfWar_UK_trade_imports.csv", header=TRUE)
# only have data to 2014 - remove last 5 years
colnames(import)[1] <- "year"
import <- import[import$year<=2009,]
import$decade <- rep(c("1905","1915","1925","1935","1945","1955","1965","1975","1985","1995","2005"),each=10)

# calculate mean temperature per decade
decade_import <- import %>% group_by(decade) %>% summarise(mean_import=mean(imports))
decade_import$decade <- as.integer(decade_import$decade)
import <- ggplot(decade_import, aes(x=decade, y=mean_import))+
  geom_line(col="black", lwd=1)+
  geom_line(data=import, aes(x=year, y=imports), col="grey")+
  labs(x="Year", y="Mean decadal imports to the UK")+
  scale_x_continuous(breaks=seq(1905,2005,by=20))+
  #scale_y_continuous(breaks=seq(8,11,by=0.2))+
  theme_classic()
import
ggsave(import, file="decadal_imports.png")

