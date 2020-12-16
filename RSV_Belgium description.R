
rm(list = ls())
setwd("/home/phuong/phuonght/gitKraken_project/RSV")
library(deSolve)
library(ggplot2)
library(scales)

#BELGIUM DATA-----------------
bel_data <- read.csv("./RSV data/RSV_cases_time_epistat.csv")
str(bel_data)
bel_data$date <- as.Date(bel_data$date, format = "%Y-%m-%d")

#Reported cases 2011-2018=====
# total cases: 63301 , less than 2yrs: 56126, 1yr:42375,2yrs:13751
pro_case_less2yrs <- 56126/63301
pro_case_2yrs <- 13751/63301
pro_case_1yr <- 42375/63301

bel_data$week <- seq(1,dim(bel_data)[1])
bel_data$week.yr <- as.numeric(strftime(bel_data$date, format = "%V"))# get week numbers from date
bel_data$year <- as.numeric(format(bel_data$date,"%Y")) # get year numbers from date
#reassign the year of week 53
bel_data$year[bel_data$week.yr==53] <- 2015 # start week:"2016-01-02"

bel_data$epi.period[bel_data$week.yr < 36] = paste0(bel_data$year[bel_data$week.yr <36]-1,"-",bel_data$year[bel_data$week.yr  <36])
bel_data$epi.period[bel_data$week.yr >= 36] = paste0(bel_data$year[bel_data$week.yr >=36 ],"-",bel_data$year[bel_data$week.yr >=36 ]+1)
bel_data$epi.period <- as.factor(bel_data$epi.period)
# week.yr = 36 ~ week.epi.period = 1 
bel_data$week.epi.period[bel_data$epi.period=="2010-2011"] <- c(18:52)# nodat in 2010
bel_data$week.epi.period[bel_data$epi.period=="2011-2012"] <- c(1:52)
bel_data$week.epi.period[bel_data$epi.period=="2012-2013"] <- c(1:52)
bel_data$week.epi.period[bel_data$epi.period=="2013-2014"] <- c(1:52)
bel_data$week.epi.period[bel_data$epi.period=="2014-2015"] <- c(1:52)
bel_data$week.epi.period[bel_data$epi.period=="2015-2016"] <- c(1:52)
bel_data$week.epi.period[bel_data$epi.period=="2016-2017"] <- c(1:53)#53 wks
bel_data$week.epi.period[bel_data$epi.period=="2017-2018"] <- c(1:30)# wk13 2018



# bel_data$cases.less2yrs <- round(bel_data$cases*pro_case_less2yrs)
bel_data$cases.1yr <- round(bel_data$cases*pro_case_1yr)
bel_data$cases.2yrs <- round(bel_data$cases*pro_case_2yrs)

break.vec <- c(as.Date("2011-01-08"),
               seq(from=as.Date("2011-01-08"), to=as.Date("2018-03-31"),
                   by="16 weeks"))
ggplot(bel_data,aes(x = date, y= cases))+
         geom_line(color ="steelblue")+
  # scale_color_manual(name = "", values = c("red","blue","black"))+
  # ggtitle (" Weekly RSV reported in Belgium")+
  scale_x_date(breaks = break.vec)+
  xlab(" Date of Diagnosis")+
  ylab(" Number of RSV cases ")+
  theme(plot.title =element_text(hjust = 0.5,size=22,face="bold"),
        title=element_text(size=15),
        strip.text.x = element_text(size=18),
        axis.text.x = element_text(angle = 60,face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"))

#Reported cases every year=====
ggplot(bel_data, aes(x = week.epi.period)) + 
  geom_line(aes(y = cases, fill=epi.period, color= epi.period)) +
  scale_x_continuous(breaks = c(1,5,9,13,17,21,25,29,33,37,41,45,49),
                   labels = c(36,40,44,48,52,4,8,12,16,20,24,28,32))+
  xlab("Weeks")+
  ylab("Confirmed RSV cases")
  theme(plot.title =element_text(hjust = 0.5,size=22,face="bold"),
        title=element_text(size=15),
        strip.text.x = element_text(size=12, face = "bold"),
        # axis.text.x = element_text(angle = 90, value = ()),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "bold",size = 12),
        legend.text = element_text(face = "bold",size = 16))

