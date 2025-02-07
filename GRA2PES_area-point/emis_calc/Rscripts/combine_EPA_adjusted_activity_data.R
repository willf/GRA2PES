#install libraries
#install.packages("plyr")
library(plyr)
#install.packages("dplyr")
library(dplyr)
#install.packages("readr")
library(readr)

#set working directory
setwd ("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/emis_calc/input")

#read prepared .csv files
eia_activity_data_new <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/emis_calc/input/use_all_btu_new_unadjusted_bBtu.csv") #US in this file means 50 states + DC

epa_activity_data_res <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/emis_calc/input/EPA_activity_data_RES_tBtu.csv") 
epa_activity_data_comm <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/emis_calc/input/EPA_activity_data_COMM_tBtu.csv") 
epa_activity_data_ind <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/emis_calc/input/EPA_activity_data_IND_tBtu.csv")

#use EIA new version and substitute EIA number with EPA adjusted number if available
colnames(epa_activity_data_res)[colnames(epa_activity_data_res) == "X2017_adj"] ="adj_2017"
colnames(epa_activity_data_comm)[colnames(epa_activity_data_comm) == "X2017_adj"] ="adj_2017"
colnames(epa_activity_data_ind)[colnames(epa_activity_data_ind) == "X2017_adj"] ="adj_2017"

#stack the three sectors
epa_activity_data = rbind(epa_activity_data_res, epa_activity_data_comm, epa_activity_data_ind)

#convert trillion Btu to billion Btu for epa_activity_data_*, 1 trillion = 1000 billion
epa_activity_data$adj_2017 = epa_activity_data$adj_2017 * 1000 #bBtu

#join with eia unadjusted data
eia_activity_data_adjusted = join(eia_activity_data_new,epa_activity_data, by=c("State","MSN"), type="left")

#fill in values in adj_2017 if it is NA
for (rr in 1:nrow(eia_activity_data_adjusted)) {
  rr_adj_2017 <- eia_activity_data_adjusted[rr, "adj_2017"]
  if (is.na(rr_adj_2017)){
    eia_activity_data_adjusted[rr, "adj_2017"] = eia_activity_data_adjusted[rr, "X2017"]
  }
}

eia_activity_data_adjusted_out = subset(eia_activity_data_adjusted, select=c("State","MSN","adj_2017"))

write_csv(eia_activity_data_adjusted_out, file = "use_all_btu_new_adjusted_bBtu.csv")
