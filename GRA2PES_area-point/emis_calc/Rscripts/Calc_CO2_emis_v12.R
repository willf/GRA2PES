#script to join MSN code, CO2 emissions factor with EIA data
#v1 change from v0
#1. got rid of motor gasoline MSN code from commercial and industrial
#2. No calculations for transportation sector anymore
#3. Divide annual total emissions of CO2 (metric tons) into weekday, Saturday, and Sunday before applying point/area ratio
#3. Not using weekday, Saturday, and Sunday temporally weighted average point/area CO emissions to calculate annual average point/area ratio
#4. using weekday, Saturday, Sunday, day-specific point/area CO emissions to calculate day-specific point/area ratio

#v2 change from v1
#Included MSN: HLICB, POINTmyCODE:PtIND_NG2. AREAmyCODE:IND_NG

#v3 change from v2
#use EIA emissions factors and BTU based MSN and usage data

#v4 change from v3
#remove wood and waste MSN code, POINTmyCODE, AREAmyCODE

#v5 change from v4
#changed directories to match new structure of the project folder
#added more comments to improve readability

#v7 change from v5
#use EPA adjusted activity data
#remove 2018 and 2019

#V8 created separate script for IND to remove matched GHGRP points, the RES, COMM part of the code is the same as V7

#v9 merge separate IND script into RES COMM script, add more EPA fuel types and do NEU correction same as EPA GHGI

#v10 
#do not subtract GHGRP ref, chem_mfg, and ng_proc points anymore
#improved code that was hard coded
#Add the choice of using sum_CO in specific fuel or lumped fuel style, use as much specific fuel sum_CO as possible, use lumped fuel as patches

#v11
#subtract FC emissions at state-level for GHGRP refineries, chemicals, minerals and metals
#only patch when state-level emissions are nonzero

#v12
#subtract fuel-specific GHGRP at state-level
#changed a couple MSN to AREAmyCODE mapping in MSN_to_AREAmyCODE: Petroleum Coke to Coal; Still Gas to Natural Gas

#install libraries
#install.packages("plyr")
library(plyr)
#install.packages("dplyr")
library(dplyr)
#install.packages("readr")
library(readr)

#set working directory
setwd ("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/emis_calc/output")

#read prepared .csv files
MSN_co2_fct <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/emis_calc/input/Selected_MSN_CO2EF_MTCO2_per_bBtu.csv")
btu_consump_states <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/emis_calc/input/use_all_btu_new_adjusted_bBtu.csv") #US in this file means 50 states + DC
AREA_sectors_CO <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Colin_SS_sumCO/sum_CO/AREA_sectors_CO_specificfuel.csv") #US in this file means US territories outside of CONUS
POINT_sectors_CO <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Colin_SS_sumCO/sum_CO/POINT_sectors_CO_specificfuel.csv") #US in this file means US territories outside of CONUS
AREA_sectors_CO_patch <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Colin_SS_sumCO/sum_CO/AREA_sectors_CO_lumpfuel.csv") #US in this file means US territories outside of CONUS
POINT_sectors_CO_patch <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Colin_SS_sumCO/sum_CO/POINT_sectors_CO_lumpfuel.csv") #US in this file means US territories outside of CONUS

#So, US (50 states + DC) total emissions are divided into point/area by point/area ratios of US territories, this is wrong. So, do not use US entry in the product. That is exactly what Colin did, so no problems!
#PtIND_NG2 is combined with PtIND_NG in the R code,
#PtIND_Oil2 is combined with PtIND_Oil in the R code,
#"STATE_FIPS" in POINT_sectors_CO are different from AREA_sectors_CO, do not use it to match point with area.
State_fullname_to_abb <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/emis_calc/input/State_fullname_to_abb.csv")
MSN_to_AREAmyCODE <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/emis_calc/input/MSN_to_AREAmyCODE.csv")
POINTmyCODE_to_AREAmyCODE <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/emis_calc/input/POINTmyCODE_to_AREAmyCODE.csv")

#############################################################################################################################
#GHGRP CO2
GHGRP_CO2_refineries_Coal <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Match_GHGRP_NEI_IND/spreadsheets/input/GHGRP_refineries_2017_CO2_FC_Coal.csv")
GHGRP_CO2_refineries_NG <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Match_GHGRP_NEI_IND/spreadsheets/input/GHGRP_refineries_2017_CO2_FC_NG.csv")
GHGRP_CO2_refineries_Petroleum <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Match_GHGRP_NEI_IND/spreadsheets/input/GHGRP_refineries_2017_CO2_FC_Petroleum.csv")
GHGRP_CO2_refineries_Other <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Match_GHGRP_NEI_IND/spreadsheets/input/GHGRP_refineries_2017_CO2_FC_Other.csv")

GHGRP_CO2_chemicals_Coal <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Match_GHGRP_NEI_IND/spreadsheets/input/GHGRP_chemicals_2017_CO2_FC_Coal.csv")
GHGRP_CO2_chemicals_NG <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Match_GHGRP_NEI_IND/spreadsheets/input/GHGRP_chemicals_2017_CO2_FC_NG.csv")
GHGRP_CO2_chemicals_Petroleum <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Match_GHGRP_NEI_IND/spreadsheets/input/GHGRP_chemicals_2017_CO2_FC_Petroleum.csv")
GHGRP_CO2_chemicals_Other <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Match_GHGRP_NEI_IND/spreadsheets/input/GHGRP_chemicals_2017_CO2_FC_Other.csv")

GHGRP_CO2_minerals_metals_Coal <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Match_GHGRP_NEI_IND/spreadsheets/input/GHGRP_minerals_metals_2017_CO2_FC_Coal.csv")
GHGRP_CO2_minerals_metals_NG <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Match_GHGRP_NEI_IND/spreadsheets/input/GHGRP_minerals_metals_2017_CO2_FC_NG.csv")
GHGRP_CO2_minerals_metals_Petroleum <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Match_GHGRP_NEI_IND/spreadsheets/input/GHGRP_minerals_metals_2017_CO2_FC_Petroleum.csv")
GHGRP_CO2_minerals_metals_Other <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Match_GHGRP_NEI_IND/spreadsheets/input/GHGRP_minerals_metals_2017_CO2_FC_Other.csv")

colnames(GHGRP_CO2_refineries_Coal)[colnames(GHGRP_CO2_refineries_Coal) == "GHGRP.ID"] ="GHGRP_ID"
colnames(GHGRP_CO2_refineries_Coal)[colnames(GHGRP_CO2_refineries_Coal) == "GHG.QUANTITY..METRIC.TONS.CO2e."] ="CO2_GHGRP" #unit: metric tons/year
GHGRP_CO2_refineries_Coal = subset(GHGRP_CO2_refineries_Coal, select=c("STATE","GHGRP_ID","CO2_GHGRP"))
colnames(GHGRP_CO2_refineries_NG)[colnames(GHGRP_CO2_refineries_NG) == "GHGRP.ID"] ="GHGRP_ID"
colnames(GHGRP_CO2_refineries_NG)[colnames(GHGRP_CO2_refineries_NG) == "GHG.QUANTITY..METRIC.TONS.CO2e."] ="CO2_GHGRP" #unit: metric tons/year
GHGRP_CO2_refineries_NG = subset(GHGRP_CO2_refineries_NG, select=c("STATE","GHGRP_ID","CO2_GHGRP"))
colnames(GHGRP_CO2_refineries_Petroleum)[colnames(GHGRP_CO2_refineries_Petroleum) == "GHGRP.ID"] ="GHGRP_ID"
colnames(GHGRP_CO2_refineries_Petroleum)[colnames(GHGRP_CO2_refineries_Petroleum) == "GHG.QUANTITY..METRIC.TONS.CO2e."] ="CO2_GHGRP" #unit: metric tons/year
GHGRP_CO2_refineries_Petroleum = subset(GHGRP_CO2_refineries_Petroleum, select=c("STATE","GHGRP_ID","CO2_GHGRP"))
colnames(GHGRP_CO2_refineries_Other)[colnames(GHGRP_CO2_refineries_Other) == "GHGRP.ID"] ="GHGRP_ID"
colnames(GHGRP_CO2_refineries_Other)[colnames(GHGRP_CO2_refineries_Other) == "GHG.QUANTITY..METRIC.TONS.CO2e."] ="CO2_GHGRP" #unit: metric tons/year
GHGRP_CO2_refineries_Other = subset(GHGRP_CO2_refineries_Other, select=c("STATE","GHGRP_ID","CO2_GHGRP"))

colnames(GHGRP_CO2_chemicals_Coal)[colnames(GHGRP_CO2_chemicals_Coal) == "GHGRP.ID"] ="GHGRP_ID"
colnames(GHGRP_CO2_chemicals_Coal)[colnames(GHGRP_CO2_chemicals_Coal) == "GHG.QUANTITY..METRIC.TONS.CO2e."] ="CO2_GHGRP" #unit: metric tons/year
GHGRP_CO2_chemicals_Coal = subset(GHGRP_CO2_chemicals_Coal, select=c("STATE","GHGRP_ID","CO2_GHGRP"))
colnames(GHGRP_CO2_chemicals_NG)[colnames(GHGRP_CO2_chemicals_NG) == "GHGRP.ID"] ="GHGRP_ID"
colnames(GHGRP_CO2_chemicals_NG)[colnames(GHGRP_CO2_chemicals_NG) == "GHG.QUANTITY..METRIC.TONS.CO2e."] ="CO2_GHGRP" #unit: metric tons/year
GHGRP_CO2_chemicals_NG = subset(GHGRP_CO2_chemicals_NG, select=c("STATE","GHGRP_ID","CO2_GHGRP"))
colnames(GHGRP_CO2_chemicals_Petroleum)[colnames(GHGRP_CO2_chemicals_Petroleum) == "GHGRP.ID"] ="GHGRP_ID"
colnames(GHGRP_CO2_chemicals_Petroleum)[colnames(GHGRP_CO2_chemicals_Petroleum) == "GHG.QUANTITY..METRIC.TONS.CO2e."] ="CO2_GHGRP" #unit: metric tons/year
GHGRP_CO2_chemicals_Petroleum = subset(GHGRP_CO2_chemicals_Petroleum, select=c("STATE","GHGRP_ID","CO2_GHGRP"))
colnames(GHGRP_CO2_chemicals_Other)[colnames(GHGRP_CO2_chemicals_Other) == "GHGRP.ID"] ="GHGRP_ID"
colnames(GHGRP_CO2_chemicals_Other)[colnames(GHGRP_CO2_chemicals_Other) == "GHG.QUANTITY..METRIC.TONS.CO2e."] ="CO2_GHGRP" #unit: metric tons/year
GHGRP_CO2_chemicals_Other = subset(GHGRP_CO2_chemicals_Other, select=c("STATE","GHGRP_ID","CO2_GHGRP"))

colnames(GHGRP_CO2_minerals_metals_Coal)[colnames(GHGRP_CO2_minerals_metals_Coal) == "GHGRP.ID"] ="GHGRP_ID"
colnames(GHGRP_CO2_minerals_metals_Coal)[colnames(GHGRP_CO2_minerals_metals_Coal) == "GHG.QUANTITY..METRIC.TONS.CO2e."] ="CO2_GHGRP" #unit: metric tons/year
GHGRP_CO2_minerals_metals_Coal = subset(GHGRP_CO2_minerals_metals_Coal, select=c("STATE","GHGRP_ID","CO2_GHGRP"))
colnames(GHGRP_CO2_minerals_metals_NG)[colnames(GHGRP_CO2_minerals_metals_NG) == "GHGRP.ID"] ="GHGRP_ID"
colnames(GHGRP_CO2_minerals_metals_NG)[colnames(GHGRP_CO2_minerals_metals_NG) == "GHG.QUANTITY..METRIC.TONS.CO2e."] ="CO2_GHGRP" #unit: metric tons/year
GHGRP_CO2_minerals_metals_NG = subset(GHGRP_CO2_minerals_metals_NG, select=c("STATE","GHGRP_ID","CO2_GHGRP"))
colnames(GHGRP_CO2_minerals_metals_Petroleum)[colnames(GHGRP_CO2_minerals_metals_Petroleum) == "GHGRP.ID"] ="GHGRP_ID"
colnames(GHGRP_CO2_minerals_metals_Petroleum)[colnames(GHGRP_CO2_minerals_metals_Petroleum) == "GHG.QUANTITY..METRIC.TONS.CO2e."] ="CO2_GHGRP" #unit: metric tons/year
GHGRP_CO2_minerals_metals_Petroleum = subset(GHGRP_CO2_minerals_metals_Petroleum, select=c("STATE","GHGRP_ID","CO2_GHGRP"))
colnames(GHGRP_CO2_minerals_metals_Other)[colnames(GHGRP_CO2_minerals_metals_Other) == "GHGRP.ID"] ="GHGRP_ID"
colnames(GHGRP_CO2_minerals_metals_Other)[colnames(GHGRP_CO2_minerals_metals_Other) == "GHG.QUANTITY..METRIC.TONS.CO2e."] ="CO2_GHGRP" #unit: metric tons/year
GHGRP_CO2_minerals_metals_Other = subset(GHGRP_CO2_minerals_metals_Other, select=c("STATE","GHGRP_ID","CO2_GHGRP"))

#############################################################################################################################
MSN_co2_fct_RESCOMM = MSN_co2_fct %>% filter(MSN %in% c("NGRCB", "DFRCB", "HLRCB", "KSRCB", "CLRCB", "NGCCB", "DFCCB", "HLCCB", "KSCCB", "RFCCB", "CLCCB", "MGCCB", "PCCCB"))
MSN_co2_fct_IND = MSN_co2_fct %>% filter(!(MSN %in% c("NGRCB", "DFRCB", "HLRCB", "KSRCB", "CLRCB", "NGCCB", "DFCCB", "HLCCB", "KSCCB", "RFCCB", "CLCCB", "MGCCB", "PCCCB")))

#join consump_states with MSN_co2_fct by MSN
calc_CO2_emis_RESCOMM = join(MSN_co2_fct_RESCOMM,btu_consump_states, by="MSN", type="left")
calc_CO2_emis_IND = join(MSN_co2_fct_IND,btu_consump_states, by="MSN", type="left")

#calculate annual total CO2 emissions in metric tons
calc_CO2_emis_RESCOMM['CO2_emis_2017'] <- calc_CO2_emis_RESCOMM['CO2_emis_fct'] * calc_CO2_emis_RESCOMM['adj_2017']
calc_CO2_emis_IND['CO2_emis_2017'] <- calc_CO2_emis_IND['CO2_emis_fct'] * calc_CO2_emis_IND['adj_2017']

#############################################################################################################################
#remove refineries, chemicals, minerals and metals GHGRP CO2 emissions from state-level

#subtract fuel-specific GHGRP at state-level
#grouping IND MSN into four categories: Coal, NG, Petroleum, and Other referencing CO2_coeffs_detailed spreadsheet from EPA
#these four groups are the same as GHGRP four fuel types
#subtract GHGRP four fuel types from these four groups of MSN CO2 emis totals at state-level
#after the GHGRP subtraction, remapping MSN to myCODE

################################################################################################################################################################
#for each state (including AK, HI) 50 states + DC + US = 52 regions, for each EF_group: Coal, NG, Petroleum, Other, sum over MSN
States <- (calc_CO2_emis_IND%>% distinct(State, .keep_all=TRUE))$State
num_states = as.numeric(length(States))

################################################################################
calc_CO2_emis_ind_States_sumMSN_Coal <- data.frame(Sector=character(num_states), 
                                                   State=character(num_states), 
                                                   EF_group=character(num_states), 
                                                   CO2_emis_2017=character(num_states)) 

count=1
for (i in States){
  calc_CO2_emis_ind_i = calc_CO2_emis_IND %>% filter(State == i & EF_group == 'Coal')
  Sector_i <- c("industrial")
  State_i <- c(i)
  EF_group_i <- c("Coal")
  CO2_emis_2017_i <- sum(calc_CO2_emis_ind_i$CO2_emis_2017)
  calc_CO2_emis_ind_States_sumMSN_Coal$Sector[count] = Sector_i
  calc_CO2_emis_ind_States_sumMSN_Coal$State[count] = State_i
  calc_CO2_emis_ind_States_sumMSN_Coal$EF_group[count] = EF_group_i
  calc_CO2_emis_ind_States_sumMSN_Coal$CO2_emis_2017[count] = CO2_emis_2017_i
  count = count + 1
}

################################################################################
calc_CO2_emis_ind_States_sumMSN_NG <- data.frame(Sector=character(num_states), 
                                                 State=character(num_states), 
                                                 EF_group=character(num_states), 
                                                 CO2_emis_2017=character(num_states)) 

count=1
for (i in States){
  calc_CO2_emis_ind_i = calc_CO2_emis_IND %>% filter(State == i & EF_group == 'NG')
  Sector_i <- c("industrial")
  State_i <- c(i)
  EF_group_i <- c("NG")
  CO2_emis_2017_i <- sum(calc_CO2_emis_ind_i$CO2_emis_2017)
  calc_CO2_emis_ind_States_sumMSN_NG$Sector[count] = Sector_i
  calc_CO2_emis_ind_States_sumMSN_NG$State[count] = State_i
  calc_CO2_emis_ind_States_sumMSN_NG$EF_group[count] = EF_group_i
  calc_CO2_emis_ind_States_sumMSN_NG$CO2_emis_2017[count] = CO2_emis_2017_i
  count = count + 1
}

################################################################################
calc_CO2_emis_ind_States_sumMSN_Petroleum <- data.frame(Sector=character(num_states), 
                                                        State=character(num_states), 
                                                        EF_group=character(num_states), 
                                                        CO2_emis_2017=character(num_states)) 

count=1
for (i in States){
  calc_CO2_emis_ind_i = calc_CO2_emis_IND %>% filter(State == i & EF_group == 'Petroleum')
  Sector_i <- c("industrial")
  State_i <- c(i)
  EF_group_i <- c("Petroleum")
  CO2_emis_2017_i <- sum(calc_CO2_emis_ind_i$CO2_emis_2017)
  calc_CO2_emis_ind_States_sumMSN_Petroleum$Sector[count] = Sector_i
  calc_CO2_emis_ind_States_sumMSN_Petroleum$State[count] = State_i
  calc_CO2_emis_ind_States_sumMSN_Petroleum$EF_group[count] = EF_group_i
  calc_CO2_emis_ind_States_sumMSN_Petroleum$CO2_emis_2017[count] = CO2_emis_2017_i
  count = count + 1
}

################################################################################
calc_CO2_emis_ind_States_sumMSN_Other <- data.frame(Sector=character(num_states), 
                                                    State=character(num_states), 
                                                    EF_group=character(num_states), 
                                                    CO2_emis_2017=character(num_states)) 

count=1
for (i in States){
  calc_CO2_emis_ind_i = calc_CO2_emis_IND %>% filter(State == i & EF_group == 'Other')
  Sector_i <- c("industrial")
  State_i <- c(i)
  EF_group_i <- c("Other")
  CO2_emis_2017_i <- sum(calc_CO2_emis_ind_i$CO2_emis_2017)
  calc_CO2_emis_ind_States_sumMSN_Other$Sector[count] = Sector_i
  calc_CO2_emis_ind_States_sumMSN_Other$State[count] = State_i
  calc_CO2_emis_ind_States_sumMSN_Other$EF_group[count] = EF_group_i
  calc_CO2_emis_ind_States_sumMSN_Other$CO2_emis_2017[count] = CO2_emis_2017_i
  count = count + 1
}

################################################################################################################################################################
################################################################################
#Sum GHGRP_CO2_refineries_Coal CO2 emissions by state
GHGRP_CO2_refineries_Coal_States_sum <- data.frame(State=character(num_states), 
                                              CO2_GHGRP_refineries=character(num_states)) 

count=1
for (i in States){
  GHGRP_CO2_refineries_Coal_i = GHGRP_CO2_refineries_Coal %>% filter(STATE == i)
  State_i <- c(i)
  CO2_GHGRP_i <- sum(GHGRP_CO2_refineries_Coal_i$CO2_GHGRP)
  GHGRP_CO2_refineries_Coal_States_sum$State[count] = State_i
  GHGRP_CO2_refineries_Coal_States_sum$CO2_GHGRP_refineries[count] = CO2_GHGRP_i
  count = count + 1
}

################################################################################
#Sum GHGRP_CO2_refineries_NG CO2 emissions by state
GHGRP_CO2_refineries_NG_States_sum <- data.frame(State=character(num_states), 
                                              CO2_GHGRP_refineries=character(num_states)) 

count=1
for (i in States){
  GHGRP_CO2_refineries_NG_i = GHGRP_CO2_refineries_NG %>% filter(STATE == i)
  State_i <- c(i)
  CO2_GHGRP_i <- sum(GHGRP_CO2_refineries_NG_i$CO2_GHGRP)
  GHGRP_CO2_refineries_NG_States_sum$State[count] = State_i
  GHGRP_CO2_refineries_NG_States_sum$CO2_GHGRP_refineries[count] = CO2_GHGRP_i
  count = count + 1
}

################################################################################
#Sum GHGRP_CO2_refineries_Petroleum CO2 emissions by state
GHGRP_CO2_refineries_Petroleum_States_sum <- data.frame(State=character(num_states), 
                                              CO2_GHGRP_refineries=character(num_states)) 

count=1
for (i in States){
  GHGRP_CO2_refineries_Petroleum_i = GHGRP_CO2_refineries_Petroleum %>% filter(STATE == i)
  State_i <- c(i)
  CO2_GHGRP_i <- sum(GHGRP_CO2_refineries_Petroleum_i$CO2_GHGRP)
  GHGRP_CO2_refineries_Petroleum_States_sum$State[count] = State_i
  GHGRP_CO2_refineries_Petroleum_States_sum$CO2_GHGRP_refineries[count] = CO2_GHGRP_i
  count = count + 1
}

################################################################################
#Sum GHGRP_CO2_refineries_Other CO2 emissions by state
GHGRP_CO2_refineries_Other_States_sum <- data.frame(State=character(num_states), 
                                              CO2_GHGRP_refineries=character(num_states)) 

count=1
for (i in States){
  GHGRP_CO2_refineries_Other_i = GHGRP_CO2_refineries_Other %>% filter(STATE == i)
  State_i <- c(i)
  CO2_GHGRP_i <- sum(GHGRP_CO2_refineries_Other_i$CO2_GHGRP)
  GHGRP_CO2_refineries_Other_States_sum$State[count] = State_i
  GHGRP_CO2_refineries_Other_States_sum$CO2_GHGRP_refineries[count] = CO2_GHGRP_i
  count = count + 1
}

################################################################################################################################################################
################################################################################
#Sum GHGRP_CO2_chemicals_Coal CO2 emissions by state
GHGRP_CO2_chemicals_Coal_States_sum <- data.frame(State=character(num_states), 
                                              CO2_GHGRP_chemicals=character(num_states)) 

count=1
for (i in States){
  GHGRP_CO2_chemicals_Coal_i = GHGRP_CO2_chemicals_Coal %>% filter(STATE == i)
  State_i <- c(i)
  CO2_GHGRP_i <- sum(GHGRP_CO2_chemicals_Coal_i$CO2_GHGRP)
  GHGRP_CO2_chemicals_Coal_States_sum$State[count] = State_i
  GHGRP_CO2_chemicals_Coal_States_sum$CO2_GHGRP_chemicals[count] = CO2_GHGRP_i
  count = count + 1
}

################################################################################
#Sum GHGRP_CO2_chemicals_NG CO2 emissions by state
GHGRP_CO2_chemicals_NG_States_sum <- data.frame(State=character(num_states), 
                                              CO2_GHGRP_chemicals=character(num_states)) 

count=1
for (i in States){
  GHGRP_CO2_chemicals_NG_i = GHGRP_CO2_chemicals_NG %>% filter(STATE == i)
  State_i <- c(i)
  CO2_GHGRP_i <- sum(GHGRP_CO2_chemicals_NG_i$CO2_GHGRP)
  GHGRP_CO2_chemicals_NG_States_sum$State[count] = State_i
  GHGRP_CO2_chemicals_NG_States_sum$CO2_GHGRP_chemicals[count] = CO2_GHGRP_i
  count = count + 1
}

################################################################################
#Sum GHGRP_CO2_chemicals_Petroleum CO2 emissions by state
GHGRP_CO2_chemicals_Petroleum_States_sum <- data.frame(State=character(num_states), 
                                              CO2_GHGRP_chemicals=character(num_states)) 

count=1
for (i in States){
  GHGRP_CO2_chemicals_Petroleum_i = GHGRP_CO2_chemicals_Petroleum %>% filter(STATE == i)
  State_i <- c(i)
  CO2_GHGRP_i <- sum(GHGRP_CO2_chemicals_Petroleum_i$CO2_GHGRP)
  GHGRP_CO2_chemicals_Petroleum_States_sum$State[count] = State_i
  GHGRP_CO2_chemicals_Petroleum_States_sum$CO2_GHGRP_chemicals[count] = CO2_GHGRP_i
  count = count + 1
}

################################################################################
#Sum GHGRP_CO2_chemicals_Other CO2 emissions by state
GHGRP_CO2_chemicals_Other_States_sum <- data.frame(State=character(num_states), 
                                              CO2_GHGRP_chemicals=character(num_states)) 

count=1
for (i in States){
  GHGRP_CO2_chemicals_Other_i = GHGRP_CO2_chemicals_Other %>% filter(STATE == i)
  State_i <- c(i)
  CO2_GHGRP_i <- sum(GHGRP_CO2_chemicals_Other_i$CO2_GHGRP)
  GHGRP_CO2_chemicals_Other_States_sum$State[count] = State_i
  GHGRP_CO2_chemicals_Other_States_sum$CO2_GHGRP_chemicals[count] = CO2_GHGRP_i
  count = count + 1
}

################################################################################################################################################################
################################################################################
#Sum GHGRP_CO2_minerals_metals_Coal CO2 emissions by state
GHGRP_CO2_minerals_metals_Coal_States_sum <- data.frame(State=character(num_states), 
                                              CO2_GHGRP_minerals_metals=character(num_states)) 

count=1
for (i in States){
  GHGRP_CO2_minerals_metals_Coal_i = GHGRP_CO2_minerals_metals_Coal %>% filter(STATE == i)
  State_i <- c(i)
  CO2_GHGRP_i <- sum(GHGRP_CO2_minerals_metals_Coal_i$CO2_GHGRP)
  GHGRP_CO2_minerals_metals_Coal_States_sum$State[count] = State_i
  GHGRP_CO2_minerals_metals_Coal_States_sum$CO2_GHGRP_minerals_metals[count] = CO2_GHGRP_i
  count = count + 1
}

################################################################################
#Sum GHGRP_CO2_minerals_metals_NG CO2 emissions by state
GHGRP_CO2_minerals_metals_NG_States_sum <- data.frame(State=character(num_states), 
                                              CO2_GHGRP_minerals_metals=character(num_states)) 

count=1
for (i in States){
  GHGRP_CO2_minerals_metals_NG_i = GHGRP_CO2_minerals_metals_NG %>% filter(STATE == i)
  State_i <- c(i)
  CO2_GHGRP_i <- sum(GHGRP_CO2_minerals_metals_NG_i$CO2_GHGRP)
  GHGRP_CO2_minerals_metals_NG_States_sum$State[count] = State_i
  GHGRP_CO2_minerals_metals_NG_States_sum$CO2_GHGRP_minerals_metals[count] = CO2_GHGRP_i
  count = count + 1
}

################################################################################
#Sum GHGRP_CO2_minerals_metals_Petroleum CO2 emissions by state
GHGRP_CO2_minerals_metals_Petroleum_States_sum <- data.frame(State=character(num_states), 
                                              CO2_GHGRP_minerals_metals=character(num_states)) 

count=1
for (i in States){
  GHGRP_CO2_minerals_metals_Petroleum_i = GHGRP_CO2_minerals_metals_Petroleum %>% filter(STATE == i)
  State_i <- c(i)
  CO2_GHGRP_i <- sum(GHGRP_CO2_minerals_metals_Petroleum_i$CO2_GHGRP)
  GHGRP_CO2_minerals_metals_Petroleum_States_sum$State[count] = State_i
  GHGRP_CO2_minerals_metals_Petroleum_States_sum$CO2_GHGRP_minerals_metals[count] = CO2_GHGRP_i
  count = count + 1
}

################################################################################
#Sum GHGRP_CO2_minerals_metals_Other CO2 emissions by state
GHGRP_CO2_minerals_metals_Other_States_sum <- data.frame(State=character(num_states), 
                                              CO2_GHGRP_minerals_metals=character(num_states)) 

count=1
for (i in States){
  GHGRP_CO2_minerals_metals_Other_i = GHGRP_CO2_minerals_metals_Other %>% filter(STATE == i)
  State_i <- c(i)
  CO2_GHGRP_i <- sum(GHGRP_CO2_minerals_metals_Other_i$CO2_GHGRP)
  GHGRP_CO2_minerals_metals_Other_States_sum$State[count] = State_i
  GHGRP_CO2_minerals_metals_Other_States_sum$CO2_GHGRP_minerals_metals[count] = CO2_GHGRP_i
  count = count + 1
}

################################################################################################################################################################
################################################################################
#join calc_CO2_emis_ind_States_sumMSN_Coal, GHGRP_CO2_refineries_Coal_States_sum, GHGRP_CO2_chemicals_Coal_States_sum, GHGRP_CO2_minerals_metals_Coal_States_sum and  by state
calc_CO2_emis_ind_States_sumMSN_Coal_rmGHGRP = join(calc_CO2_emis_ind_States_sumMSN_Coal,GHGRP_CO2_refineries_Coal_States_sum,
                                     by=c("State"), type="left")

calc_CO2_emis_ind_States_sumMSN_Coal_rmGHGRP = join(calc_CO2_emis_ind_States_sumMSN_Coal_rmGHGRP,GHGRP_CO2_chemicals_Coal_States_sum,
                                               by=c("State"), type="left")

calc_CO2_emis_ind_States_sumMSN_Coal_rmGHGRP = join(calc_CO2_emis_ind_States_sumMSN_Coal_rmGHGRP,GHGRP_CO2_minerals_metals_Coal_States_sum,
                                               by=c("State"), type="left")

calc_CO2_emis_ind_States_sumMSN_Coal_rmGHGRP$CO2_nonGHGRP <- as.numeric(calc_CO2_emis_ind_States_sumMSN_Coal_rmGHGRP$CO2_emis_2017) - 
                                                        as.numeric(calc_CO2_emis_ind_States_sumMSN_Coal_rmGHGRP$CO2_GHGRP_refineries) -
                                                        as.numeric(calc_CO2_emis_ind_States_sumMSN_Coal_rmGHGRP$CO2_GHGRP_chemicals) -
                                                        as.numeric(calc_CO2_emis_ind_States_sumMSN_Coal_rmGHGRP$CO2_GHGRP_minerals_metals)

calc_CO2_emis_ind_States_sumMSN_Coal_rmGHGRP$CO2_nonGHGRP_frac <- as.numeric(calc_CO2_emis_ind_States_sumMSN_Coal_rmGHGRP$CO2_nonGHGRP)/
                                                             as.numeric(calc_CO2_emis_ind_States_sumMSN_Coal_rmGHGRP$CO2_emis_2017)

write.csv(calc_CO2_emis_ind_States_sumMSN_Coal_rmGHGRP, file = "calc_CO2_emis_ind_States_sumMSN_Coal_rmGHGRP_v12.csv")

#If for a state the CO2_nonGHGRP_frac is negative, make it zero
calc_CO2_emis_ind_States_sumMSN_Coal_rmGHGRP$CO2_nonGHGRP_frac[calc_CO2_emis_ind_States_sumMSN_Coal_rmGHGRP$CO2_nonGHGRP_frac<0] = 0

#If for a state the CO2_nonGHGRP_frac is NaN due to 0/0, make it 1
calc_CO2_emis_ind_States_sumMSN_Coal_rmGHGRP$CO2_nonGHGRP_frac[is.nan(calc_CO2_emis_ind_States_sumMSN_Coal_rmGHGRP$CO2_nonGHGRP_frac) &
                                                               calc_CO2_emis_ind_States_sumMSN_Coal_rmGHGRP$CO2_emis_2017==0 & calc_CO2_emis_ind_States_sumMSN_Coal_rmGHGRP$CO2_nonGHGRP==0] = 1

state_CO2_nonGHGRP_frac_Coal = subset(calc_CO2_emis_ind_States_sumMSN_Coal_rmGHGRP, select=c("EF_group", "State", "CO2_nonGHGRP_frac" ))

################################################################################
#join calc_CO2_emis_ind_States_sumMSN_NG, GHGRP_CO2_refineries_NG_States_sum, GHGRP_CO2_chemicals_NG_States_sum, GHGRP_CO2_minerals_metals_NG_States_sum and  by state
calc_CO2_emis_ind_States_sumMSN_NG_rmGHGRP = join(calc_CO2_emis_ind_States_sumMSN_NG,GHGRP_CO2_refineries_NG_States_sum,
                                     by=c("State"), type="left")

calc_CO2_emis_ind_States_sumMSN_NG_rmGHGRP = join(calc_CO2_emis_ind_States_sumMSN_NG_rmGHGRP,GHGRP_CO2_chemicals_NG_States_sum,
                                               by=c("State"), type="left")

calc_CO2_emis_ind_States_sumMSN_NG_rmGHGRP = join(calc_CO2_emis_ind_States_sumMSN_NG_rmGHGRP,GHGRP_CO2_minerals_metals_NG_States_sum,
                                               by=c("State"), type="left")

calc_CO2_emis_ind_States_sumMSN_NG_rmGHGRP$CO2_nonGHGRP <- as.numeric(calc_CO2_emis_ind_States_sumMSN_NG_rmGHGRP$CO2_emis_2017) - 
                                                        as.numeric(calc_CO2_emis_ind_States_sumMSN_NG_rmGHGRP$CO2_GHGRP_refineries) -
                                                        as.numeric(calc_CO2_emis_ind_States_sumMSN_NG_rmGHGRP$CO2_GHGRP_chemicals) -
                                                        as.numeric(calc_CO2_emis_ind_States_sumMSN_NG_rmGHGRP$CO2_GHGRP_minerals_metals)

calc_CO2_emis_ind_States_sumMSN_NG_rmGHGRP$CO2_nonGHGRP_frac <- as.numeric(calc_CO2_emis_ind_States_sumMSN_NG_rmGHGRP$CO2_nonGHGRP)/
                                                             as.numeric(calc_CO2_emis_ind_States_sumMSN_NG_rmGHGRP$CO2_emis_2017)

write.csv(calc_CO2_emis_ind_States_sumMSN_NG_rmGHGRP, file = "calc_CO2_emis_ind_States_sumMSN_NG_rmGHGRP_v12.csv")

#If for a state the CO2_nonGHGRP_frac is negative, make it zero
calc_CO2_emis_ind_States_sumMSN_NG_rmGHGRP$CO2_nonGHGRP_frac[calc_CO2_emis_ind_States_sumMSN_NG_rmGHGRP$CO2_nonGHGRP_frac<0] = 0

#If for a state the CO2_nonGHGRP_frac is NaN due to 0/0, make it 1
calc_CO2_emis_ind_States_sumMSN_NG_rmGHGRP$CO2_nonGHGRP_frac[is.nan(calc_CO2_emis_ind_States_sumMSN_NG_rmGHGRP$CO2_nonGHGRP_frac) &
                                                               calc_CO2_emis_ind_States_sumMSN_NG_rmGHGRP$CO2_emis_2017==0 & calc_CO2_emis_ind_States_sumMSN_NG_rmGHGRP$CO2_nonGHGRP==0] = 1

state_CO2_nonGHGRP_frac_NG = subset(calc_CO2_emis_ind_States_sumMSN_NG_rmGHGRP, select=c("EF_group", "State", "CO2_nonGHGRP_frac" ))

################################################################################
#join calc_CO2_emis_ind_States_sumMSN_Petroleum, GHGRP_CO2_refineries_Petroleum_States_sum, GHGRP_CO2_chemicals_Petroleum_States_sum, GHGRP_CO2_minerals_metals_Petroleum_States_sum and  by state
calc_CO2_emis_ind_States_sumMSN_Petroleum_rmGHGRP = join(calc_CO2_emis_ind_States_sumMSN_Petroleum,GHGRP_CO2_refineries_Petroleum_States_sum,
                                     by=c("State"), type="left")

calc_CO2_emis_ind_States_sumMSN_Petroleum_rmGHGRP = join(calc_CO2_emis_ind_States_sumMSN_Petroleum_rmGHGRP,GHGRP_CO2_chemicals_Petroleum_States_sum,
                                               by=c("State"), type="left")

calc_CO2_emis_ind_States_sumMSN_Petroleum_rmGHGRP = join(calc_CO2_emis_ind_States_sumMSN_Petroleum_rmGHGRP,GHGRP_CO2_minerals_metals_Petroleum_States_sum,
                                               by=c("State"), type="left")

calc_CO2_emis_ind_States_sumMSN_Petroleum_rmGHGRP$CO2_nonGHGRP <- as.numeric(calc_CO2_emis_ind_States_sumMSN_Petroleum_rmGHGRP$CO2_emis_2017) - 
                                                        as.numeric(calc_CO2_emis_ind_States_sumMSN_Petroleum_rmGHGRP$CO2_GHGRP_refineries) -
                                                        as.numeric(calc_CO2_emis_ind_States_sumMSN_Petroleum_rmGHGRP$CO2_GHGRP_chemicals) -
                                                        as.numeric(calc_CO2_emis_ind_States_sumMSN_Petroleum_rmGHGRP$CO2_GHGRP_minerals_metals)

calc_CO2_emis_ind_States_sumMSN_Petroleum_rmGHGRP$CO2_nonGHGRP_frac <- as.numeric(calc_CO2_emis_ind_States_sumMSN_Petroleum_rmGHGRP$CO2_nonGHGRP)/
                                                             as.numeric(calc_CO2_emis_ind_States_sumMSN_Petroleum_rmGHGRP$CO2_emis_2017)

write.csv(calc_CO2_emis_ind_States_sumMSN_Petroleum_rmGHGRP, file = "calc_CO2_emis_ind_States_sumMSN_Petroleum_rmGHGRP_v12.csv")

#If for a state the CO2_nonGHGRP_frac is negative, make it zero
calc_CO2_emis_ind_States_sumMSN_Petroleum_rmGHGRP$CO2_nonGHGRP_frac[calc_CO2_emis_ind_States_sumMSN_Petroleum_rmGHGRP$CO2_nonGHGRP_frac<0] = 0

#If for a state the CO2_nonGHGRP_frac is NaN due to 0/0, make it 1
calc_CO2_emis_ind_States_sumMSN_Petroleum_rmGHGRP$CO2_nonGHGRP_frac[is.nan(calc_CO2_emis_ind_States_sumMSN_Petroleum_rmGHGRP$CO2_nonGHGRP_frac) &
                                                               calc_CO2_emis_ind_States_sumMSN_Petroleum_rmGHGRP$CO2_emis_2017==0 & calc_CO2_emis_ind_States_sumMSN_Petroleum_rmGHGRP$CO2_nonGHGRP==0] = 1

state_CO2_nonGHGRP_frac_Petroleum = subset(calc_CO2_emis_ind_States_sumMSN_Petroleum_rmGHGRP, select=c("EF_group", "State", "CO2_nonGHGRP_frac" ))

################################################################################
#join calc_CO2_emis_ind_States_sumMSN_Other, GHGRP_CO2_refineries_Other_States_sum, GHGRP_CO2_chemicals_Other_States_sum, GHGRP_CO2_minerals_metals_Other_States_sum and  by state
calc_CO2_emis_ind_States_sumMSN_Other_rmGHGRP = join(calc_CO2_emis_ind_States_sumMSN_Other,GHGRP_CO2_refineries_Other_States_sum,
                                     by=c("State"), type="left")

calc_CO2_emis_ind_States_sumMSN_Other_rmGHGRP = join(calc_CO2_emis_ind_States_sumMSN_Other_rmGHGRP,GHGRP_CO2_chemicals_Other_States_sum,
                                               by=c("State"), type="left")

calc_CO2_emis_ind_States_sumMSN_Other_rmGHGRP = join(calc_CO2_emis_ind_States_sumMSN_Other_rmGHGRP,GHGRP_CO2_minerals_metals_Other_States_sum,
                                               by=c("State"), type="left")

calc_CO2_emis_ind_States_sumMSN_Other_rmGHGRP$CO2_nonGHGRP <- as.numeric(calc_CO2_emis_ind_States_sumMSN_Other_rmGHGRP$CO2_emis_2017) - 
                                                        as.numeric(calc_CO2_emis_ind_States_sumMSN_Other_rmGHGRP$CO2_GHGRP_refineries) -
                                                        as.numeric(calc_CO2_emis_ind_States_sumMSN_Other_rmGHGRP$CO2_GHGRP_chemicals) -
                                                        as.numeric(calc_CO2_emis_ind_States_sumMSN_Other_rmGHGRP$CO2_GHGRP_minerals_metals)

calc_CO2_emis_ind_States_sumMSN_Other_rmGHGRP$CO2_nonGHGRP_frac <- as.numeric(calc_CO2_emis_ind_States_sumMSN_Other_rmGHGRP$CO2_nonGHGRP)/
                                                             as.numeric(calc_CO2_emis_ind_States_sumMSN_Other_rmGHGRP$CO2_emis_2017)

write.csv(calc_CO2_emis_ind_States_sumMSN_Other_rmGHGRP, file = "calc_CO2_emis_ind_States_sumMSN_Other_rmGHGRP_v12.csv")

#If for a state the CO2_nonGHGRP_frac is negative, make it zero
calc_CO2_emis_ind_States_sumMSN_Other_rmGHGRP$CO2_nonGHGRP_frac[calc_CO2_emis_ind_States_sumMSN_Other_rmGHGRP$CO2_nonGHGRP_frac<0] = 0

#If for a state the CO2_nonGHGRP_frac is NaN due to 0/0, make it 1
calc_CO2_emis_ind_States_sumMSN_Other_rmGHGRP$CO2_nonGHGRP_frac[is.nan(calc_CO2_emis_ind_States_sumMSN_Other_rmGHGRP$CO2_nonGHGRP_frac) &
                                                               calc_CO2_emis_ind_States_sumMSN_Other_rmGHGRP$CO2_emis_2017==0 & calc_CO2_emis_ind_States_sumMSN_Other_rmGHGRP$CO2_nonGHGRP==0] = 1

state_CO2_nonGHGRP_frac_Other = subset(calc_CO2_emis_ind_States_sumMSN_Other_rmGHGRP, select=c("EF_group", "State", "CO2_nonGHGRP_frac" ))

################################################################################
#stack state_CO2_nonGHGRP_frac_*
state_CO2_nonGHGRP_frac_all = rbind(state_CO2_nonGHGRP_frac_Coal,state_CO2_nonGHGRP_frac_NG,state_CO2_nonGHGRP_frac_Petroleum,state_CO2_nonGHGRP_frac_Other)

#apply state_CO2_nonGHGRP_frac_all to calc_CO2_emis_IND['CO2_emis_2017']
calc_CO2_emis_IND = join(calc_CO2_emis_IND,state_CO2_nonGHGRP_frac_all,by=c("EF_group", "State"), type="left")

calc_CO2_emis_IND$CO2_emis_2017_nonGHGRP = calc_CO2_emis_IND$CO2_emis_2017 * calc_CO2_emis_IND$CO2_nonGHGRP_frac

#subset
EIA_CO2_emis_ind = subset(calc_CO2_emis_IND, select=c("MSN","Description", "State", "CO2_emis_2017_nonGHGRP" ))
colnames(EIA_CO2_emis_ind)[colnames(EIA_CO2_emis_ind) == "CO2_emis_2017_nonGHGRP"] ="CO2_emis_2017"

EIA_CO2_emis_rescomm = subset(calc_CO2_emis_RESCOMM, select=c("MSN","Description", "State", "CO2_emis_2017"))

EIA_CO2_emis <- rbind(EIA_CO2_emis_rescomm, EIA_CO2_emis_ind)

###############################################################################################################################
#divide EIA_CO2_emis into sectors
#residential NGRCB DFRCB HLRCB KSRCB CLRCB
EIA_CO2_emis_res1 = EIA_CO2_emis %>% filter(MSN == "NGRCB")
EIA_CO2_emis_res2 = EIA_CO2_emis %>% filter(MSN == "DFRCB")
EIA_CO2_emis_res3 = EIA_CO2_emis %>% filter(MSN == "HLRCB")
EIA_CO2_emis_res4 = EIA_CO2_emis %>% filter(MSN == "KSRCB")
EIA_CO2_emis_res5 = EIA_CO2_emis %>% filter(MSN == "CLRCB")

EIA_CO2_emis_res <- rbind(EIA_CO2_emis_res1, EIA_CO2_emis_res2, 
                          EIA_CO2_emis_res3, EIA_CO2_emis_res4, 
                          EIA_CO2_emis_res5)

#commercial NGCCB DFCCB HLCCB KSCCB RFCCB CLCCB MGCCB PCCCB
EIA_CO2_emis_com1 = EIA_CO2_emis %>% filter(MSN == "NGCCB")
EIA_CO2_emis_com2 = EIA_CO2_emis %>% filter(MSN == "DFCCB")
EIA_CO2_emis_com3 = EIA_CO2_emis %>% filter(MSN == "HLCCB")
EIA_CO2_emis_com4 = EIA_CO2_emis %>% filter(MSN == "KSCCB")
EIA_CO2_emis_com5 = EIA_CO2_emis %>% filter(MSN == "RFCCB")
EIA_CO2_emis_com6 = EIA_CO2_emis %>% filter(MSN == "CLCCB")
EIA_CO2_emis_com7 = EIA_CO2_emis %>% filter(MSN == "MGCCB")
EIA_CO2_emis_com8 = EIA_CO2_emis %>% filter(MSN == "PCCCB")

EIA_CO2_emis_com <- rbind(EIA_CO2_emis_com1, EIA_CO2_emis_com2, 
                          EIA_CO2_emis_com3, EIA_CO2_emis_com4, 
                          EIA_CO2_emis_com5, EIA_CO2_emis_com6,
                          EIA_CO2_emis_com7, EIA_CO2_emis_com8)

###############################################################################################################################################
#for each sector, and for each state (including AK, HI) 50 states + DC + US = 52 regions, sum over all the MSN categories
States <- EIA_CO2_emis_res1$State
num_states = as.numeric(length(States))

#residential
EIA_CO2_emis_res_States_sumMSN <- data.frame(Sector=character(num_states), 
                                             State=character(num_states), 
                                             CO2_emis_2017=character(num_states)) 

count=1
for (i in States){
  EIA_CO2_emis_res_i = EIA_CO2_emis_res %>% filter(State == i)
  Sector_i <- c("Residential")
  State_i <- c(i)
  CO2_emis_2017_i <- sum(EIA_CO2_emis_res_i$CO2_emis_2017)
  EIA_CO2_emis_res_States_sumMSN$Sector[count] = Sector_i
  EIA_CO2_emis_res_States_sumMSN$State[count] = State_i
  EIA_CO2_emis_res_States_sumMSN$CO2_emis_2017[count] = CO2_emis_2017_i
  count = count + 1
}

#commercial
EIA_CO2_emis_com_States_sumMSN <- data.frame(Sector=character(num_states), 
                                             State=character(num_states), 
                                             CO2_emis_2017=character(num_states)) 

count=1
for (i in States){
  EIA_CO2_emis_com_i = EIA_CO2_emis_com %>% filter(State == i)
  Sector_i <- c("commerical")
  State_i <- c(i)
  CO2_emis_2017_i <- sum(EIA_CO2_emis_com_i$CO2_emis_2017)
  EIA_CO2_emis_com_States_sumMSN$Sector[count] = Sector_i
  EIA_CO2_emis_com_States_sumMSN$State[count] = State_i
  EIA_CO2_emis_com_States_sumMSN$CO2_emis_2017[count] = CO2_emis_2017_i
  count = count + 1
}

#industrial
EIA_CO2_emis_ind_States_sumMSN <- data.frame(Sector=character(num_states), 
                                             State=character(num_states), 
                                             CO2_emis_2017=character(num_states)) 

count=1
for (i in States){
  EIA_CO2_emis_ind_i = EIA_CO2_emis_ind %>% filter(State == i)
  Sector_i <- c("industrial")
  State_i <- c(i)
  CO2_emis_2017_i <- sum(EIA_CO2_emis_ind_i$CO2_emis_2017)
  EIA_CO2_emis_ind_States_sumMSN$Sector[count] = Sector_i
  EIA_CO2_emis_ind_States_sumMSN$State[count] = State_i
  EIA_CO2_emis_ind_States_sumMSN$CO2_emis_2017[count] = CO2_emis_2017_i
  count = count + 1
}

#########################################################################################################################
#Divide annual total emissions of CO2 (metric tons) into weekday, Saturday, and Sunday before applying point/area ratio
#preparations
# 1. rename "sector" in AREA_sectors_CO and POINT_sectors_CO to AREAmyCODE and POINTmyCODE respectively
colnames(AREA_sectors_CO)[colnames(AREA_sectors_CO) == "sector"] ="AREAmyCODE"
colnames(POINT_sectors_CO)[colnames(POINT_sectors_CO) == "sector"] ="POINTmyCODE"
colnames(AREA_sectors_CO_patch)[colnames(AREA_sectors_CO_patch) == "sector"] ="AREAmyCODE"
colnames(POINT_sectors_CO_patch)[colnames(POINT_sectors_CO_patch) == "sector"] ="POINTmyCODE"

#2.0 get point myCODEs real
point_myCODEs_real <- (POINT_sectors_CO%>% distinct(POINTmyCODE, .keep_all=TRUE))$POINTmyCODE
num_point_myCODEs_real = as.numeric(length(point_myCODEs_real))

# 2.1 combine PtIND_NG and PtIND_NG2 in POINT_sectors_CO as PtIND_NG
POINT_sectors_CO_IND_NG = POINT_sectors_CO%>% filter(POINTmyCODE == 'PtIND_NG')
POINT_sectors_CO_IND_NG2 = POINT_sectors_CO%>% filter(POINTmyCODE == 'PtIND_NG2')

POINT_sectors_CO_IND_NG2 = subset(POINT_sectors_CO_IND_NG2, 
                                  select=c("day","STATE_FIPS","STATE_NAME","emis_sum"))

colnames(POINT_sectors_CO_IND_NG2)[colnames(POINT_sectors_CO_IND_NG2) == "emis_sum"] ="emis_sum2"

POINT_sectors_CO_IND_NG_both = join(POINT_sectors_CO_IND_NG,POINT_sectors_CO_IND_NG2, 
                                    by=c("day","STATE_FIPS","STATE_NAME"), type="left")

POINT_sectors_CO_IND_NG_both['emis_sum_both'] <- as.numeric(POINT_sectors_CO_IND_NG_both$emis_sum) + 
                                                 as.numeric(POINT_sectors_CO_IND_NG_both$emis_sum2)

POINT_sectors_CO_IND_NG_both = subset(POINT_sectors_CO_IND_NG_both, 
                                      select=c("POINTmyCODE","day","STATE_FIPS","STATE_NAME","emis_sum_both"))

#remove PtIND_NG2 rows from POINT_sectors_CO
POINT_sectors_CO = subset(POINT_sectors_CO, POINTmyCODE != "PtIND_NG2")

#put the combined PtIND_NG back
POINT_sectors_CO = join(POINT_sectors_CO,POINT_sectors_CO_IND_NG_both, 
                        by=c("POINTmyCODE","day","STATE_FIPS","STATE_NAME"), type="left")
POINT_sectors_CO[is.na(POINT_sectors_CO)] = 0
POINT_sectors_CO['emis_sum_larger'] <- pmax(POINT_sectors_CO$emis_sum, POINT_sectors_CO$emis_sum_both)
POINT_sectors_CO = subset(POINT_sectors_CO, select=c("POINTmyCODE","day","STATE_FIPS","STATE_NAME","emis_sum_larger"))
colnames(POINT_sectors_CO)[colnames(POINT_sectors_CO) == "emis_sum_larger"] ="emis_sum"

# 2.2 combine PtIND_Oil and PtIND_Oil2 in POINT_sectors_CO as PtIND_Oil
POINT_sectors_CO_IND_Oil = POINT_sectors_CO%>% filter(POINTmyCODE == 'PtIND_Oil')
POINT_sectors_CO_IND_Oil2 = POINT_sectors_CO%>% filter(POINTmyCODE == 'PtIND_Oil2')

POINT_sectors_CO_IND_Oil2 = subset(POINT_sectors_CO_IND_Oil2, select=c("day","STATE_FIPS","STATE_NAME","emis_sum"))

colnames(POINT_sectors_CO_IND_Oil2)[colnames(POINT_sectors_CO_IND_Oil2) == "emis_sum"] ="emis_sum2"

POINT_sectors_CO_IND_Oil_both = join(POINT_sectors_CO_IND_Oil,POINT_sectors_CO_IND_Oil2, 
                                     by=c("day","STATE_FIPS","STATE_NAME"), type="left")

POINT_sectors_CO_IND_Oil_both['emis_sum_both'] <- as.numeric(POINT_sectors_CO_IND_Oil_both$emis_sum) + 
                                                  as.numeric(POINT_sectors_CO_IND_Oil_both$emis_sum2)

POINT_sectors_CO_IND_Oil_both = subset(POINT_sectors_CO_IND_Oil_both, 
                                       select=c("POINTmyCODE","day","STATE_FIPS","STATE_NAME","emis_sum_both"))

#remove PtIND_Oil2 rows from POINT_sectors_CO
POINT_sectors_CO = subset(POINT_sectors_CO, POINTmyCODE != "PtIND_Oil2")

#put the combined PtIND_Oil back
POINT_sectors_CO = join(POINT_sectors_CO,POINT_sectors_CO_IND_Oil_both, 
                        by=c("POINTmyCODE","day","STATE_FIPS","STATE_NAME"), type="left")
POINT_sectors_CO[is.na(POINT_sectors_CO)] = 0
POINT_sectors_CO['emis_sum_larger'] <- pmax(POINT_sectors_CO$emis_sum, POINT_sectors_CO$emis_sum_both)
POINT_sectors_CO = subset(POINT_sectors_CO, select=c("POINTmyCODE","day","STATE_FIPS","STATE_NAME","emis_sum_larger"))
colnames(POINT_sectors_CO)[colnames(POINT_sectors_CO) == "emis_sum_larger"] ="emis_sum"

colnames(POINT_sectors_CO)[colnames(POINT_sectors_CO) == "emis_sum"] ="emis_sum_POINT"
colnames(POINT_sectors_CO)[colnames(POINT_sectors_CO) == "STATE_FIPS"] ="STATE_FIPS_POINT"
colnames(AREA_sectors_CO)[colnames(AREA_sectors_CO) == "emis_sum"] ="emis_sum_AREA"
colnames(AREA_sectors_CO)[colnames(AREA_sectors_CO) == "STATE_FIPS"] ="STATE_FIPS_AREA"

POINT_sectors_CO_patch = subset(POINT_sectors_CO_patch, select=c("POINTmyCODE","day","STATE_FIPS","STATE_NAME","emis_sum"))
colnames(POINT_sectors_CO_patch)[colnames(POINT_sectors_CO_patch) == "emis_sum"] ="emis_sum_POINT"
colnames(POINT_sectors_CO_patch)[colnames(POINT_sectors_CO_patch) == "STATE_FIPS"] ="STATE_FIPS_POINT"
colnames(AREA_sectors_CO_patch)[colnames(AREA_sectors_CO_patch) == "emis_sum"] ="emis_sum_AREA"
colnames(AREA_sectors_CO_patch)[colnames(AREA_sectors_CO_patch) == "STATE_FIPS"] ="STATE_FIPS_AREA"

# 3. add AREAmyCODE to POINT_sectors_CO
POINT_sectors_CO = join(POINT_sectors_CO,POINTmyCODE_to_AREAmyCODE, by="POINTmyCODE", type="left")
POINT_sectors_CO_patch = join(POINT_sectors_CO_patch,POINTmyCODE_to_AREAmyCODE, by="POINTmyCODE", type="left")

# 4. join POINT_sectors_CO with AREA_sectors_CO to sum POINT and AREA emissions of CO
TOTAL_sectors_CO = join(AREA_sectors_CO, POINT_sectors_CO, by=c("AREAmyCODE","day","STATE_NAME"), type="left")

# 5. make NaN = 0
TOTAL_sectors_CO[is.na(TOTAL_sectors_CO)] = 0

# 6. calculate area+point emissions of CO
TOTAL_sectors_CO['emis_sum_TOTAL'] <- as.numeric(TOTAL_sectors_CO$emis_sum_AREA) + as.numeric(TOTAL_sectors_CO$emis_sum_POINT)

# 7. subset TOTAL_sectors_CO
TOTAL_sectors_CO = subset(TOTAL_sectors_CO, select=c("AREAmyCODE","day","STATE_FIPS_AREA","STATE_NAME","emis_sum_TOTAL"))

# 8. Add state abbreviation
POINT_sectors_CO_Stateabb = join(POINT_sectors_CO,State_fullname_to_abb, by="STATE_NAME", type="left")
AREA_sectors_CO_Stateabb = join(AREA_sectors_CO,State_fullname_to_abb, by="STATE_NAME", type="left")
TOTAL_sectors_CO_Stateabb = join(TOTAL_sectors_CO,State_fullname_to_abb, by="STATE_NAME", type="left")

POINT_sectors_CO_patch_Stateabb = join(POINT_sectors_CO_patch,State_fullname_to_abb, by="STATE_NAME", type="left")
AREA_sectors_CO_patch_Stateabb = join(AREA_sectors_CO_patch,State_fullname_to_abb, by="STATE_NAME", type="left")

# 9. get area myCODEs
AREA_sectors_CO_myCODEs <- AREA_sectors_CO_Stateabb%>% filter(State == 'AL' & day == 'weekdy')
area_myCODEs <- AREA_sectors_CO_myCODEs$AREAmyCODE
num_area_myCODEs = as.numeric(length(area_myCODEs))

# 10. get point myCODEs after IND_NG2 and IND_Oil2 are combined into IND_NG and IND_Oil
POINT_sectors_CO_myCODEs <- POINT_sectors_CO_Stateabb%>% filter(State == 'AL' & day == 'weekdy')
point_myCODEs <- POINT_sectors_CO_myCODEs$POINTmyCODE
num_point_myCODEs = as.numeric(length(point_myCODEs))

# 11.get CONUS states
POINT_sectors_CO_States <- POINT_sectors_CO_Stateabb%>% filter(POINTmyCODE == 'PtCOMM_NG' & day == 'weekdy')
States_CONUS <- POINT_sectors_CO_States$State
num_states_CONUS = as.numeric(length(States_CONUS))

TOTAL_CO_emis_days_frac <- data.frame(myCODE=character(num_states_CONUS*num_area_myCODEs), 
                                      State=character(num_states_CONUS*num_area_myCODEs), 
                                      weekdy_emis_frac=character(num_states_CONUS*num_area_myCODEs), 
                                      satdy_emis_frac=character(num_states_CONUS*num_area_myCODEs), 
                                      sundy_emis_frac=character(num_states_CONUS*num_area_myCODEs)) 
count = 1
for (j in area_myCODEs){
  TOTAL_sectors_CO_Stateabb_j = TOTAL_sectors_CO_Stateabb %>% filter(AREAmyCODE == j) #50*3
  myCODE_j <- c(j)
  for (i in States_CONUS){
    TOTAL_sectors_CO_Stateabb_ji = TOTAL_sectors_CO_Stateabb_j %>% filter(State == i) #3
    weekdy_emis = TOTAL_sectors_CO_Stateabb_ji%>% filter(day == 'weekdy')
    satdy_emis = TOTAL_sectors_CO_Stateabb_ji%>% filter(day == 'satdy')
    sundy_emis = TOTAL_sectors_CO_Stateabb_ji%>% filter(day == 'sundy')
    weekdy_emis = weekdy_emis$emis_sum_TOTAL #metric tons/day
    satdy_emis = satdy_emis$emis_sum_TOTAL
    sundy_emis = sundy_emis$emis_sum_TOTAL
    week_emis_sum = weekdy_emis*5 + satdy_emis*1 + sundy_emis*1 #metric tons/week
    weekdy_frac = (weekdy_emis*5)/week_emis_sum
    satdy_frac = (satdy_emis*1)/week_emis_sum
    sundy_frac = (sundy_emis*1)/week_emis_sum
    State_i <-c(i)
    #save calculated day frac
    TOTAL_CO_emis_days_frac$myCODE[count] = myCODE_j
    TOTAL_CO_emis_days_frac$State[count]= State_i
    TOTAL_CO_emis_days_frac$weekdy_emis_frac[count] = weekdy_frac
    TOTAL_CO_emis_days_frac$satdy_emis_frac[count] = satdy_frac
    TOTAL_CO_emis_days_frac$sundy_emis_frac[count] = sundy_frac
    count = count + 1
  }
}

#if *dy_frac in TOTAL_CO_emis_days_frac are NaNs (caused by week_emis_sum = 0)
#use average values of the same myCODE to replace NaNs
count = 1
for (j in area_myCODEs){
  TOTAL_CO_emis_days_frac_j = TOTAL_CO_emis_days_frac %>% filter(myCODE == j) #50
  #calculate the 50 states avg.
  weekdy_emis_frac_avg = mean(as.numeric(TOTAL_CO_emis_days_frac_j$weekdy_emis_frac), na.rm=TRUE)
  satdy_emis_frac_avg = mean(as.numeric(TOTAL_CO_emis_days_frac_j$satdy_emis_frac), na.rm=TRUE)
  sundy_emis_frac_avg = mean(as.numeric(TOTAL_CO_emis_days_frac_j$sundy_emis_frac), na.rm=TRUE)
  for (i in States_CONUS){
    TOTAL_CO_emis_days_frac_ji = TOTAL_CO_emis_days_frac_j %>% filter(State == i)
    if (is.nan(as.numeric(TOTAL_CO_emis_days_frac_ji$weekdy_emis_frac))){
      TOTAL_CO_emis_days_frac$weekdy_emis_frac[count] = weekdy_emis_frac_avg
      TOTAL_CO_emis_days_frac$satdy_emis_frac[count] = satdy_emis_frac_avg
      TOTAL_CO_emis_days_frac$sundy_emis_frac[count] = sundy_emis_frac_avg
    } 
    count = count + 1
  }
}

#do a sanity check of the fractions by sum them up
TOTAL_CO_emis_days_frac['emis_frac_tot'] <- as.numeric(TOTAL_CO_emis_days_frac$weekdy_emis_frac) + 
                                            as.numeric(TOTAL_CO_emis_days_frac$satdy_emis_frac) + 
                                            as.numeric(TOTAL_CO_emis_days_frac$sundy_emis_frac)

#get rid of AK and HI for all three sectors
EIA_CO2_emis_CONUS = subset(EIA_CO2_emis, State != c("AK", "HI"))
EIA_CO2_emis_res_CONUS = subset(EIA_CO2_emis_res, State != c("AK", "HI"))
EIA_CO2_emis_com_CONUS = subset(EIA_CO2_emis_com, State != c("AK", "HI"))
EIA_CO2_emis_ind_CONUS = subset(EIA_CO2_emis_ind, State != c("AK", "HI"))

#add myCODE to EIA_CO2_emis_*
EIA_CO2_emis_CONUS = join(EIA_CO2_emis_CONUS, MSN_to_AREAmyCODE, by="MSN", type="left")
colnames(EIA_CO2_emis_CONUS)[colnames(EIA_CO2_emis_CONUS) == "AREAmyCODE"] ="myCODE"

#join EIA_CO2_emis_CONUS with TOTAL_CO_emis_days_frac
EIA_CO2_emis_CONUS = join(EIA_CO2_emis_CONUS,TOTAL_CO_emis_days_frac, by=c("myCODE","State"), type="left")

#divide annual total emissions of CO2 in (metric tons/year) to 52*5 weekdays, 52*1 Saturdays, and 52*1 Sundays
EIA_CO2_emis_CONUS = subset(EIA_CO2_emis_CONUS, select=c("MSN","myCODE","State","CO2_emis_2017","weekdy_emis_frac", "satdy_emis_frac", "sundy_emis_frac"))

EIA_CO2_emis_CONUS['CO2_emis_2017_weekdys'] <- as.numeric(EIA_CO2_emis_CONUS$CO2_emis_2017) * as.numeric(EIA_CO2_emis_CONUS$weekdy_emis_frac)
EIA_CO2_emis_CONUS['CO2_emis_2017_satdys'] <- as.numeric(EIA_CO2_emis_CONUS$CO2_emis_2017) * as.numeric(EIA_CO2_emis_CONUS$satdy_emis_frac)
EIA_CO2_emis_CONUS['CO2_emis_2017_sundys'] <- as.numeric(EIA_CO2_emis_CONUS$CO2_emis_2017) * as.numeric(EIA_CO2_emis_CONUS$sundy_emis_frac)

EIA_CO2_emis_CONUS_days = subset(EIA_CO2_emis_CONUS, select=c("MSN","myCODE","State",
                                                              "CO2_emis_2017_weekdys","CO2_emis_2017_satdys","CO2_emis_2017_sundys"))

#Convert the unit of CO2 emissions in (metric tons/year on the type of day) to (metric tons/day)
num_weekdys = 52*5
num_satdys = 52*1
num_sundys = 52*1

EIA_CO2_emis_CONUS_days['CO2_emis_2017_perweekdy'] = EIA_CO2_emis_CONUS_days$CO2_emis_2017_weekdys / num_weekdys
EIA_CO2_emis_CONUS_days['CO2_emis_2017_persatdy'] = EIA_CO2_emis_CONUS_days$CO2_emis_2017_satdys / num_satdys
EIA_CO2_emis_CONUS_days['CO2_emis_2017_persundy'] = EIA_CO2_emis_CONUS_days$CO2_emis_2017_sundys / num_sundys

EIA_CO2_emis_CONUS_days = subset(EIA_CO2_emis_CONUS_days, select=c("MSN","myCODE","State",
                                                              "CO2_emis_2017_perweekdy","CO2_emis_2017_persatdy","CO2_emis_2017_persundy"))

#########################################################################################################################
#Calculate Area/Point ratio for each Commercial and Industrial myCODE, in each state, for each type of day
#Join AREA with POINT CO emissions
AREA_POINT_sectors_CO_Stateabb = join(AREA_sectors_CO_Stateabb, POINT_sectors_CO_Stateabb, 
                                      by=c("AREAmyCODE","day","STATE_NAME","State"), type="left")

AREA_POINT_sectors_CO_patch_Stateabb = join(AREA_sectors_CO_patch_Stateabb, POINT_sectors_CO_patch_Stateabb, 
                                      by=c("AREAmyCODE","day","STATE_NAME","State"), type="left")
#reorder columns
AREA_POINT_sectors_CO_Stateabb = subset(AREA_POINT_sectors_CO_Stateabb, 
                                        select=c("AREAmyCODE", "POINTmyCODE", 
                                                 "day", "STATE_FIPS_AREA", "STATE_FIPS_POINT",
                                                 "STATE_NAME", "State", "emis_sum_AREA", "emis_sum_POINT"))

AREA_POINT_sectors_CO_patch_Stateabb = subset(AREA_POINT_sectors_CO_patch_Stateabb, 
                                        select=c("AREAmyCODE", "POINTmyCODE", 
                                                 "day", "STATE_FIPS_AREA", "STATE_FIPS_POINT",
                                                 "STATE_NAME", "State", "emis_sum_AREA", "emis_sum_POINT"))

#Make residential point emissions, which are NaNs to zero
AREA_POINT_sectors_CO_Stateabb$emis_sum_POINT[is.na(AREA_POINT_sectors_CO_Stateabb$emis_sum_POINT)] = 0

AREA_POINT_sectors_CO_patch_Stateabb$emis_sum_POINT[is.na(AREA_POINT_sectors_CO_patch_Stateabb$emis_sum_POINT)] = 0

#calculate area/point ratio
AREA_POINT_sectors_CO_Stateabb$emis_sum_TOTAL=AREA_POINT_sectors_CO_Stateabb$emis_sum_AREA + AREA_POINT_sectors_CO_Stateabb$emis_sum_POINT
AREA_POINT_sectors_CO_Stateabb$emis_frac_AREA = AREA_POINT_sectors_CO_Stateabb$emis_sum_AREA/AREA_POINT_sectors_CO_Stateabb$emis_sum_TOTAL
AREA_POINT_sectors_CO_Stateabb$emis_frac_POINT = AREA_POINT_sectors_CO_Stateabb$emis_sum_POINT/AREA_POINT_sectors_CO_Stateabb$emis_sum_TOTAL

AREA_POINT_sectors_CO_patch_Stateabb$emis_sum_TOTAL=AREA_POINT_sectors_CO_patch_Stateabb$emis_sum_AREA + AREA_POINT_sectors_CO_patch_Stateabb$emis_sum_POINT
AREA_POINT_sectors_CO_patch_Stateabb$emis_frac_AREA = AREA_POINT_sectors_CO_patch_Stateabb$emis_sum_AREA/AREA_POINT_sectors_CO_patch_Stateabb$emis_sum_TOTAL
AREA_POINT_sectors_CO_patch_Stateabb$emis_frac_POINT = AREA_POINT_sectors_CO_patch_Stateabb$emis_sum_POINT/AREA_POINT_sectors_CO_patch_Stateabb$emis_sum_TOTAL

###########################################################################
#if state-level emissions > 0, and NaN for both Area_frac and Point_frac (caused by area+point emis = 0) use area/point frac from patch
#reformat AREA_POINT_sectors_CO_patch_Stateabb so it can be joined with AREA_POINT_sectors_CO_Stateabb
#change *_Allfuel to *_Coal _NG _Oil
#RES
AREA_POINT_sectors_CO_patch_Stateabb_RES = AREA_POINT_sectors_CO_patch_Stateabb  %>% filter(AREAmyCODE == 'RES_Allfuel')
AREA_POINT_sectors_CO_patch_Stateabb_RES_Coal = AREA_POINT_sectors_CO_patch_Stateabb_RES
AREA_POINT_sectors_CO_patch_Stateabb_RES_NG = AREA_POINT_sectors_CO_patch_Stateabb_RES
AREA_POINT_sectors_CO_patch_Stateabb_RES_Oil = AREA_POINT_sectors_CO_patch_Stateabb_RES

AREA_POINT_sectors_CO_patch_Stateabb_RES_Coal$AREAmyCODE = 'RES_Coal'
AREA_POINT_sectors_CO_patch_Stateabb_RES_NG$AREAmyCODE = 'RES_NG'
AREA_POINT_sectors_CO_patch_Stateabb_RES_Oil$AREAmyCODE = 'RES_Oil'

AREA_POINT_sectors_CO_patch_Stateabb_RES_specificfuel = rbind(AREA_POINT_sectors_CO_patch_Stateabb_RES_Coal,AREA_POINT_sectors_CO_patch_Stateabb_RES_NG,AREA_POINT_sectors_CO_patch_Stateabb_RES_Oil)

#COMM
AREA_POINT_sectors_CO_patch_Stateabb_COMM = AREA_POINT_sectors_CO_patch_Stateabb  %>% filter(AREAmyCODE == 'COMM_Allfuel')
AREA_POINT_sectors_CO_patch_Stateabb_COMM_Coal = AREA_POINT_sectors_CO_patch_Stateabb_COMM
AREA_POINT_sectors_CO_patch_Stateabb_COMM_NG = AREA_POINT_sectors_CO_patch_Stateabb_COMM
AREA_POINT_sectors_CO_patch_Stateabb_COMM_Oil = AREA_POINT_sectors_CO_patch_Stateabb_COMM

AREA_POINT_sectors_CO_patch_Stateabb_COMM_Coal$AREAmyCODE = 'COMM_Coal'
AREA_POINT_sectors_CO_patch_Stateabb_COMM_NG$AREAmyCODE = 'COMM_NG'
AREA_POINT_sectors_CO_patch_Stateabb_COMM_Oil$AREAmyCODE = 'COMM_Oil'

AREA_POINT_sectors_CO_patch_Stateabb_COMM_Coal$POINTmyCODE = 'PtCOMM_Coal'
AREA_POINT_sectors_CO_patch_Stateabb_COMM_NG$POINTmyCODE = 'PtCOMM_NG'
AREA_POINT_sectors_CO_patch_Stateabb_COMM_Oil$POINTmyCODE = 'PtCOMM_Oil'

AREA_POINT_sectors_CO_patch_Stateabb_COMM_specificfuel = rbind(AREA_POINT_sectors_CO_patch_Stateabb_COMM_Coal,AREA_POINT_sectors_CO_patch_Stateabb_COMM_NG,AREA_POINT_sectors_CO_patch_Stateabb_COMM_Oil)

#IND
AREA_POINT_sectors_CO_patch_Stateabb_IND = AREA_POINT_sectors_CO_patch_Stateabb  %>% filter(AREAmyCODE == 'IND_Allfuel')
AREA_POINT_sectors_CO_patch_Stateabb_IND_Coal = AREA_POINT_sectors_CO_patch_Stateabb_IND
AREA_POINT_sectors_CO_patch_Stateabb_IND_NG = AREA_POINT_sectors_CO_patch_Stateabb_IND
AREA_POINT_sectors_CO_patch_Stateabb_IND_Oil = AREA_POINT_sectors_CO_patch_Stateabb_IND

AREA_POINT_sectors_CO_patch_Stateabb_IND_Coal$AREAmyCODE = 'IND_Coal'
AREA_POINT_sectors_CO_patch_Stateabb_IND_NG$AREAmyCODE = 'IND_NG'
AREA_POINT_sectors_CO_patch_Stateabb_IND_Oil$AREAmyCODE = 'IND_Oil'

AREA_POINT_sectors_CO_patch_Stateabb_IND_Coal$POINTmyCODE = 'PtIND_Coal'
AREA_POINT_sectors_CO_patch_Stateabb_IND_NG$POINTmyCODE = 'PtIND_NG'
AREA_POINT_sectors_CO_patch_Stateabb_IND_Oil$POINTmyCODE = 'PtIND_Oil'

AREA_POINT_sectors_CO_patch_Stateabb_IND_specificfuel = rbind(AREA_POINT_sectors_CO_patch_Stateabb_IND_Coal,AREA_POINT_sectors_CO_patch_Stateabb_IND_NG,AREA_POINT_sectors_CO_patch_Stateabb_IND_Oil)

AREA_POINT_sectors_CO_patch_Stateabb_specificfuel = rbind(AREA_POINT_sectors_CO_patch_Stateabb_RES_specificfuel,AREA_POINT_sectors_CO_patch_Stateabb_COMM_specificfuel,AREA_POINT_sectors_CO_patch_Stateabb_IND_specificfuel)

colnames(AREA_POINT_sectors_CO_patch_Stateabb_specificfuel)[colnames(AREA_POINT_sectors_CO_patch_Stateabb_specificfuel) == "emis_sum_AREA"] ="emis_sum_AREA_patch"
colnames(AREA_POINT_sectors_CO_patch_Stateabb_specificfuel)[colnames(AREA_POINT_sectors_CO_patch_Stateabb_specificfuel) == "emis_sum_POINT"] ="emis_sum_POINT_patch"
colnames(AREA_POINT_sectors_CO_patch_Stateabb_specificfuel)[colnames(AREA_POINT_sectors_CO_patch_Stateabb_specificfuel) == "emis_sum_TOTAL"] ="emis_sum_TOTAL_patch"
colnames(AREA_POINT_sectors_CO_patch_Stateabb_specificfuel)[colnames(AREA_POINT_sectors_CO_patch_Stateabb_specificfuel) == "emis_frac_AREA"] ="emis_frac_AREA_patch"
colnames(AREA_POINT_sectors_CO_patch_Stateabb_specificfuel)[colnames(AREA_POINT_sectors_CO_patch_Stateabb_specificfuel) == "emis_frac_POINT"] ="emis_frac_POINT_patch"

#Join AREA_POINT_sectors_CO_patch_Stateabb_specificfuel with AREA_POINT_sectors_CO_Stateabb
AREA_POINT_sectors_CO_Stateabb = join(AREA_POINT_sectors_CO_Stateabb,AREA_POINT_sectors_CO_patch_Stateabb_specificfuel, by=c("AREAmyCODE","POINTmyCODE","day","STATE_FIPS_AREA","STATE_FIPS_POINT","STATE_NAME","State"), type="left")
AREA_POINT_sectors_CO_Stateabb = subset(AREA_POINT_sectors_CO_Stateabb, select=c("AREAmyCODE","POINTmyCODE","day","STATE_FIPS_AREA","STATE_FIPS_POINT","STATE_NAME","State","emis_frac_AREA","emis_frac_AREA_patch","emis_frac_POINT","emis_frac_POINT_patch"))

#reformat EIA_CO2_emis_CONUS_days so it can join with AREA_POINT_sectors_CO_Stateabb
colnames(EIA_CO2_emis_CONUS_days)[colnames(EIA_CO2_emis_CONUS_days) == "myCODE"] ="AREAmyCODE"

#Get number of MSN code
num_MSNs = as.numeric(length(MSN_co2_fct$MSN))

EIA_CO2_emis_CONUS_days_reform <- data.frame(MSN=character(num_states_CONUS*num_MSNs*3), 
                                             AREAmyCODE=character(num_states_CONUS*num_MSNs*3), 
                                             State=character(num_states_CONUS*num_MSNs*3),
                                             day=character(num_states_CONUS*num_MSNs*3), 
                                             CO2_emis_2017_perdy=character(num_states_CONUS*num_MSNs*3))
                                             
#2017
EIA_CO2_emis_CONUS_days_reform$MSN[1:(num_states_CONUS*num_MSNs)] = EIA_CO2_emis_CONUS_days$MSN
EIA_CO2_emis_CONUS_days_reform$AREAmyCODE[1:(num_states_CONUS*num_MSNs)] = EIA_CO2_emis_CONUS_days$AREAmyCODE
EIA_CO2_emis_CONUS_days_reform$State[1:(num_states_CONUS*num_MSNs)] = EIA_CO2_emis_CONUS_days$State
EIA_CO2_emis_CONUS_days_reform$day[1:(num_states_CONUS*num_MSNs)] = 'weekdy'
EIA_CO2_emis_CONUS_days_reform$CO2_emis_2017_perdy[1:(num_states_CONUS*num_MSNs)] = EIA_CO2_emis_CONUS_days$CO2_emis_2017_perweekdy

EIA_CO2_emis_CONUS_days_reform$MSN[(num_states_CONUS*num_MSNs+1):(num_states_CONUS*num_MSNs*2)] = EIA_CO2_emis_CONUS_days$MSN
EIA_CO2_emis_CONUS_days_reform$AREAmyCODE[(num_states_CONUS*num_MSNs+1):(num_states_CONUS*num_MSNs*2)] = EIA_CO2_emis_CONUS_days$AREAmyCODE
EIA_CO2_emis_CONUS_days_reform$State[(num_states_CONUS*num_MSNs+1):(num_states_CONUS*num_MSNs*2)] = EIA_CO2_emis_CONUS_days$State
EIA_CO2_emis_CONUS_days_reform$day[(num_states_CONUS*num_MSNs+1):(num_states_CONUS*num_MSNs*2)] = 'satdy'
EIA_CO2_emis_CONUS_days_reform$CO2_emis_2017_perdy[(num_states_CONUS*num_MSNs+1):(num_states_CONUS*num_MSNs*2)] = EIA_CO2_emis_CONUS_days$CO2_emis_2017_persatdy

EIA_CO2_emis_CONUS_days_reform$MSN[(num_states_CONUS*num_MSNs*2+1):(num_states_CONUS*num_MSNs*3)] = EIA_CO2_emis_CONUS_days$MSN
EIA_CO2_emis_CONUS_days_reform$AREAmyCODE[(num_states_CONUS*num_MSNs*2+1):(num_states_CONUS*num_MSNs*3)] = EIA_CO2_emis_CONUS_days$AREAmyCODE
EIA_CO2_emis_CONUS_days_reform$State[(num_states_CONUS*num_MSNs*2+1):(num_states_CONUS*num_MSNs*3)] = EIA_CO2_emis_CONUS_days$State
EIA_CO2_emis_CONUS_days_reform$day[(num_states_CONUS*num_MSNs*2+1):(num_states_CONUS*num_MSNs*3)] = 'sundy'
EIA_CO2_emis_CONUS_days_reform$CO2_emis_2017_perdy[(num_states_CONUS*num_MSNs*2+1):(num_states_CONUS*num_MSNs*3)] = EIA_CO2_emis_CONUS_days$CO2_emis_2017_persundy

#join EIA_CO2_emis_CONUS_days_reform with AREA_POINT_CO_emis_frac
EIA_CO2_emis_CONUS_days_reform = join(EIA_CO2_emis_CONUS_days_reform,AREA_POINT_sectors_CO_Stateabb, by=c("AREAmyCODE","State", "day"), type="left")

#Create a record of what are patched for creating the patched spatial surrogates
AREA_POINT_sumCO_what_need_patch = EIA_CO2_emis_CONUS_days_reform  %>% filter((CO2_emis_2017_perdy!=0)&is.nan(emis_frac_AREA)&is.nan(emis_frac_POINT))

#patch
for (r in  1:nrow(EIA_CO2_emis_CONUS_days_reform))
{
  if ((EIA_CO2_emis_CONUS_days_reform$CO2_emis_2017_perdy[r]!=0) & is.nan(EIA_CO2_emis_CONUS_days_reform$emis_frac_AREA[r]) & is.nan(EIA_CO2_emis_CONUS_days_reform$emis_frac_POINT[r]))
  { 
    print(EIA_CO2_emis_CONUS_days_reform$AREAmyCODE[r])
    print(EIA_CO2_emis_CONUS_days_reform$State[r])
    EIA_CO2_emis_CONUS_days_reform$emis_frac_AREA[r] = EIA_CO2_emis_CONUS_days_reform$emis_frac_AREA_patch[r]
    EIA_CO2_emis_CONUS_days_reform$emis_frac_POINT[r] = EIA_CO2_emis_CONUS_days_reform$emis_frac_POINT_patch[r]
  }
}

#reorder columns
EIA_CO2_emis_CONUS_days_reform = subset(EIA_CO2_emis_CONUS_days_reform, select=c("MSN", "AREAmyCODE", "POINTmyCODE", 
                                                                                 "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME", "State",
                                                                                 "day", 
                                                                                 "CO2_emis_2017_perdy",
                                                                                 "emis_frac_AREA", "emis_frac_POINT"))

#make NaN to zero, then NaN point area frac only exist for zero emissions, which do not need to be distributed into area/point
EIA_CO2_emis_CONUS_days_reform$emis_frac_POINT[is.na(EIA_CO2_emis_CONUS_days_reform$emis_frac_POINT)] = 0
EIA_CO2_emis_CONUS_days_reform$emis_frac_AREA[is.na(EIA_CO2_emis_CONUS_days_reform$emis_frac_AREA)] = 0

#divide emissions into area and point
EIA_CO2_emis_CONUS_days_reform$CO2_emis_2017_perdy_area = as.numeric(EIA_CO2_emis_CONUS_days_reform$CO2_emis_2017_perdy) * as.numeric(EIA_CO2_emis_CONUS_days_reform$emis_frac_AREA)
EIA_CO2_emis_CONUS_days_reform$CO2_emis_2017_perdy_point = as.numeric(EIA_CO2_emis_CONUS_days_reform$CO2_emis_2017_perdy) * as.numeric(EIA_CO2_emis_CONUS_days_reform$emis_frac_POINT)

#do a sanity check
CO2_emis_perdy <- sum(as.numeric(EIA_CO2_emis_CONUS_days_reform$CO2_emis_2017_perdy))
CO2_emis_perdy_area <- sum(as.numeric(EIA_CO2_emis_CONUS_days_reform$CO2_emis_2017_perdy_area))
CO2_emis_perdy_point <- sum(as.numeric(EIA_CO2_emis_CONUS_days_reform$CO2_emis_2017_perdy_point))
CO2_emis_perdy_total = CO2_emis_perdy_area + CO2_emis_perdy_point
#no problem

#simplify EIA_CO2_emis_CONUS_days_reform
EIA_CO2_emis_CONUS_days_reform = subset(EIA_CO2_emis_CONUS_days_reform, select=c("MSN", "AREAmyCODE", "POINTmyCODE", 
                                                                                 "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                 "day", 
                                                                                 "CO2_emis_2017_perdy_area",
                                                                                 "CO2_emis_2017_perdy_point"))

#################################################################################################################################################
#make final result in the format as Colin's CO files
#sector(myCODE), day, STATE_FIPS, STATE_NAME, emis_sum2017

#sum 35 MSN into 12 AREAmyCODE and 8 POINTmyCODE
AREA_sectors_CO2 <- data.frame(sector=character(num_states_CONUS*num_area_myCODEs*3), day=character(num_states_CONUS*num_area_myCODEs*3), STATE_FIPS=character(num_states_CONUS*num_area_myCODEs*3),
                                           STATE_NAME=character(num_states_CONUS*num_area_myCODEs*3), emis_sum2017=character(num_states_CONUS*num_area_myCODEs*3))

POINT_sectors_CO2 <- data.frame(sector=character(num_states_CONUS*num_point_myCODEs*3), day=character(num_states_CONUS*num_point_myCODEs*3), STATE_FIPS=character(num_states_CONUS*num_point_myCODEs*3),
                               STATE_NAME=character(num_states_CONUS*num_point_myCODEs*3), emis_sum2017=character(num_states_CONUS*num_point_myCODEs*3))

areamyCODE_num = 1
ptmyCODE_num = 1
for (j in area_myCODEs){
  print(j)
  
  if (j == 'COMM_NG'){
    EIA_CO2_emis_CONUS_days_MSN1 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'NGCCB')
    EIA_CO2_emis_CONUS_days_MSN2 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'HLCCB')
    
    names(EIA_CO2_emis_CONUS_days_MSN1)[1] <- 'MSN1'
    names(EIA_CO2_emis_CONUS_days_MSN1)[8] <- 'CO2_emis_2017_perdy_area1'
    names(EIA_CO2_emis_CONUS_days_MSN1)[9] <- 'CO2_emis_2017_perdy_point1'
    
    names(EIA_CO2_emis_CONUS_days_MSN2)[1] <- 'MSN2'
    names(EIA_CO2_emis_CONUS_days_MSN2)[8] <- 'CO2_emis_2017_perdy_area2'
    names(EIA_CO2_emis_CONUS_days_MSN2)[9] <- 'CO2_emis_2017_perdy_point2'
    
    EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_MSN1,EIA_CO2_emis_CONUS_days_MSN2, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                      "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                       "day"), type="left")
    
    AREA_sectors_CO2$sector[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$AREAmyCODE
    AREA_sectors_CO2$day[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$day
    AREA_sectors_CO2$STATE_FIPS[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_FIPS_AREA
    AREA_sectors_CO2$STATE_NAME[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_NAME
    AREA_sectors_CO2$emis_sum2017[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area1 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area2
    
    POINT_sectors_CO2$sector[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$POINTmyCODE
    POINT_sectors_CO2$day[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$day
    POINT_sectors_CO2$STATE_FIPS[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_FIPS_POINT
    POINT_sectors_CO2$STATE_NAME[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_NAME
    POINT_sectors_CO2$emis_sum2017[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point1 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point2
    
    areamyCODE_num = areamyCODE_num + 1
    ptmyCODE_num = ptmyCODE_num + 1
  }
  
  if (j == 'COMM_Coal'){
    EIA_CO2_emis_CONUS_days_MSN1 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'CLCCB')
    EIA_CO2_emis_CONUS_days_MSN2 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'PCCCB')

    names(EIA_CO2_emis_CONUS_days_MSN1)[1] <- 'MSN1'
    names(EIA_CO2_emis_CONUS_days_MSN1)[8] <- 'CO2_emis_2017_perdy_area1'
    names(EIA_CO2_emis_CONUS_days_MSN1)[9] <- 'CO2_emis_2017_perdy_point1'
    
    names(EIA_CO2_emis_CONUS_days_MSN2)[1] <- 'MSN2'
    names(EIA_CO2_emis_CONUS_days_MSN2)[8] <- 'CO2_emis_2017_perdy_area2'
    names(EIA_CO2_emis_CONUS_days_MSN2)[9] <- 'CO2_emis_2017_perdy_point2'
    
    EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_MSN1,EIA_CO2_emis_CONUS_days_MSN2, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                      "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                       "day"), type="left")

    AREA_sectors_CO2$sector[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$AREAmyCODE
    AREA_sectors_CO2$day[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$day
    AREA_sectors_CO2$STATE_FIPS[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_FIPS_AREA
    AREA_sectors_CO2$STATE_NAME[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_NAME
    AREA_sectors_CO2$emis_sum2017[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area1 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area2
    
    POINT_sectors_CO2$sector[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$POINTmyCODE
    POINT_sectors_CO2$day[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$day
    POINT_sectors_CO2$STATE_FIPS[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_FIPS_POINT
    POINT_sectors_CO2$STATE_NAME[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_NAME
    POINT_sectors_CO2$emis_sum2017[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point1 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point2
    
    areamyCODE_num = areamyCODE_num + 1
    ptmyCODE_num = ptmyCODE_num + 1
  }
  
  if (j == 'COMM_Oil'){
    EIA_CO2_emis_CONUS_days_MSN1 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'DFCCB')
    EIA_CO2_emis_CONUS_days_MSN2 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'RFCCB')
    EIA_CO2_emis_CONUS_days_MSN3 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'KSCCB')
    EIA_CO2_emis_CONUS_days_MSN4 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'MGCCB')

    names(EIA_CO2_emis_CONUS_days_MSN1)[1] <- 'MSN1'
    names(EIA_CO2_emis_CONUS_days_MSN1)[8] <- 'CO2_emis_2017_perdy_area1'
    names(EIA_CO2_emis_CONUS_days_MSN1)[9] <- 'CO2_emis_2017_perdy_point1'
    
    names(EIA_CO2_emis_CONUS_days_MSN2)[1] <- 'MSN2'
    names(EIA_CO2_emis_CONUS_days_MSN2)[8] <- 'CO2_emis_2017_perdy_area2'
    names(EIA_CO2_emis_CONUS_days_MSN2)[9] <- 'CO2_emis_2017_perdy_point2'

    names(EIA_CO2_emis_CONUS_days_MSN3)[1] <- 'MSN3'
    names(EIA_CO2_emis_CONUS_days_MSN3)[8] <- 'CO2_emis_2017_perdy_area3'
    names(EIA_CO2_emis_CONUS_days_MSN3)[9] <- 'CO2_emis_2017_perdy_point3'

    names(EIA_CO2_emis_CONUS_days_MSN4)[1] <- 'MSN4'
    names(EIA_CO2_emis_CONUS_days_MSN4)[8] <- 'CO2_emis_2017_perdy_area4'
    names(EIA_CO2_emis_CONUS_days_MSN4)[9] <- 'CO2_emis_2017_perdy_point4'
    
    EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_MSN1,EIA_CO2_emis_CONUS_days_MSN2, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                     "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                     "day"), type="left")
    EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_j,EIA_CO2_emis_CONUS_days_MSN3, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                     "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                     "day"), type="left")
    EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_j,EIA_CO2_emis_CONUS_days_MSN4, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                     "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                     "day"), type="left")
    
    AREA_sectors_CO2$sector[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$AREAmyCODE
    AREA_sectors_CO2$day[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$day
    AREA_sectors_CO2$STATE_FIPS[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_FIPS_AREA
    AREA_sectors_CO2$STATE_NAME[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_NAME
    AREA_sectors_CO2$emis_sum2017[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area1 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area2 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area3 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area4
    
    POINT_sectors_CO2$sector[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$POINTmyCODE
    POINT_sectors_CO2$day[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$day
    POINT_sectors_CO2$STATE_FIPS[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_FIPS_POINT
    POINT_sectors_CO2$STATE_NAME[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_NAME
    POINT_sectors_CO2$emis_sum2017[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point1 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point2 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point3 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point4 
    
    areamyCODE_num = areamyCODE_num + 1
    ptmyCODE_num = ptmyCODE_num + 1
  }

  if (j == 'IND_NG'){
    EIA_CO2_emis_CONUS_days_MSN1 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'NGICB')
    EIA_CO2_emis_CONUS_days_MSN2 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'HLICB')
    EIA_CO2_emis_CONUS_days_MSN3 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'PPICB')
    EIA_CO2_emis_CONUS_days_MSN4 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'SGICB')

    names(EIA_CO2_emis_CONUS_days_MSN1)[1] <- 'MSN1'
    names(EIA_CO2_emis_CONUS_days_MSN1)[8] <- 'CO2_emis_2017_perdy_area1'
    names(EIA_CO2_emis_CONUS_days_MSN1)[9] <- 'CO2_emis_2017_perdy_point1'
    
    names(EIA_CO2_emis_CONUS_days_MSN2)[1] <- 'MSN2'
    names(EIA_CO2_emis_CONUS_days_MSN2)[8] <- 'CO2_emis_2017_perdy_area2'
    names(EIA_CO2_emis_CONUS_days_MSN2)[9] <- 'CO2_emis_2017_perdy_point2'

    names(EIA_CO2_emis_CONUS_days_MSN3)[1] <- 'MSN3'
    names(EIA_CO2_emis_CONUS_days_MSN3)[8] <- 'CO2_emis_2017_perdy_area3'
    names(EIA_CO2_emis_CONUS_days_MSN3)[9] <- 'CO2_emis_2017_perdy_point3'

    names(EIA_CO2_emis_CONUS_days_MSN4)[1] <- 'MSN4'
    names(EIA_CO2_emis_CONUS_days_MSN4)[8] <- 'CO2_emis_2017_perdy_area4'
    names(EIA_CO2_emis_CONUS_days_MSN4)[9] <- 'CO2_emis_2017_perdy_point4'
    
    EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_MSN1,EIA_CO2_emis_CONUS_days_MSN2, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                     "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                     "day"), type="left")
    EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_j,EIA_CO2_emis_CONUS_days_MSN3, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                     "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                     "day"), type="left")
    EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_j,EIA_CO2_emis_CONUS_days_MSN4, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                     "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                     "day"), type="left")
    
    AREA_sectors_CO2$sector[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$AREAmyCODE
    AREA_sectors_CO2$day[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$day
    AREA_sectors_CO2$STATE_FIPS[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_FIPS_AREA
    AREA_sectors_CO2$STATE_NAME[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_NAME
    AREA_sectors_CO2$emis_sum2017[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area1 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area2 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area3 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area4
    
    POINT_sectors_CO2$sector[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$POINTmyCODE
    POINT_sectors_CO2$day[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$day
    POINT_sectors_CO2$STATE_FIPS[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_FIPS_POINT
    POINT_sectors_CO2$STATE_NAME[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_NAME
    POINT_sectors_CO2$emis_sum2017[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point1 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point2 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point3 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point4
    
    areamyCODE_num = areamyCODE_num + 1
    ptmyCODE_num = ptmyCODE_num + 1
  }
  
  if (j == 'IND_Oil'){
    EIA_CO2_emis_CONUS_days_MSN1 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'DFICB')
    EIA_CO2_emis_CONUS_days_MSN2 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'RFICB')
    EIA_CO2_emis_CONUS_days_MSN3 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'KSICB')
    EIA_CO2_emis_CONUS_days_MSN4 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'ABICB')
    EIA_CO2_emis_CONUS_days_MSN5 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'ARICB')
    EIA_CO2_emis_CONUS_days_MSN6 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'COICB')
    EIA_CO2_emis_CONUS_days_MSN7 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'FNICB')
    EIA_CO2_emis_CONUS_days_MSN8 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'FOICB')
    EIA_CO2_emis_CONUS_days_MSN9 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'LUICB')
    EIA_CO2_emis_CONUS_days_MSN10 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'MBICB')
    EIA_CO2_emis_CONUS_days_MSN11 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'MGICB')
    EIA_CO2_emis_CONUS_days_MSN12 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'MSICB')
    EIA_CO2_emis_CONUS_days_MSN13 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'SNICB')
    EIA_CO2_emis_CONUS_days_MSN14 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'UOICB')
    EIA_CO2_emis_CONUS_days_MSN15 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'WXICB')

    names(EIA_CO2_emis_CONUS_days_MSN1)[1] <- 'MSN1'
    names(EIA_CO2_emis_CONUS_days_MSN1)[8] <- 'CO2_emis_2017_perdy_area1'
    names(EIA_CO2_emis_CONUS_days_MSN1)[9] <- 'CO2_emis_2017_perdy_point1'   
    
    names(EIA_CO2_emis_CONUS_days_MSN2)[1] <- 'MSN2'
    names(EIA_CO2_emis_CONUS_days_MSN2)[8] <- 'CO2_emis_2017_perdy_area2'
    names(EIA_CO2_emis_CONUS_days_MSN2)[9] <- 'CO2_emis_2017_perdy_point2'

    names(EIA_CO2_emis_CONUS_days_MSN3)[1] <- 'MSN3'
    names(EIA_CO2_emis_CONUS_days_MSN3)[8] <- 'CO2_emis_2017_perdy_area3'
    names(EIA_CO2_emis_CONUS_days_MSN3)[9] <- 'CO2_emis_2017_perdy_point3'

    names(EIA_CO2_emis_CONUS_days_MSN4)[1] <- 'MSN4'
    names(EIA_CO2_emis_CONUS_days_MSN4)[8] <- 'CO2_emis_2017_perdy_area4'
    names(EIA_CO2_emis_CONUS_days_MSN4)[9] <- 'CO2_emis_2017_perdy_point4'

    names(EIA_CO2_emis_CONUS_days_MSN5)[1] <- 'MSN5'
    names(EIA_CO2_emis_CONUS_days_MSN5)[8] <- 'CO2_emis_2017_perdy_area5'
    names(EIA_CO2_emis_CONUS_days_MSN5)[9] <- 'CO2_emis_2017_perdy_point5'

    names(EIA_CO2_emis_CONUS_days_MSN6)[1] <- 'MSN6'
    names(EIA_CO2_emis_CONUS_days_MSN6)[8] <- 'CO2_emis_2017_perdy_area6'
    names(EIA_CO2_emis_CONUS_days_MSN6)[9] <- 'CO2_emis_2017_perdy_point6'

    names(EIA_CO2_emis_CONUS_days_MSN7)[1] <- 'MSN7'
    names(EIA_CO2_emis_CONUS_days_MSN7)[8] <- 'CO2_emis_2017_perdy_area7'
    names(EIA_CO2_emis_CONUS_days_MSN7)[9] <- 'CO2_emis_2017_perdy_point7'

    names(EIA_CO2_emis_CONUS_days_MSN8)[1] <- 'MSN8'
    names(EIA_CO2_emis_CONUS_days_MSN8)[8] <- 'CO2_emis_2017_perdy_area8'
    names(EIA_CO2_emis_CONUS_days_MSN8)[9] <- 'CO2_emis_2017_perdy_point8'

    names(EIA_CO2_emis_CONUS_days_MSN9)[1] <- 'MSN9'
    names(EIA_CO2_emis_CONUS_days_MSN9)[8] <- 'CO2_emis_2017_perdy_area9'
    names(EIA_CO2_emis_CONUS_days_MSN9)[9] <- 'CO2_emis_2017_perdy_point9'

    names(EIA_CO2_emis_CONUS_days_MSN10)[1] <- 'MSN10'
    names(EIA_CO2_emis_CONUS_days_MSN10)[8] <- 'CO2_emis_2017_perdy_area10'
    names(EIA_CO2_emis_CONUS_days_MSN10)[9] <- 'CO2_emis_2017_perdy_point10'

    names(EIA_CO2_emis_CONUS_days_MSN11)[1] <- 'MSN11'
    names(EIA_CO2_emis_CONUS_days_MSN11)[8] <- 'CO2_emis_2017_perdy_area11'
    names(EIA_CO2_emis_CONUS_days_MSN11)[9] <- 'CO2_emis_2017_perdy_point11'

    names(EIA_CO2_emis_CONUS_days_MSN12)[1] <- 'MSN12'
    names(EIA_CO2_emis_CONUS_days_MSN12)[8] <- 'CO2_emis_2017_perdy_area12'
    names(EIA_CO2_emis_CONUS_days_MSN12)[9] <- 'CO2_emis_2017_perdy_point12'

    names(EIA_CO2_emis_CONUS_days_MSN13)[1] <- 'MSN13'
    names(EIA_CO2_emis_CONUS_days_MSN13)[8] <- 'CO2_emis_2017_perdy_area13'
    names(EIA_CO2_emis_CONUS_days_MSN13)[9] <- 'CO2_emis_2017_perdy_point13'

    names(EIA_CO2_emis_CONUS_days_MSN14)[1] <- 'MSN14'
    names(EIA_CO2_emis_CONUS_days_MSN14)[8] <- 'CO2_emis_2017_perdy_area14'
    names(EIA_CO2_emis_CONUS_days_MSN14)[9] <- 'CO2_emis_2017_perdy_point14'

    names(EIA_CO2_emis_CONUS_days_MSN15)[1] <- 'MSN15'
    names(EIA_CO2_emis_CONUS_days_MSN15)[8] <- 'CO2_emis_2017_perdy_area15'
    names(EIA_CO2_emis_CONUS_days_MSN15)[9] <- 'CO2_emis_2017_perdy_point15'
    
    EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_MSN1,EIA_CO2_emis_CONUS_days_MSN2, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                     "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                     "day"), type="left")
    EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_j,EIA_CO2_emis_CONUS_days_MSN3, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                  "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                  "day"), type="left")
    EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_j,EIA_CO2_emis_CONUS_days_MSN4, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                  "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                  "day"), type="left")
    EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_j,EIA_CO2_emis_CONUS_days_MSN5, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                  "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                  "day"), type="left")
    EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_j,EIA_CO2_emis_CONUS_days_MSN6, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                  "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                  "day"), type="left")
    EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_j,EIA_CO2_emis_CONUS_days_MSN7, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                  "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                  "day"), type="left")
    EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_j,EIA_CO2_emis_CONUS_days_MSN8, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                  "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                  "day"), type="left")
    EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_j,EIA_CO2_emis_CONUS_days_MSN9, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                  "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                  "day"), type="left")
    EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_j,EIA_CO2_emis_CONUS_days_MSN10, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                  "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                  "day"), type="left")
     EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_j,EIA_CO2_emis_CONUS_days_MSN11, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                  "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                  "day"), type="left")
     EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_j,EIA_CO2_emis_CONUS_days_MSN12, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                  "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                  "day"), type="left")
     EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_j,EIA_CO2_emis_CONUS_days_MSN13, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                  "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                  "day"), type="left")
     EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_j,EIA_CO2_emis_CONUS_days_MSN14, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                  "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                  "day"), type="left")
     EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_j,EIA_CO2_emis_CONUS_days_MSN15, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                  "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                  "day"), type="left")                                                            
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
    AREA_sectors_CO2$sector[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$AREAmyCODE
    AREA_sectors_CO2$day[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$day
    AREA_sectors_CO2$STATE_FIPS[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_FIPS_AREA
    AREA_sectors_CO2$STATE_NAME[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_NAME
    AREA_sectors_CO2$emis_sum2017[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area1 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area2 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area3 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area4 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area5 +
                                                                             EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area6 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area7 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area8 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area9 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area10 +
                                                                             EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area11 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area12 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area13 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area14 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area15
    
    POINT_sectors_CO2$sector[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$POINTmyCODE
    POINT_sectors_CO2$day[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$day
    POINT_sectors_CO2$STATE_FIPS[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_FIPS_POINT
    POINT_sectors_CO2$STATE_NAME[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_NAME
    POINT_sectors_CO2$emis_sum2017[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point1 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point2 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point3 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point4 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point5 +
                                                                              EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point6 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point7 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point8 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point9 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point10 +
                                                                              EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point11 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point12 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point13 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point14 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point15
    
    areamyCODE_num = areamyCODE_num + 1
    ptmyCODE_num = ptmyCODE_num + 1
  }

  if (j == 'IND_Coal'){
    EIA_CO2_emis_CONUS_days_MSN1 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'CLKCB')
    EIA_CO2_emis_CONUS_days_MSN2 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'CLOCB')
    EIA_CO2_emis_CONUS_days_MSN3 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'PCICB')

    names(EIA_CO2_emis_CONUS_days_MSN1)[1] <- 'MSN1'
    names(EIA_CO2_emis_CONUS_days_MSN1)[8] <- 'CO2_emis_2017_perdy_area1'
    names(EIA_CO2_emis_CONUS_days_MSN1)[9] <- 'CO2_emis_2017_perdy_point1'   
    
    names(EIA_CO2_emis_CONUS_days_MSN2)[1] <- 'MSN2'
    names(EIA_CO2_emis_CONUS_days_MSN2)[8] <- 'CO2_emis_2017_perdy_area2'
    names(EIA_CO2_emis_CONUS_days_MSN2)[9] <- 'CO2_emis_2017_perdy_point2'

    names(EIA_CO2_emis_CONUS_days_MSN3)[1] <- 'MSN3'
    names(EIA_CO2_emis_CONUS_days_MSN3)[8] <- 'CO2_emis_2017_perdy_area3'
    names(EIA_CO2_emis_CONUS_days_MSN3)[9] <- 'CO2_emis_2017_perdy_point3'

    EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_MSN1,EIA_CO2_emis_CONUS_days_MSN2, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                     "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                     "day"), type="left")
    EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_j,EIA_CO2_emis_CONUS_days_MSN3, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                  "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                  "day"), type="left")
    
    AREA_sectors_CO2$sector[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$AREAmyCODE
    AREA_sectors_CO2$day[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$day
    AREA_sectors_CO2$STATE_FIPS[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_FIPS_AREA
    AREA_sectors_CO2$STATE_NAME[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_NAME
    AREA_sectors_CO2$emis_sum2017[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area1 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area2 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area3
    
    POINT_sectors_CO2$sector[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$POINTmyCODE
    POINT_sectors_CO2$day[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$day
    POINT_sectors_CO2$STATE_FIPS[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_FIPS_POINT
    POINT_sectors_CO2$STATE_NAME[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_NAME
    POINT_sectors_CO2$emis_sum2017[((ptmyCODE_num-1)*150+1):(ptmyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point1 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point2 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_point3
  
    areamyCODE_num = areamyCODE_num + 1
    ptmyCODE_num = ptmyCODE_num + 1
  }

  if (j == 'RES_Coal'){
    EIA_CO2_emis_CONUS_days_MSN = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'CLRCB')

    EIA_CO2_emis_CONUS_days_j = EIA_CO2_emis_CONUS_days_MSN
    
    AREA_sectors_CO2$sector[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$AREAmyCODE
    AREA_sectors_CO2$day[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$day
    AREA_sectors_CO2$STATE_FIPS[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_FIPS_AREA
    AREA_sectors_CO2$STATE_NAME[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_NAME
    AREA_sectors_CO2$emis_sum2017[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area
    
    areamyCODE_num = areamyCODE_num + 1
  }
  
  if (j == 'RES_NG'){
    EIA_CO2_emis_CONUS_days_MSN1 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'NGRCB')
    EIA_CO2_emis_CONUS_days_MSN2 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'HLRCB')
    
    names(EIA_CO2_emis_CONUS_days_MSN1)[1] <- 'MSN1'
    names(EIA_CO2_emis_CONUS_days_MSN1)[8] <- 'CO2_emis_2017_perdy_area1'
    names(EIA_CO2_emis_CONUS_days_MSN1)[9] <- 'CO2_emis_2017_perdy_point1'
    
    names(EIA_CO2_emis_CONUS_days_MSN2)[1] <- 'MSN2'
    names(EIA_CO2_emis_CONUS_days_MSN2)[8] <- 'CO2_emis_2017_perdy_area2'
    names(EIA_CO2_emis_CONUS_days_MSN2)[9] <- 'CO2_emis_2017_perdy_point2'

    EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_MSN1,EIA_CO2_emis_CONUS_days_MSN2, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                     "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                     "day"), type="left")

    AREA_sectors_CO2$sector[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$AREAmyCODE
    AREA_sectors_CO2$day[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$day
    AREA_sectors_CO2$STATE_FIPS[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_FIPS_AREA
    AREA_sectors_CO2$STATE_NAME[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_NAME
    AREA_sectors_CO2$emis_sum2017[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area1 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area2
  
    areamyCODE_num = areamyCODE_num + 1
  }

  if (j == 'RES_Oil'){
    EIA_CO2_emis_CONUS_days_MSN1 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'DFRCB')
    EIA_CO2_emis_CONUS_days_MSN2 = EIA_CO2_emis_CONUS_days_reform  %>% filter(MSN == 'KSRCB')

    names(EIA_CO2_emis_CONUS_days_MSN1)[1] <- 'MSN1'
    names(EIA_CO2_emis_CONUS_days_MSN1)[8] <- 'CO2_emis_2017_perdy_area1'
    names(EIA_CO2_emis_CONUS_days_MSN1)[9] <- 'CO2_emis_2017_perdy_point1'
    
    names(EIA_CO2_emis_CONUS_days_MSN2)[1] <- 'MSN2'
    names(EIA_CO2_emis_CONUS_days_MSN2)[8] <- 'CO2_emis_2017_perdy_area2'
    names(EIA_CO2_emis_CONUS_days_MSN2)[9] <- 'CO2_emis_2017_perdy_point2'
    
    EIA_CO2_emis_CONUS_days_j = join(EIA_CO2_emis_CONUS_days_MSN1,EIA_CO2_emis_CONUS_days_MSN2, by=c("AREAmyCODE", "POINTmyCODE",
                                                                                                     "STATE_FIPS_AREA", "STATE_FIPS_POINT", "STATE_NAME",
                                                                                                     "day"), type="left")
    
    AREA_sectors_CO2$sector[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$AREAmyCODE
    AREA_sectors_CO2$day[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$day
    AREA_sectors_CO2$STATE_FIPS[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_FIPS_AREA
    AREA_sectors_CO2$STATE_NAME[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$STATE_NAME
    AREA_sectors_CO2$emis_sum2017[((areamyCODE_num-1)*150+1):(areamyCODE_num*150)] = EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area1 + EIA_CO2_emis_CONUS_days_j$CO2_emis_2017_perdy_area2
    
    areamyCODE_num = areamyCODE_num + 1
  }
}

#in POINT_sectors_CO2 separate PtIND_NG and PtIND_NG2 by CO emis ratio of the two
#in POINT_sectors_CO2 separate PtIND_Oil and PtIND_Oil2 by CO emis ratio of the two
POINT_sectors_CO_org <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Colin_SS_sumCO/sum_CO/POINT_sectors_CO_specificfuel.csv") #deleted all *_Wood #US in this file means US territories outside of CONUS

POINT_sectors_CO_IND_NG_org = POINT_sectors_CO_org%>% filter(sector == 'PtIND_NG')
POINT_sectors_CO_IND_NG2_org = POINT_sectors_CO_org%>% filter(sector == 'PtIND_NG2')
POINT_sectors_CO_IND_Oil_org = POINT_sectors_CO_org%>% filter(sector == 'PtIND_Oil')
POINT_sectors_CO_IND_Oil2_org = POINT_sectors_CO_org%>% filter(sector == 'PtIND_Oil2')

POINT_sectors_CO_IND_NG2_org = subset(POINT_sectors_CO_IND_NG2_org, select=c("day","STATE_FIPS","STATE_NAME","emis_sum"))
POINT_sectors_CO_IND_Oil2_org = subset(POINT_sectors_CO_IND_Oil2_org, select=c("day","STATE_FIPS","STATE_NAME","emis_sum"))

colnames(POINT_sectors_CO_IND_NG2_org)[colnames(POINT_sectors_CO_IND_NG2_org) == "emis_sum"] ="emis_sum2"
colnames(POINT_sectors_CO_IND_Oil2_org)[colnames(POINT_sectors_CO_IND_Oil2_org) == "emis_sum"] ="emis_sum2"

POINT_sectors_CO_IND_NG_both_org = join(POINT_sectors_CO_IND_NG_org,POINT_sectors_CO_IND_NG2_org, by=c("day","STATE_FIPS","STATE_NAME"), type="left")
POINT_sectors_CO_IND_Oil_both_org = join(POINT_sectors_CO_IND_Oil_org,POINT_sectors_CO_IND_Oil2_org, by=c("day","STATE_FIPS","STATE_NAME"), type="left")

POINT_sectors_CO_IND_NG_both_org['emis_sum_both'] <- as.numeric(POINT_sectors_CO_IND_NG_both_org$emis_sum) + as.numeric(POINT_sectors_CO_IND_NG_both_org$emis_sum2)
POINT_sectors_CO_IND_NG_both_org['frac_IND_NG'] <- as.numeric(POINT_sectors_CO_IND_NG_both_org$emis_sum) / as.numeric(POINT_sectors_CO_IND_NG_both_org$emis_sum_both)
POINT_sectors_CO_IND_NG_both_org['frac_IND_NG2'] <- as.numeric(POINT_sectors_CO_IND_NG_both_org$emis_sum2) / as.numeric(POINT_sectors_CO_IND_NG_both_org$emis_sum_both)

POINT_sectors_CO_IND_Oil_both_org['emis_sum_both'] <- as.numeric(POINT_sectors_CO_IND_Oil_both_org$emis_sum) + as.numeric(POINT_sectors_CO_IND_Oil_both_org$emis_sum2)
POINT_sectors_CO_IND_Oil_both_org['frac_IND_Oil'] <- as.numeric(POINT_sectors_CO_IND_Oil_both_org$emis_sum) / as.numeric(POINT_sectors_CO_IND_Oil_both_org$emis_sum_both)
POINT_sectors_CO_IND_Oil_both_org['frac_IND_Oil2'] <- as.numeric(POINT_sectors_CO_IND_Oil_both_org$emis_sum2) / as.numeric(POINT_sectors_CO_IND_Oil_both_org$emis_sum_both)

#if NaN for both frac_IND_NG and frac_IND_NG2 (caused by IND_NG + IND_NG2 emis = 0) make frac_IND_NG = 1 and frac_IND_NG2 = 0
for (r in  1:nrow(POINT_sectors_CO_IND_NG_both_org))
{
  if (is.nan(POINT_sectors_CO_IND_NG_both_org$frac_IND_NG[r]) & is.nan(POINT_sectors_CO_IND_NG_both_org$frac_IND_NG2[r]))
  {
    POINT_sectors_CO_IND_NG_both_org$frac_IND_NG[r] = 1
    POINT_sectors_CO_IND_NG_both_org$frac_IND_NG2[r] = 0
  }
}

#if NaN for both frac_IND_Oil and frac_IND_Oil2 (caused by IND_Oil + IND_Oil2 emis = 0) make frac_IND_Oil = 1 and frac_IND_Oil2 = 0
for (r in  1:nrow(POINT_sectors_CO_IND_Oil_both_org))
{
  if (is.nan(POINT_sectors_CO_IND_Oil_both_org$frac_IND_Oil[r]) & is.nan(POINT_sectors_CO_IND_Oil_both_org$frac_IND_Oil2[r]))
  {
    POINT_sectors_CO_IND_Oil_both_org$frac_IND_Oil[r] = 1
    POINT_sectors_CO_IND_Oil_both_org$frac_IND_Oil2[r] = 0
  }
}

#simplify POINT_sectors_CO_IND_NG_both_org
POINT_sectors_CO_IND_NG_2frac = subset(POINT_sectors_CO_IND_NG_both_org, select=c("sector", "day", "STATE_FIPS", "STATE_NAME", "frac_IND_NG", "frac_IND_NG2")) 

#simplify POINT_sectors_CO_IND_Oil_both_org
POINT_sectors_CO_IND_Oil_2frac = subset(POINT_sectors_CO_IND_Oil_both_org, select=c("sector", "day", "STATE_FIPS", "STATE_NAME", "frac_IND_Oil", "frac_IND_Oil2")) 

#separate IND_NG and IND_NG2
POINT_sectors_CO2 = join(POINT_sectors_CO2,POINT_sectors_CO_IND_NG_2frac, by=c("sector", "day","STATE_FIPS","STATE_NAME"), type="left")

POINT_sectors_CO2['emis_sum2017_IND_NG'] <- as.numeric(POINT_sectors_CO2$emis_sum2017) * as.numeric(POINT_sectors_CO2$frac_IND_NG)

POINT_sectors_CO2['emis_sum2017_IND_NG2'] <- as.numeric(POINT_sectors_CO2$emis_sum2017) * as.numeric(POINT_sectors_CO2$frac_IND_NG2)

#separate IND_Oil and IND_Oil2
POINT_sectors_CO2 = join(POINT_sectors_CO2,POINT_sectors_CO_IND_Oil_2frac, by=c("sector", "day","STATE_FIPS","STATE_NAME"), type="left")

POINT_sectors_CO2['emis_sum2017_IND_Oil'] <- as.numeric(POINT_sectors_CO2$emis_sum2017) * as.numeric(POINT_sectors_CO2$frac_IND_Oil)

POINT_sectors_CO2['emis_sum2017_IND_Oil2'] <- as.numeric(POINT_sectors_CO2$emis_sum2017) * as.numeric(POINT_sectors_CO2$frac_IND_Oil2)

#create a new POINT array
POINT_sectors_CO2_new <- data.frame(sector=character(num_states_CONUS*num_point_myCODEs_real*3), 
                                    day=character(num_states_CONUS*num_point_myCODEs_real*3), 
                                    STATE_FIPS=character(num_states_CONUS*num_point_myCODEs_real*3),
                                    STATE_NAME=character(num_states_CONUS*num_point_myCODEs_real*3), 
                                    emis_sum2017=character(num_states_CONUS*num_point_myCODEs_real*3))

PtCOMM_Coal_CO2 = POINT_sectors_CO2 %>% filter(sector == "PtCOMM_Coal")
PtCOMM_NG_CO2 = POINT_sectors_CO2 %>% filter(sector == "PtCOMM_NG")
PtCOMM_Oil_CO2 = POINT_sectors_CO2 %>% filter(sector == "PtCOMM_Oil")
PtIND_Coal_CO2 = POINT_sectors_CO2 %>% filter(sector == "PtIND_Coal")
PtIND_NG_CO2 = POINT_sectors_CO2 %>% filter(sector == "PtIND_NG")
PtIND_Oil_CO2 = POINT_sectors_CO2 %>% filter(sector == "PtIND_Oil")

loc_ind = 1
POINT_sectors_CO2_new$sector[((loc_ind-1)*150+1):(loc_ind*150)] = PtCOMM_Coal_CO2$sector
POINT_sectors_CO2_new$day[((loc_ind-1)*150+1):(loc_ind*150)] = PtCOMM_Coal_CO2$day
POINT_sectors_CO2_new$STATE_FIPS[((loc_ind-1)*150+1):(loc_ind*150)] = PtCOMM_Coal_CO2$STATE_FIPS
POINT_sectors_CO2_new$STATE_NAME[((loc_ind-1)*150+1):(loc_ind*150)] = PtCOMM_Coal_CO2$STATE_NAME
POINT_sectors_CO2_new$emis_sum2017[((loc_ind-1)*150+1):(loc_ind*150)] = PtCOMM_Coal_CO2$emis_sum2017

loc_ind = 2
POINT_sectors_CO2_new$sector[((loc_ind-1)*150+1):(loc_ind*150)] = PtCOMM_NG_CO2$sector
POINT_sectors_CO2_new$day[((loc_ind-1)*150+1):(loc_ind*150)] = PtCOMM_NG_CO2$day
POINT_sectors_CO2_new$STATE_FIPS[((loc_ind-1)*150+1):(loc_ind*150)] = PtCOMM_NG_CO2$STATE_FIPS
POINT_sectors_CO2_new$STATE_NAME[((loc_ind-1)*150+1):(loc_ind*150)] = PtCOMM_NG_CO2$STATE_NAME
POINT_sectors_CO2_new$emis_sum2017[((loc_ind-1)*150+1):(loc_ind*150)] = PtCOMM_NG_CO2$emis_sum2017

loc_ind = 3
POINT_sectors_CO2_new$sector[((loc_ind-1)*150+1):(loc_ind*150)] = PtCOMM_Oil_CO2$sector
POINT_sectors_CO2_new$day[((loc_ind-1)*150+1):(loc_ind*150)] = PtCOMM_Oil_CO2$day
POINT_sectors_CO2_new$STATE_FIPS[((loc_ind-1)*150+1):(loc_ind*150)] = PtCOMM_Oil_CO2$STATE_FIPS
POINT_sectors_CO2_new$STATE_NAME[((loc_ind-1)*150+1):(loc_ind*150)] = PtCOMM_Oil_CO2$STATE_NAME
POINT_sectors_CO2_new$emis_sum2017[((loc_ind-1)*150+1):(loc_ind*150)] = PtCOMM_Oil_CO2$emis_sum2017

loc_ind = 4
POINT_sectors_CO2_new$sector[((loc_ind-1)*150+1):(loc_ind*150)] = PtIND_Coal_CO2$sector
POINT_sectors_CO2_new$day[((loc_ind-1)*150+1):(loc_ind*150)] = PtIND_Coal_CO2$day
POINT_sectors_CO2_new$STATE_FIPS[((loc_ind-1)*150+1):(loc_ind*150)] = PtIND_Coal_CO2$STATE_FIPS
POINT_sectors_CO2_new$STATE_NAME[((loc_ind-1)*150+1):(loc_ind*150)] = PtIND_Coal_CO2$STATE_NAME
POINT_sectors_CO2_new$emis_sum2017[((loc_ind-1)*150+1):(loc_ind*150)] = PtIND_Coal_CO2$emis_sum2017

loc_ind = 5
POINT_sectors_CO2_new$sector[((loc_ind-1)*150+1):(loc_ind*150)] = 'PtIND_NG'
POINT_sectors_CO2_new$day[((loc_ind-1)*150+1):(loc_ind*150)] = PtIND_NG_CO2$day
POINT_sectors_CO2_new$STATE_FIPS[((loc_ind-1)*150+1):(loc_ind*150)] = PtIND_NG_CO2$STATE_FIPS
POINT_sectors_CO2_new$STATE_NAME[((loc_ind-1)*150+1):(loc_ind*150)] = PtIND_NG_CO2$STATE_NAME
POINT_sectors_CO2_new$emis_sum2017[((loc_ind-1)*150+1):(loc_ind*150)] = PtIND_NG_CO2$emis_sum2017_IND_NG

loc_ind = 6
POINT_sectors_CO2_new$sector[((loc_ind-1)*150+1):(loc_ind*150)] = 'PtIND_NG2'
POINT_sectors_CO2_new$day[((loc_ind-1)*150+1):(loc_ind*150)] = PtIND_NG_CO2$day
POINT_sectors_CO2_new$STATE_FIPS[((loc_ind-1)*150+1):(loc_ind*150)] = PtIND_NG_CO2$STATE_FIPS
POINT_sectors_CO2_new$STATE_NAME[((loc_ind-1)*150+1):(loc_ind*150)] = PtIND_NG_CO2$STATE_NAME
POINT_sectors_CO2_new$emis_sum2017[((loc_ind-1)*150+1):(loc_ind*150)] = PtIND_NG_CO2$emis_sum2017_IND_NG2

loc_ind = 7
POINT_sectors_CO2_new$sector[((loc_ind-1)*150+1):(loc_ind*150)] = 'PtIND_Oil'
POINT_sectors_CO2_new$day[((loc_ind-1)*150+1):(loc_ind*150)] = PtIND_Oil_CO2$day
POINT_sectors_CO2_new$STATE_FIPS[((loc_ind-1)*150+1):(loc_ind*150)] = PtIND_Oil_CO2$STATE_FIPS
POINT_sectors_CO2_new$STATE_NAME[((loc_ind-1)*150+1):(loc_ind*150)] = PtIND_Oil_CO2$STATE_NAME
POINT_sectors_CO2_new$emis_sum2017[((loc_ind-1)*150+1):(loc_ind*150)] = PtIND_Oil_CO2$emis_sum2017_IND_Oil

loc_ind = 8
POINT_sectors_CO2_new$sector[((loc_ind-1)*150+1):(loc_ind*150)] = 'PtIND_Oil2'
POINT_sectors_CO2_new$day[((loc_ind-1)*150+1):(loc_ind*150)] = PtIND_Oil_CO2$day
POINT_sectors_CO2_new$STATE_FIPS[((loc_ind-1)*150+1):(loc_ind*150)] = PtIND_Oil_CO2$STATE_FIPS
POINT_sectors_CO2_new$STATE_NAME[((loc_ind-1)*150+1):(loc_ind*150)] = PtIND_Oil_CO2$STATE_NAME
POINT_sectors_CO2_new$emis_sum2017[((loc_ind-1)*150+1):(loc_ind*150)] = PtIND_Oil_CO2$emis_sum2017_IND_Oil2

#write
write.csv(EIA_CO2_emis_res_States_sumMSN, file = "EIA_CO2_emis_res_States_sumMSN_v12.csv")
write.csv(EIA_CO2_emis_com_States_sumMSN, file = "EIA_CO2_emis_com_States_sumMSN_v12.csv")
write.csv(EIA_CO2_emis_ind_States_sumMSN, file = "EIA_CO2_emis_ind_States_sumMSN_v12.csv")
write.csv(EIA_CO2_emis_res_CONUS, file = "EIA_CO2_emis_res_CONUS_v12.csv")
write.csv(EIA_CO2_emis_com_CONUS, file = "EIA_CO2_emis_com_CONUS_v12.csv")
write.csv(EIA_CO2_emis_ind_CONUS, file = "EIA_CO2_emis_ind_CONUS_v12.csv")
write.csv(AREA_sectors_CO2, file = "AREA_sectors_CO2_v12.csv")
write.csv(POINT_sectors_CO2_new, file = "POINT_sectors_CO2_v12.csv")
write.csv(AREA_POINT_sumCO_what_need_patch, file = "AREA_POINT_sumCO_what_need_patch_v12.csv")


