#script to patch spatial surrogates

#install libraries
#install.packages("plyr")
library(plyr)
#install.packages("dplyr")
library(dplyr)
#install.packages("readr")
library(readr)

#set working directory
setwd ("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Patch_SS/CO_spatial_surrogates")
AREA_POINT_sumCO_what_need_patch <- read.csv("C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/emis_calc/output/AREA_POINT_sumCO_what_need_patch_v12.csv")
SS_base = 'C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Colin_SS_sumCO/CO_spatial_surrogate'
out_dir = 'C:/Users/clyu/Desktop/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Patch_SS/CO_spatial_surrogates/'

colnames = c('dayav','HR01','HR02','HR03','HR04','HR05','HR06','HR07','HR08','HR09',
             'HR10','HR11','HR12','HR13','HR14','HR15','HR16','HR17','HR18','HR19','HR20',
             'HR21','HR22','HR23','HR24')

#AREA
AREA_sumCO_what_need_patch = subset(AREA_POINT_sumCO_what_need_patch, select=c("AREAmyCODE","day","STATE_NAME"))
ss_base.area = '_CO_surr.rds'

#get uniq AREAmyCODE
AREAmyCODEs_uniq <- (AREA_sumCO_what_need_patch%>% distinct(AREAmyCODE, .keep_all=TRUE))$AREAmyCODE

for (cur_myCODE in AREAmyCODEs_uniq)
{
  print(cur_myCODE)
  cur_myCODE_what_need_patch = AREA_sumCO_what_need_patch %>% filter(AREAmyCODE == cur_myCODE)
  
  #get uniq day
  day_uniq <- (cur_myCODE_what_need_patch%>% distinct(day, .keep_all=TRUE))$day
  
  for (cur_dd in day_uniq)
  {
    print(cur_dd)
    cur_myCODE_dd_what_need_patch = cur_myCODE_what_need_patch %>% filter(day == cur_dd)
    
    FIRST_STATE_FLG1 = TRUE
    FIRST_STATE_FLG2 = TRUE
    for (rr in 1:nrow(cur_myCODE_dd_what_need_patch)) #all these rows use the same unpatched and patch file
    {
      #read cur_state
      cur_state = cur_myCODE_dd_what_need_patch$STATE_NAME[rr]
      print(cur_state)
      
      if (FIRST_STATE_FLG1)
      {
        #read unpatched ss file
        fn_unpatched = paste(SS_base,'/',cur_dd,'/',cur_myCODE,ss_base.area,sep='')
        print(fn_unpatched)
        ss_unpatched = readRDS(fn_unpatched)
        
        #decide what patch ss file to read and read it
        if (cur_myCODE %in% c('RES_Coal','RES_NG','RES_Oil'))
        {
          fn_patch = paste(SS_base,'/',cur_dd,'/RES_Allfuel',ss_base.area,sep='')
          print(fn_patch)
          ss_patch = readRDS(fn_patch)
        }
        else if (cur_myCODE %in% c('COMM_Coal','COMM_NG','COMM_Oil'))
        {
          fn_patch = paste(SS_base,'/',cur_dd,'/COMM_Allfuel',ss_base.area,sep='')
          print(fn_patch)
          ss_patch = readRDS(fn_patch)
        }
        else if (cur_myCODE %in% c('IND_Coal','IND_NG','IND_Oil'))
        {
          fn_patch = paste(SS_base,'/',cur_dd,'/IND_Allfuel',ss_base.area,sep='')
          print(fn_patch)
          ss_patch = readRDS(fn_patch)
        }
        FIRST_STATE_FLG1 = FALSE
      }
      
      #Get patch of the current state
      ss_patch_state = ss_patch %>% filter(STATE_NAME == cur_state)
      #state_ss_total = sum(ss_patch_state$dayav)
      #print("state_ss_total is")
      #print(state_ss_total) #1
      
      #rbind individual state patch together
      if (FIRST_STATE_FLG2)
      {
        ss_patch_allstates = ss_patch_state
        FIRST_STATE_FLG2 = FALSE
      }
      else
      {
        ss_patch_allstates = rbind(ss_patch_allstates,ss_patch_state)  
      }
    }  
    
    #rename ss_patch_allstates columns so it can be joined with unpatched ss
    for (colname in colnames)
    {
      colnames(ss_patch_allstates)[colnames(ss_patch_allstates) == colname] = paste(colname,"_patch",sep='')
    }
      
    #join patch with unpatched ss file
    ss_unpatched_join_patch = join(ss_unpatched,ss_patch_allstates, 
                              by=c("Id","Row","Col","LON","LAT","TZ","STATE_FIPS","STATE_NAME","URB_RUR","VCP_Reg","PADD","CZ"), type="left")
    
    #apply patch
    for (colname in colnames)
    {
      colname_patch = paste(colname,"_patch",sep='')
      colname_patched = paste(colname,"_patched",sep='')
      ss_unpatched_join_patch[colname_patched] = pmax(ss_unpatched_join_patch[colname],ss_unpatched_join_patch[colname_patch],na.rm = TRUE)
    }
    ss_patched = subset(ss_unpatched_join_patch,
                 select=c("Id","Row","Col","LON","LAT","TZ","STATE_FIPS","STATE_NAME","URB_RUR","VCP_Reg","PADD","CZ",
                 "dayav_patched","HR01_patched","HR02_patched","HR03_patched","HR04_patched","HR05_patched","HR06_patched","HR07_patched","HR08_patched","HR09_patched",
                 "HR10_patched","HR11_patched","HR12_patched","HR13_patched","HR14_patched","HR15_patched","HR16_patched","HR17_patched","HR18_patched","HR19_patched",
                 "HR20_patched","HR21_patched","HR22_patched","HR23_patched","HR24_patched"))
    
    #rename ss_patched columns so it can be written out
    for (colname in colnames)
    {
      colname_patched = paste(colname,"_patched",sep='')
      colnames(ss_patched)[colnames(ss_patched) == colname_patched] = colname
    }
    
    #write the patched ss file
    fn_out = paste(out_dir,cur_dd,'/',cur_myCODE,ss_base.area,sep='')
    print(fn_out)
    saveRDS(ss_patched,file = fn_out)
  }
}

#POINT
POINT_sumCO_what_need_patch = subset(AREA_POINT_sumCO_what_need_patch, select=c("POINTmyCODE","day","STATE_NAME"))
POINT_sumCO_what_need_patch = POINT_sumCO_what_need_patch %>% filter(POINTmyCODE!='NA')
ss_base.point = '_CO_surr_'

#get uniq POINTmyCODE
POINTmyCODEs_uniq <- (POINT_sumCO_what_need_patch%>% distinct(POINTmyCODE, .keep_all=TRUE))$POINTmyCODE

for (cur_myCODE in POINTmyCODEs_uniq)
{
  print(cur_myCODE)
  cur_myCODE_what_need_patch = POINT_sumCO_what_need_patch %>% filter(POINTmyCODE == cur_myCODE)
  
  #get uniq day
  day_uniq <- (cur_myCODE_what_need_patch%>% distinct(day, .keep_all=TRUE))$day
  
  for (cur_dd in day_uniq)
  {
    print(cur_dd)
    cur_myCODE_dd_what_need_patch = cur_myCODE_what_need_patch %>% filter(day == cur_dd)
    
    FIRST_STATE_FLG1 = TRUE
    FIRST_STATE_FLG2 = TRUE
    for (rr in 1:nrow(cur_myCODE_dd_what_need_patch)) #all these rows use the same unpatched and patch file
    {
      #read cur_state
      cur_state = cur_myCODE_dd_what_need_patch$STATE_NAME[rr]
      print(cur_state)
      
      if (FIRST_STATE_FLG1)
      {
        #read unpatched ss file
        fn_unpatched = paste(SS_base,'/',cur_dd,'/',cur_myCODE,ss_base.point,cur_dd,'.rds',sep='')
        print(fn_unpatched)
        ss_unpatched = readRDS(fn_unpatched)
        
        #decide what patch ss file to read and read it
        if (cur_myCODE %in% c('PtCOMM_Coal','PtCOMM_NG','PtCOMM_Oil'))
        {
          fn_patch = paste(SS_base,'/',cur_dd,'/PtCOMM_Allfuel',ss_base.point,cur_dd,'.rds',sep='')
          print(fn_patch)
          ss_patch = readRDS(fn_patch)
        }
        else if (cur_myCODE %in% c('PtIND_Coal','PtIND_NG','PtIND_NG2','PtIND_Oil','PtIND_Oil2'))
        {
          fn_patch = paste(SS_base,'/',cur_dd,'/PtIND_Allfuel',ss_base.point,cur_dd,'.rds',sep='')
          print(fn_patch)
          ss_patch = readRDS(fn_patch)
        }
        FIRST_STATE_FLG1 = FALSE
      }
      
      #add unique identifier to the points
      ss_unpatched$uniq_ID = seq(1, nrow(ss_unpatched), by=1)
      ss_patch$uniq_ID = seq(1, nrow(ss_patch), by=1)
      
      #Get patch of the current state
      ss_patch_state = ss_patch %>% filter(STATE == cur_state)
      #state_ss_total = sum(ss_patch_state$dayav)
      #print("state_ss_total is")
      #print(state_ss_total) #1
      
      #rbind individual state patch together
      if (FIRST_STATE_FLG2)
      {
        ss_patch_allstates = ss_patch_state
        FIRST_STATE_FLG2 = FALSE
      }
      else
      {
        ss_patch_allstates = rbind(ss_patch_allstates,ss_patch_state)  
      }
    }  
    
    #rename ss_patch_allstates columns so it can be joined with unpatched ss
    for (colname in colnames)
    {
      colnames(ss_patch_allstates)[colnames(ss_patch_allstates) == colname] = paste(colname,"_patch",sep='')
    }
    
    #join patch with unpatched ss file
    ss_unpatched_join_patch = join(ss_unpatched,ss_patch_allstates, 
                                   by=c("uniq_ID","LON","LAT","STATE","STATE_FIPS"), type="left")
    
    #apply patch
    for (colname in colnames)
    {
      colname_patch = paste(colname,"_patch",sep='')
      colname_patched = paste(colname,"_patched",sep='')
      ss_unpatched_join_patch[colname_patched] = pmax(ss_unpatched_join_patch[colname],ss_unpatched_join_patch[colname_patch],na.rm = TRUE)
    }
    ss_patched = subset(ss_unpatched_join_patch,
                 select=c("LON","LAT","STATE","STATE_FIPS",
                 "dayav_patched","HR01_patched","HR02_patched","HR03_patched","HR04_patched","HR05_patched","HR06_patched","HR07_patched","HR08_patched","HR09_patched",
                 "HR10_patched","HR11_patched","HR12_patched","HR13_patched","HR14_patched","HR15_patched","HR16_patched","HR17_patched","HR18_patched","HR19_patched",
                 "HR20_patched","HR21_patched","HR22_patched","HR23_patched","HR24_patched"))
    
    #rename ss_patched columns so it can be written out
    for (colname in colnames)
    {
      colname_patched = paste(colname,"_patched",sep='')
      colnames(ss_patched)[colnames(ss_patched) == colname_patched] = colname
    }
    
    #write the patched ss file
    fn_out = paste(out_dir,cur_dd,'/',cur_myCODE,ss_base.point,cur_dd,'.rds',sep='')
    print(fn_out)
    saveRDS(ss_patched,file = fn_out)
  }
}
