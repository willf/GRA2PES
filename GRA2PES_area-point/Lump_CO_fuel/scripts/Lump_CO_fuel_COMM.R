#Combine fuel types
library(plyr)

sector = 'COMM' #CHANGE
base_dir = '/wrk/d2/bmcdonald/NEI17/area/Month00/' #CHANGE
out_dir = '/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Lump_CO_fuel/area/Month00/' #CHANGE

colnames = c('dayav','HR01','HR02','HR03','HR04','HR05','HR06','HR07','HR08','HR09',
             'HR10','HR11','HR12','HR13','HR14','HR15','HR16','HR17','HR18','HR19','HR20',
             'HR21','HR22','HR23','HR24')

days =  c('weekdy','satdy','sundy')

for(j in 1:length(days))
{ 
  #read Month00 CO files, the order of grid cells are the same #CHANGE
  Coal_CO_fn = paste(base_dir,sector,'_Coal/',days[j],'/',sector,'_Coal_CO.rds',sep='')
  NG_CO_fn = paste(base_dir,sector,'_NG/',days[j],'/',sector,'_NG_CO.rds',sep='')
  Oil_CO_fn = paste(base_dir,sector,'_Oil/',days[j],'/',sector,'_Oil_CO.rds',sep='')

  #CHANGE
  Coal_CO = readRDS(Coal_CO_fn)
  NG_CO = readRDS(NG_CO_fn)
  Oil_CO = readRDS(Oil_CO_fn)

  #create combine array
  Allfuel_CO = subset(Coal_CO, select=c("Id","Row","Col","LON","LAT","TZ","STATE_FIPS","STATE_NAME","URB_RUR","VCP_Reg","PADD","CZ"))

  #combine fuels
  for (col in 1:length(colnames))
  { 
    Allfuel_CO[colnames[col]] = Coal_CO[colnames[col]] + NG_CO[colnames[col]] + Oil_CO[colnames[col]]#CHANGE
  }
  
  fn_out = paste(out_dir,sector,'_Allfuel/',days[j],'/',sector,'_Allfuel_CO_',days[j],'.rds',sep='')
  print(fn_out)
  saveRDS(Allfuel_CO,file = fn_out)
}

