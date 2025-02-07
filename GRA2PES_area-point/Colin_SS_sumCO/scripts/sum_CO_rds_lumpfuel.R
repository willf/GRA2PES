rm(list=ls()) #Clear memory
ptm <- proc.time()
#install.packages(data.table)
library(data.table)

options(stringsAsFactors=FALSE)

# Inputs
base_dir.area = '/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Lump_CO_fuel/area/Month00/'
base_dir.point = '/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Lump_CO_fuel/point/Month00/'
out_dir = '/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Colin_SS_sumCO/sum_CO/'

emis.sector.area = c('RES_Allfuel','COMM_Allfuel','IND_Allfuel')
emis.sector.point = c('PtCOMM_Allfuel','PtIND_Allfuel')

out_fn.area = 'AREA_sectors_CO_lumpfuel.csv'
out_fn.point = 'POINT_sectors_CO_lumpfuel.csv'
fn_base = '_CO'
days =  c('weekdy','satdy','sundy')


out.area.df <- data.frame(sector=character(),
                    day=character(),
                    STATE_FIPS=integer(),
                    STATE_NAME=character(),
                    emis_sum=double(),
                    stringsAsFactors=FALSE)

out.point.df <- data.frame(sector=character(),
                          day=character(),
                          STATE_FIPS=integer(),
                          STATE_NAME=character(),
                          emis_sum=double(),
                          stringsAsFactors=FALSE)

# Sum for area
for(i in 1:length(emis.sector.area))
{
  for(j in 1:length(days))
  {
    fn = paste(base_dir.area,emis.sector.area[i],'/',days[j],'/',emis.sector.area[i],fn_base,'_',days[j],'.rds',sep='')
    print(fn)
    df = readRDS(fn)
    df$STATE_NAME = as.character(df$STATE_NAME)
    df$STATE_NAME[is.na(df$STATE_NAME)]="US"
    list = aggregate(df$dayav,by=list(STATE_FIPS=df$STATE_FIPS,STATE_NAME=df$STATE_NAME),FUN=sum)
    list$sector[1:length(list$STATE_FIPS)] = emis.sector.area[i]
    list$day[1:length(list$STATE_FIPS)] = days[j]
    
    out.area.df = rbind(out.area.df,data.frame(sector=list$sector,
                                     day=list$day,
                                     STATE_FIPS=list$STATE_FIPS,
                                     STATE_NAME=list$STATE_NAME,
                                     emis_sum=list$x,stringsAsFactors=FALSE))
    
  }
}

write.csv(out.area.df, file=paste(out_dir,out_fn.area,sep=''))

# Sum for point
for(i in 1:length(emis.sector.point))
{
  for(j in 1:length(days))
  {
    fn = paste(base_dir.point,emis.sector.point[i],'/',days[j],'/',emis.sector.point[i],fn_base,'_',days[j],'.rds',sep='')
    print(fn)
    df = readRDS(fn)
    df$STATE = as.character(df$STATE)
    df$STATE[is.na(df$STATE)]="US"
    df$STATE_FIPS = as.numeric(df$STATE_FIPS)
    df$STATE_FIPS[is.na(df$STATE_FIPS)]=0
    list = aggregate(df$dayav,by=list(STATE_FIPS=df$STATE_FIPS,STATE_NAME=df$STATE),FUN=sum)
    list$sector[1:length(list$STATE_FIPS)] = emis.sector.point[i]
    list$day[1:length(list$STATE_FIPS)] = days[j]
    
    out.point.df = rbind(out.point.df,data.frame(sector=list$sector,
                                     day=list$day,
                                     STATE_FIPS=list$STATE_FIPS,
                                     STATE_NAME=list$STATE_NAME,
                                     emis_sum=list$x,stringsAsFactors=FALSE))
    
  }
}

write.csv(out.point.df, file=paste(out_dir,out_fn.point,sep=''))




