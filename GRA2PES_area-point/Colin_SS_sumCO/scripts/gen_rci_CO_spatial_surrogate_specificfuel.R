rm(list=ls()) #Clear memory
ptm <- proc.time()
#install.packages(data.table)
library(data.table)

options(stringsAsFactors=FALSE)

# Inputs
base_dir.area = '/wrk/d2/bmcdonald/NEI17/area/Month00/'
base_dir.point = '/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Stu_RELPT/output_RELPT_emis/Month00/'
out_dir = '/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Colin_SS_sumCO/CO_spatial_surrogate/'

emis.sector.area = c('RES_Coal','RES_NG','RES_Oil','COMM_Coal','COMM_NG','COMM_Oil','IND_Coal','IND_NG','IND_Oil')
emis.sector.point = c('PtCOMM_Coal','PtCOMM_NG','PtCOMM_Oil','PtIND_Coal','PtIND_NG','PtIND_NG2','PtIND_Oil','PtIND_Oil2')

fn_base.area = '_CO.rds'
fn_out_base.area = '_CO_surr.rds'
fn_base.point = '_CO'
fn_out_base.point = '_CO_surr'

days =  c('weekdy','satdy','sundy')
colnames = c('dayav','HR01','HR02','HR03','HR04','HR05','HR06','HR07','HR08','HR09',
            'HR10','HR11','HR12','HR13','HR14','HR15','HR16','HR17','HR18','HR19','HR20',
            'HR21','HR22','HR23','HR24')

# Sum for area
first = TRUE
for(i in 1:length(emis.sector.area))
{
  for(j in 1:length(days))
  {
    fn = paste(base_dir.area,emis.sector.area[i],'/',days[j],'/',emis.sector.area[i],fn_base.area,sep='')
    print(fn)
    df = readRDS(fn)
    df$STATE_NAME = as.character(df$STATE_NAME)
    df$STATE_NAME[is.na(df$STATE_NAME)]="US"
    
    states = unique(df$STATE_NAME)

    list = aggregate(df[colnames],by=list(STATE_FIPS=df$STATE_FIPS,STATE_NAME=df$STATE_NAME),FUN=sum)
    
    for (s in 1:length(states))
    {
        print(states[s])
        for (col in 1:length(colnames))
        {
            if (list[list$STATE_NAME==states[s],'dayav'] > 0 )
            {
                df[df$STATE_NAME==states[s],colnames[col]] = df[df$STATE_NAME==states[s],colnames[col]]/list[list$STATE_NAME==states[s],'dayav']
            }
            else
            {
                df[df$STATE_NAME==states[s],colnames[col]] = 0 
            }
        }
        
    } 
        
    
    
    
    fn_out = paste(out_dir,days[j],'/',emis.sector.area[i],fn_out_base.area,sep='')
    print(fn_out)
    saveRDS(df,file = fn_out)
  }
}



# Sum for point
for(i in 1:length(emis.sector.point))
{
  for(j in 1:length(days))
  {
    fn = paste(base_dir.point,emis.sector.point[i],'/',days[j],'/',emis.sector.point[i],fn_base.point,'_',days[j],'.rds',sep='')
    print(fn)
    df = readRDS(fn)
    df$STATE = as.character(df$STATE)
    df$STATE[is.na(df$STATE)]="US"
    df$STATE_FIPS = as.numeric(df$STATE_FIPS)
    df$STATE_FIPS[is.na(df$STATE_FIPS)]=0
    list = aggregate(df[colnames],by=list(STATE_FIPS=df$STATE_FIPS,STATE_NAME=df$STATE),FUN=sum)

    
    states = unique(df$STATE)
    

    for (s in 1:length(states))
    {
        print(states[s])
        for (col in 1:length(colnames))
        {
            if (list[list$STATE_NAME==states[s],'dayav'] > 0 )
            {
                df[df$STATE==states[s],colnames[col]] = df[df$STATE==states[s],colnames[col]]/list[list$STATE_NAME==states[s],'dayav']
            }
            else
            {
                df[df$STATE==states[s],colnames[col]] = 0 
            }
        }

    } 

    fn_out = paste(out_dir,days[j],'/',emis.sector.point[i],fn_out_base.point,'_',days[j],'.rds',sep='')
    print(fn_out)
    saveRDS(df,file = fn_out)
  }
}





