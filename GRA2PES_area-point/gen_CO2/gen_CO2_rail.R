rm(list=ls()) #Clear memory
ptm <- proc.time()
#install.packages(data.table)
library(data.table)

options(stringsAsFactors=FALSE)

# Inputs
base_dir.area = '/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Colin_SS_sumCO/CO_spatial_surrogate/'
base_dir.point = '/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Colin_SS_sumCO/CO_spatial_surrogate/'
in_dir = '/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_ss_bug_fix/gen_CO2/'
out_dir = '/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_ss_bug_fix/gen_CO2/CO2_rds/'
emis.sector.area = c('RAIL')

emis.sector.point = c('PtRAIL')

fn_co2.area = 'AREA_rail_CO2.csv'
fn_co2.point = 'POINT_rail_CO2.csv'
var = 'emis_sum2017'
fn_base.area = '_CO_surr.rds'
fn_out_base.area = '_CO2.rds'
fn_base.point = '_CO_surr'
fn_out_base.point = '_CO2'
days =  c('weekdy','satdy','sundy')

colnames = c('dayav','HR01','HR02','HR03','HR04','HR05','HR06','HR07','HR08','HR09',
            'HR10','HR11','HR12','HR13','HR14','HR15','HR16','HR17','HR18','HR19','HR20',
            'HR21','HR22','HR23','HR24')

CO2.area = read.csv(file=paste(in_dir,fn_co2.area,sep=''),stringsAsFactors = FALSE)
CO2.point = read.csv(file=paste(in_dir,fn_co2.point,sep=''),stringsAsFactors = FALSE)

# Sum for area
for(i in 1:length(emis.sector.area))
{
  for(j in 1:length(days))
  {
    fn = paste(base_dir.area,days[j],'/',emis.sector.area[i],fn_base.area,sep='')
    print(fn)
    df = readRDS(fn)
    
    co2_value = CO2.area[CO2.area$sector==emis.sector.area[i] & CO2.area$day==days[j] ,var]

    for (col in 1:length(colnames))
    {
            
            df[,colnames[col]] = df[,colnames[col]]*co2_value

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
    fn = paste(base_dir.point,days[j],'/',emis.sector.point[i],fn_base.point,'_',days[j],'.rds',sep='')
    print(fn)
    df = readRDS(fn)

    co2_value = CO2.point[CO2.point$sector==emis.sector.point[i] & CO2.point$day==days[j] ,var]
    for (col in 1:length(colnames))
    {
        df[,colnames[col]] = df[,colnames[col]]*co2_value
    } 

    fn_out = paste(out_dir,days[j],'/',emis.sector.point[i],fn_out_base.point,'_',days[j],'.rds',sep='')
    print(fn_out)
    saveRDS(df,file = fn_out)
  }
}





