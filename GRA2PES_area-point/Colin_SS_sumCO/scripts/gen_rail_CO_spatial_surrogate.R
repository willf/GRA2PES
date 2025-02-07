rm(list=ls()) #Clear memory
ptm <- proc.time()
#install.packages(data.table)
library(data.table)

options(stringsAsFactors=FALSE)

# Inputs
base_dir.area = '/wrk/d2/bmcdonald/NEI17/area/Month00/'
base_dir.point = '/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Stu_RELPT/output_RELPT_emis/Month00/'
out_dir = '/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_ss_bug_fix/Colin_SS_sumCO/CO_spatial_surrogate/'
emis.sector.area = c('RAIL')

emis.sector.point = c('PtRAIL')
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
    sum = sum(df$dayav)
        for (col in 1:length(colnames))
        {
            df[,colnames[col]] = df[,colnames[col]]/sum
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
    sum = sum(df$dayav)

        for (col in 1:length(colnames))
        {
            df[,colnames[col]] = df[,colnames[col]]/sum
        }

    fn_out = paste(out_dir,days[j],'/',emis.sector.point[i],fn_out_base.point,'_',days[j],'.rds',sep='')
    print(fn_out)
    saveRDS(df,file = fn_out)
  }
}





