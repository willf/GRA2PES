rm(list=ls()) #Clear memory
ptm <- proc.time()
library(foreign)
#library(rgdal)
#library(sp)
#library(maptools)

#--------------------------------------------------------
# Set input variables
#--------------------------------------------------------
#Select month and day-of-week
dow.list     <- c("weekdy","satdy","sundy")
month.list   <- c(8)

base.dir   <- "/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_GRA2PES2021/scal_sum_ncf/scaled2021_rds/area"
output.dir <- "/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_GRA2PES2021/scal_sum_ncf/summed2021_rds/area"

#Define sectors following CAMS
sector.list  <- c("AreaINDF")
out.name     <- "TotlArea_INDF"

#Variables to process
crit.list    <- c("CO2","CO","NH3","NOX","PM10-PRI","PM25-PRI","SO2","VOC")
hc.list      <- c("HC01","HC02","HC03","HC04","HC05","HC06","HC07","HC08","HC09","HC10",
                  "HC11","HC12","HC13","HC14","HC15","HC16","HC17","HC18","HC19","HC20",
                  "HC21","HC22","HC23","HC24","HC25","HC26","HC27","HC28","HC29","HC30",
                  "HC31","HC32","HC33","HC34","HC35","HC36","HC37","HC38","HC39","HC40",
                  "HC41","HC42","HC43","HC44","HC45","HC46","HC47","HC48","HC49","HC50",
                  "HC51","HC52","HC53","HC54","HC55","HC56","HC57","HC58","HC59","HC60",
                  "HC61","HC62","HC63","HC64","HC65","HC66","HC67","HC68")
pm.list      <- c("PM01","PM02","PM03","PM04","PM05","PM06","PM07","PM08","PM09","PM10",
                  "PM11","PM12","PM13","PM14","PM15","PM16","PM17","PM18","PM19")

var.list     <- c(crit.list, hc.list, pm.list)

#Set flags to turn on functions
out.rds.flag <- TRUE
output.flag  <- TRUE

#--------------------------------------------------------
# Summarize Emissions
#--------------------------------------------------------
#Cycle by month
for(m in 1:length(month.list))
{
  
  month <- month.list[m]
  
  #Cycle by day-of-week
  for(d in 1:length(dow.list))
  {
    
    dow <- dow.list[d]
    
    #Cycle by pollutant
    for(v in 1:length(var.list))
    {
      
      pollutant <- var.list[v]
      
      exist.count <- 0
      
      #Cycle by sector
      for(s in 1:length(sector.list))
      {
        
        sector <- sector.list[s]
        
        #---------------------------------------
        # Sum daily emissions
        #---------------------------------------
        setwd(paste(base.dir,"/Month",sprintf("%02d",month),"/",sector,"/",dow,sep=""))
        
        print(paste("Reading ",pollutant," for ",sector," on month ",month," ",dow,"...",sep=""))
        
        emis.file <- grep(paste(sector,"_",pollutant,".rds",sep=""),list.files(),value = TRUE)
        
        if(length(emis.file) > 0)
        {
          
          emis.df <- readRDS(emis.file)

          if(exist.count == 0)
          {
            
            sum.df <- emis.df
            
          }else{
            
            hr.cols <- grep("HR", colnames(sum.df))
            day.col <- grep("dayav", colnames(sum.df))
            
            sum.df[,c(day.col, hr.cols)] <- sum.df[,c(day.col, hr.cols)] + emis.df[,c(day.col, hr.cols)]
            
          }
          
          exist.count <- exist.count + 1
          
        } else if (length(emis.file) == 0 & exist.count == 0  ) {
          ex.file <- grep(paste(sector,"_","CO",".rds",sep=""),list.files(),value = TRUE)
          emis.df  = readRDS(ex.file)
          sum.df <- emis.df
          hr.cols <- grep("HR", colnames(sum.df))
          day.col <- grep("dayav", colnames(sum.df))
          sum.df[,c(day.col, hr.cols)] <- sum.df[,c(day.col, hr.cols)]*0 #End if file exists
        }
      
      }#End sector for loop
      
      #---------------------------------------
      # Write Emissions Output
      #---------------------------------------
      month.dir <- paste(output.dir,"/","Month",sprintf("%02d",month),sep="")
      
      if(output.flag)
      {
        setwd(month.dir)
        
        #Create sector directory
        sector.dir <- paste(month.dir,"/",out.name,sep="")
        
        if(dir.exists(sector.dir))
        {
          
          setwd(sector.dir)
          
        }else{
          
          print(paste("Creating ",out.name," output directory...",sep=""))
          dir.create(sector.dir)
          setwd(sector.dir)
          
        }
        
        #Create dow directory
        dow.dir <- paste(sector.dir,"/",dow,sep="")
        
        if(dir.exists(dow.dir))
        {
          
          setwd(dow.dir)
          
        }else{
          
          print(paste("Creating ",dow," output directory...",sep=""))
          dir.create(dow.dir)
          setwd(dow.dir)
          
        }
        
        #Output total mass emissions (in metric tons)
        print(paste("Writing sum of ", pollutant, " on ", dow, "...", sep=""))
        
        if(out.rds.flag)
        {
          
          saveRDS(sum.df, file = paste(out.name,"_",pollutant,".rds",sep=""))
          
        }else{
          
          write.csv(sum.df, file = paste(out.name,"_",pollutant,".csv",sep=""),
                    row.names = FALSE)
          
        }
        
      } #End output flag if statement
      
    } #End pollutant for loop
    
  } #End day-of-week loop
  
} #End month for loop

print(proc.time()-ptm)
