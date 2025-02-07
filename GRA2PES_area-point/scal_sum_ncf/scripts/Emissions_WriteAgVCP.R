rm(list=ls()) #Clear memory
ptm <- proc.time()
library(foreign)
#library(rgdal)
#library(sp)
#library(maptools)

options(warn = 2)

#--------------------------------------------------------
# Set input variables
#--------------------------------------------------------
#Set observation files and directories
#base.dir      <- "//ozone10g/csd4_bmcdonald/NYC_Modeling/Emissions/n11_bySector/area4k"
#input.dir     <- "//ozone10g/csd4_bmcdonald/NYC_Modeling/Emissions/vcp18_agVOC/input"
#output.dir    <- "//ozone10g/csd4_bmcdonald/NYC_Modeling/Emissions/vcp18_agVOC"
#base.dir       <- "/wrk/csd4/bmcdonald/NYC_Modeling/Emissions/n11_bySector/area4k"
#input.dir      <- "/wrk/csd4/bmcdonald/NYC_Modeling/Emissions/vcp18_agVOC_v3_alk/input"
#output.dir     <- "/wrk/csd4/bmcdonald/NYC_Modeling/Emissions/vcp18_agVOC_v3_alk"

#Select month and day-of-week
month.list   <- c(1,2,3,4,5,6,7,8,9,10,11,12)
dow.list     <- c("weekdy","satdy","sundy")

base.dir     <- paste("/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_GRA2PES2021/scal_sum_ncf/base2017_rds/area/Month00",sep="")
input.dir <- "/wrk/charkins/emissions/GRA2PES/V7_NRT_scaling/VCP21_202404/input"
output.dir   <- "/wrk/users/charkins/emissions/V7_GRA2PES/VCP21"




#Set ArcGIS shapefile with domain
domain.dir   <- "/wrk/csd4/bmcdonald/COVID_Modeling/Domains"
domain.name  <- "nei04k_domain"

#Sectors to process
sector.list    <- c("AgPesticide")
scale.factors  <- c(4.5)

out.name = "AgVCP"

#Variables to process
var.list       <- c("VOC")

#Speciation files (number of speciation profiles must match number of sectors)
voc.wt.file    <- "agVCP_wt.csv"
voc.mw.file    <- "agVCP_mw.csv"

#Set flags
read.domain.flag <- FALSE
speciate.flag    <- TRUE
output.flag      <- TRUE

#--------------------------------------------------------
# Read Domain File
#--------------------------------------------------------
setwd(domain.dir)

if(read.domain.flag)
{
  print("Reading domain file...")
  
  #Read domain polygon shapefile
  domain.sp  <- read.dbf(paste(domain.name,".dbf",sep=""))
  
  #Save shapefile as .rds file (faster for future)
  saveRDS(domain.sp, file = paste(domain.name,".rds",sep=""))
  
}else{
  
  #Read domain polygon shapefile
  domain.sp <- readRDS(paste(domain.name,".rds",sep=""))
  
}

#--------------------------------------------------------
# Read Speciation Profiles
#--------------------------------------------------------
setwd(input.dir)

voc.wt.df <- read.csv(voc.wt.file, header = TRUE)
voc.mw.df <- read.csv(voc.mw.file, header = TRUE)

#--------------------------------------------------------
# Process Emissions
#--------------------------------------------------------
#Cycle by month
for(m in 1:length(month.list))
{
  
  month <- month.list[m]
  #Cycle by day-of-week
  for(d in 1:length(dow.list))
  {
    dow <- dow.list[d]
  
    #Cycle by variable
    for(v in 1:length(var.list))
    {
      pollutant <- var.list[v]
      
      #Cycle by sector
      for(s in 1:length(sector.list))
      {
  
        sector <- sector.list[s]
      
        #---------------------------------------
        # Sum emissions data across sectors
        #---------------------------------------
        setwd(paste(base.dir,'/',sector,"/",dow,sep=""))
        
        print(paste("Processing ",pollutant," for ",sector," on ",dow,"...",sep=""))
        
        emis.file <- grep(paste(sector,"_",pollutant,".rds",sep=""),list.files(),value = TRUE)
        input.df  <- readRDS(emis.file)
        
        coord.df  <- domain.sp[,c(1:6)]
        
        hr.cols_emis <- grep("HR", colnames(input.df))
        day.col_emis <- grep("dayav", colnames(input.df))
        
        emis.df   <- input.df[,c(day.col_emis , hr.cols_emis)] #Exclude coordinate data
        
        #Initialize summation array
        if(s == 1)
        {
          sum.emis.df <- array(0, c(nrow(emis.df),ncol(emis.df)))
        }
        
        sum.emis.df <- sum.emis.df + emis.df * scale.factors[s]
        
        #---------------------------------------
        # Speciate emissions
        #---------------------------------------
        if(speciate.flag & pollutant == "VOC"){
          
          #Initialize VOC speciation array
  
          #Cycle through each hydrocarbon bin
          for(bin in 1:nrow(voc.wt.df))
          {
            
            bin.name <- voc.wt.df$Bin[bin]
            
            print(paste("Processing ",bin.name," for ",sector," on ",dow,"...",sep=""))
            
            #Obtain VOC weight fraction (in percent) and MW
            sector.col    <- grep(sector, colnames(voc.wt.df))
            sector.voc.wt <- voc.wt.df[bin,sector.col]
            sector.voc.mw <- voc.mw.df[bin,sector.col]
            
            #Estimate VOC emissions (in moles)
            if(is.na(sector.voc.mw))
            {
              #Create blank array if MW does not exist
              sector.voc.df <- data.frame(matrix(0, nrow=nrow(emis.df), ncol=ncol(emis.df)))
              colnames(sector.voc.df) <- colnames(emis.df)
              
            }else{
              
              #Save molar emissions if MW exists
              sector.voc.df <- emis.df * scale.factors[s] * (sector.voc.wt / 100) * 10^6 / sector.voc.mw
              
            }
            
            #Store sum of emissions in array
            if(bin == 1 & s == 1)
            {
              
              sum.voc.list <- list(sector.voc.df)
  
            }else if(bin > 1 & s == 1){
              
              sum.voc.list[[bin]] <- sector.voc.df
  
            }else{
              
              sum.voc.list[[bin]] <- sum.voc.list[[bin]] + sector.voc.df
  
            }
            
          } #End hydrocarbon bin for loop
          
        } #End speciation flag if statement
  
      } #End sector loop
  
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
        
        emis.out <- cbind(coord.df, sum.emis.df)
        
        saveRDS(emis.out, file = paste("AgVCP_",pollutant,".rds",sep=""))
        # write.csv(emis.out, file = paste("AgVCP_",pollutant,".csv",sep=""),
        #             row.names = FALSE)
        
      } #End output flag if statement
  
      #---------------------------------------
      # Write Speciated VOC Emissions Output
      #---------------------------------------    
      #Output speciated VOC emissions (in moles)
      if(output.flag & speciate.flag & pollutant == "VOC")
      {
        setwd(dow.dir)
        
        #Cycle through each speciated bin
        for(bin in 1:nrow(voc.wt.df))
        {
          voc.name <- voc.wt.df$Bin[bin]
          
          print(paste("Writing sum of ", voc.name, " on ", dow, "...", sep=""))
          
          voc.out <- cbind(coord.df, sum.voc.list[[bin]])
          
          saveRDS(voc.out, file = paste("AgVCP_", voc.name, ".rds",sep=""))
          # write.csv(voc.out, file = paste("AgVCP_", voc.name, ".csv",sep=""),
          #             row.names = FALSE)
          
        } #End output for loop (speciation)
        
      } #End output flag if statement (speciation)
      
    } #End pollutant for loop
  
  } #End day-of-week for loop

} #End month for loop
print(proc.time()-ptm)
