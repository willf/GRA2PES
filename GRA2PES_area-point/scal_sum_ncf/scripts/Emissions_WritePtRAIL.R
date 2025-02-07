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
month.list        <- c(1,2,3,4,5,6,7,8,9,10,11,12) ###CHANGE
dow.list          <- c("weekdy","satdy","sundy")

#Set observation files and directories
base.in     <- paste("/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_GRA2PES2021/scal_sum_ncf/base2017_rds/point/Month00",sep="") ###CHANGE
#temporal scaling factors, these don't have the later months in the year, clyu rearranged some files from charkins original files to match new sectoral definitions 
input.dir <- "/wrk/charkins/emissions/GRA2PES/V7_NRT_scaling/POINT21_202404/input" ###CHANGE
output.dir   <- "/wrk/users/charkins/emissions/V7_GRA2PES/POINT21" ###CHANGE

#Sectors to process
egu.list     <- c("PtEGU_Coal","PtEGU_NG","PtEGU_Oil","PtEGU_BIO")
comm.list    <- c("PtCOMM_Coal","PtCOMM_NG","PtCOMM_Oil","PtCOMM_BIO")
indf.list    <- c("PtIND_Coal","PtIND_NG","PtIND_NG2","PtIND_Oil","PtIND_Oil2","PtIND_BIO")
indp.list    <- c("PtCHEM","PtFOOD","PtMETAL","PtREFINE","PtPULP","PtELECT","PtMOTOR","PtAPPAREL","PtPHOTO","PtDRUG","PtMISC","PtMISC2")
vcp.list     <- c("PtAdhesives","PtCoatings","PtDegreasing","PtInks")
og.list      <- c("PtOnG","PtEVAPGAS","PtSTORAGE") #sectoral definition updated
fug.list     <- c("PtCONST") #sectoral definition updated
rail.list    <- c("PtRAIL") #separate out from MOB
air.list     <- c("PtAVIATION") #separate out from MOB

sector.list  <- c(rail.list)

out.name          <- "PtRAIL"

scale.file    <- paste(out.name,"_monthly.csv",sep="")
cems.CO2.file <- paste("PtEGU_cemsCO2.csv",sep="")
cems.NOX.file <- paste("PtEGU_cemsNOX.csv",sep="")
cems.SO2.file <- paste("PtEGU_cemsSO2.csv",sep="")

#Variables to process
crit.list         <- c("CO2","CO","NH3","NOX","PM10-PRI","PM25-PRI","SO2","VOC")
hc.list           <- c("HC01","HC02","HC03","HC04","HC05","HC06","HC07","HC08","HC09","HC10",
                       "HC11","HC12","HC13","HC14","HC15","HC16","HC17","HC18","HC19","HC20",
                       "HC21","HC22","HC23","HC24","HC25","HC26","HC27","HC28","HC29","HC30",
                       "HC31","HC32","HC33","HC34","HC35","HC36","HC37","HC38","HC39","HC40",
                       "HC41","HC42","HC43","HC44","HC45","HC46","HC47","HC48","HC49","HC50")
pm.list           <- c("PM01","PM02","PM03","PM04","PM05","PM06","PM07","PM08","PM09","PM10",
                       "PM11","PM12","PM13","PM14","PM15","PM16","PM17","PM18","PM19")

var.list          <- c(crit.list, hc.list, pm.list)

#Set flags
speciate.flag  <- FALSE
cems.flag      <- FALSE
out.rds.flag   <- TRUE
output.flag    <- TRUE

#--------------------------------------------------------
# Process Emissions
#--------------------------------------------------------
setwd(base.in)

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
      
      exist.count <- 0
      
      #Cycle by sector
      for(s in 1:length(sector.list))
      {
        sector <- sector.list[s]
        
        #---------------------------------------
        # Sum emissions data across sectors
        #---------------------------------------
        setwd(paste(base.in,"/",sector,"/",dow,sep=""))
        
        emis.file <- grep(paste(sector,"_",pollutant,"_",dow,".rds",sep=""),list.files(),value = TRUE)
        
        if(length(emis.file) > 0)
        {
          
          print(paste("Processing ",pollutant," for ",sector," on ",dow,"...",sep=""))
          
        }else{
          
          print(paste("Missing ",pollutant," for ",sector," on ",dow,"...",sep=""))
          
        }
        
        if(length(emis.file) > 0)
        {
          
          input.df  <- readRDS(emis.file)
          
          coord.df  <- input.df[,c(1:4)]
          emis.df   <- input.df[,c(5:ncol(input.df))] #Exclude coordinate data
          
          #Initialize summation array
          if(exist.count == 0)
          {
            sum.emis.df <- array(0, c(nrow(emis.df),ncol(emis.df)))
            exist.count <- exist.count + 1
          }
          
          #scale EGU CO2 SO2 NOX by relpoint-specific scaling factors generated from monthly and annual CEMS emissions relpoint templates
          if(sector == "PtEGU_Coal" | sector == "PtEGU_NG" | sector == "PtEGU_Oil" | sector == "PtEGU_BIO")
          {
            
            if(pollutant == "NOX" & cems.flag)
            {
              
              #Scale factors provided by point
              scale.data  <- read.csv(paste(input.dir,"/",cems.NOX.file,sep=""))
              scale.col   <- grep(paste("Month",sprintf("%02d",month),sep=""), colnames(scale.data))
              scale.factors <- scale.data[,scale.col]
              
            }else if(pollutant == "SO2" & cems.flag){
              
              #Scale factors provided by point
              scale.data  <- read.csv(paste(input.dir,"/",cems.SO2.file,sep=""))
              scale.col   <- grep(paste("Month",sprintf("%02d",month),sep=""), colnames(scale.data))
              scale.factors <- scale.data[,scale.col]
            
            }else if(pollutant == "CO2" & cems.flag){
              
              #Scale factors provided by point
              scale.data  <- read.csv(paste(input.dir,"/",cems.CO2.file,sep=""))
              scale.col   <- grep(paste("Month",sprintf("%02d",month),sep=""), colnames(scale.data))
              scale.factors <- scale.data[,scale.col]
              
              
            }else{
              
              #Otherwise scale nationally by fuel type
              scale.data  <- read.csv(paste(input.dir,"/",scale.file,sep=""))
              scale.col   <- grep(paste("^",sector,"$",sep=""), colnames(scale.data))
              scale.value <- scale.data[month,scale.col]
              
              scale.factors <- rep(scale.value,nrow(emis.df))
              
            }
            
          }else{
            
            #Otherwise scale nationally by fuel type
            scale.data  <- read.csv(paste(input.dir,"/",scale.file,sep=""))
            scale.col   <- grep(paste("^",sector,"$",sep=""), colnames(scale.data))
            scale.value <- scale.data[month,scale.col]
            
            scale.factors <- rep(scale.value,nrow(emis.df))
            
          }
          
          sum.emis.df <- sum.emis.df + emis.df * scale.factors
          
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
                sector.voc.df <- emis.df * scale.factors * (sector.voc.wt / 100) * 10^6 / sector.voc.mw
                
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
          
        } #End if file exists
        
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
        
        if(out.rds.flag)
        {
          
          saveRDS(emis.out, file = paste(out.name,"_",pollutant,".rds",sep=""))
          
        }else{
          
          write.csv(emis.out, file = paste(out.name,"_",pollutant,".csv",sep=""),
                    row.names = FALSE)
          
        }
        
      } #End output flag if statement
      
      #---------------------------------------
      # Write Speciated VOC Emissions Output
      #---------------------------------------    
      #Output speciated VOC emissions (in moles)
      if(output.flag & speciate.flag & pollutant == "VOC")
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
        
        #Cycle through each speciated bin
        for(bin in 1:nrow(voc.wt.df))
        {
          voc.name <- voc.wt.df$Bin[bin]
          
          print(paste("Writing sum of ", voc.name, " on ", dow, "...", sep=""))
          
          voc.out <- cbind(coord.df, sum.voc.list[[bin]])
          
          if(out.rds.flag)
          {
            
            saveRDS(voc.out, file = paste(out.name,"_",voc.name,".rds",sep=""),
                    row.names = FALSE)
            
          }else{
            
            write.csv(voc.out, file = paste(out.name,"_",voc.name,".csv",sep=""),
                      row.names = FALSE)
            
          }
          
        } #End output for loop (speciation)
        
      } #End output flag if statement (speciation)
      
    } #End pollutant for loop
    
  } #End day-of-week for loop 
  
} #End month for loop

print(proc.time()-ptm)
