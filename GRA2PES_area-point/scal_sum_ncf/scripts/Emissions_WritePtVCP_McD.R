rm(list=ls()) #Clear memory
ptm <- proc.time()
library(foreign)
#library(rgdal)
#library(sp)
#library(maptools)

#--------------------------------------------------------
# Set input variables
#--------------------------------------------------------
month.list        <- c(1,2,3,4,5,6,7,8,9,10,11,12) #August only #CHANGE

#Set observation files and directories
base.in     <- paste("/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_GRA2PES2021/scal_sum_ncf/base2017_rds/point/Month00",sep="") ###CHANGE
input.dir <- "/wrk/charkins/emissions/GRA2PES/V7_NRT_scaling/POINT21_202404/input" #temporal scaling factors ###CHANGE
output.dir   <- "/wrk/users/charkins/emissions/V7_GRA2PES/POINT21" ###CHANGE

#Day-of-week subdirectories
dow.list          <- c("weekdy","satdy","sundy")

#Sectors to process
out.name          <- "PtVCP"
sector.list       <- c("PtAdhesives","PtCoatings","PtDegreasing","PtInks")

vcp.factors.us    <- c(164*1.00, 2.9*1.00, 0.48, 0.76*1.00)
vcp.factors.ca    <- c(164*0.70, 2.9*0.65, 0.48, 0.76*1.00)
vcp.factors.ny    <- c(164*0.90, 2.9*0.80, 0.48, 0.76*1.00)
vcp.factors.nj    <- c(164*0.90, 2.9*0.80, 0.48, 0.76*1.00)
vcp.factors.ct    <- c(164*0.90, 2.9*0.80, 0.48, 0.76*1.00)

scale.file        <- paste(out.name,"_monthly.csv",sep="")

#Variables to process
crit.list         <- c("CO","NH3","NOX","PM10-PRI","PM25-PRI","SO2","VOC")
hc.list           <- c("HC01","HC02","HC03","HC04","HC05","HC06","HC07","HC08","HC09","HC10",
                       "HC11","HC12","HC13","HC14","HC15","HC16","HC17","HC18","HC19","HC20",
                       "HC21","HC22","HC23","HC24","HC25","HC26","HC27","HC28","HC29","HC30",
                       "HC31","HC32","HC33","HC34","HC35","HC36","HC37","HC38","HC39","HC40",
                       "HC41","HC42","HC43","HC44","HC45","HC46","HC47","HC48","HC49","HC50",
                       "HC51","HC52","HC53","HC54","HC55","HC56","HC57","HC58","HC59","HC60",
                       "HC61","HC62","HC63","HC64","HC65","HC66","HC67","HC68")
pm.list           <- c("PM01","PM02","PM03","PM04","PM05","PM06","PM07","PM08","PM09","PM10",
                       "PM11","PM12","PM13","PM14","PM15","PM16","PM17","PM18","PM19")

var.list          <- c(crit.list, pm.list)

#Speciation files (number of speciation profiles must match number of sectors)
voc.wt.file       <- "ptOVCP_wt.csv"
voc.mw.file       <- "ptOVCP_mw.csv"

#Set flags
speciate.flag  <- TRUE
out.rds.flag   <- TRUE
output.flag    <- TRUE

#--------------------------------------------------------
# Initialize Files/Variables
#--------------------------------------------------------
setwd(input.dir)

voc.wt.df <- read.csv(voc.wt.file, header = TRUE)
voc.mw.df <- read.csv(voc.mw.file, header = TRUE)

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
      
      #Cycle by sector
      for(s in 1:length(sector.list))
      {
        sector <- sector.list[s]
        
        #---------------------------------------
        # Sum emissions data across sectors
        #---------------------------------------
        setwd(paste(base.in,"/",sector,"/",dow,sep=""))
        
        print(paste("Processing ",pollutant," for ",sector," on ",dow,"...",sep=""))
        
        emis.file <- grep(paste(sector,"_",pollutant,"_",dow,".rds",sep=""),list.files(),value = TRUE)
        input.df  <- readRDS(emis.file)
        
        coord.df  <- input.df[,c(1:4)]
        emis.df   <- input.df[,c(5:ncol(input.df))] #Exclude coordinate data
        
        ca.indices <- grep("California",  coord.df$STATE)
        ny.indices <- grep("New York",    coord.df$STATE)
        nj.indices <- grep("New Jersey",  coord.df$STATE)
        ct.indices <- grep("Connecticut", coord.df$STATE)
        
        #Initialize summation array
        if(s == 1)
        {
          sum.emis.df <- array(0, c(nrow(emis.df),ncol(emis.df)))
        }
        
        #Scale by monthly activity data
        scale.data  <- read.csv(paste(input.dir,"/",scale.file,sep=""))
        scale.col   <- grep(paste("^",sector,"$",sep=""), colnames(scale.data))
        scale.value <- scale.data[month,scale.col]
        
        month.factors <- rep(scale.value,nrow(emis.df))
        
        #Scale VCP emissions by McDonald et al. (2018)
        if(pollutant == "VOC")
        {
          
          #Account for regional regulations
          vcp.factors <- rep(vcp.factors.us[s],nrow(emis.df))
          vcp.factors[ca.indices] <- vcp.factors.ca[s]
          vcp.factors[ny.indices] <- vcp.factors.ny[s]
          vcp.factors[nj.indices] <- vcp.factors.nj[s]
          vcp.factors[ct.indices] <- vcp.factors.ct[s]
          
          sum.emis.df <- sum.emis.df + emis.df * month.factors * vcp.factors
          
        }else{
          
          sum.emis.df <- sum.emis.df + emis.df * month.factors
          
        }
        
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
              sector.voc.df <- emis.df * month.factors * vcp.factors * (sector.voc.wt / 100) * 10^6 / sector.voc.mw
              
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
        
        #Output total mass emissions (in metric tons)
        print(paste("Writing sum of ", pollutant, " on ", dow, "...", sep=""))
        
        emis.out <- cbind(coord.df, sum.emis.df)
        
        #Create sector directory
        if(dir.exists(out.name))
        {
          
          setwd(out.name)
          
        }else{
          
          print(paste("Creating ",out.name," output directory...",sep=""))
          dir.create(out.name)
          setwd(out.name)
          
        }
        
        #Create dow directory
        if(dir.exists(dow))
        {
          
          setwd(dow)
          
        }else{
          
          print(paste("Creating ",dow," output directory...",sep=""))
          dir.create(dow)
          setwd(dow)
          
        }
        
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
        setwd(paste(month.dir,"/",out.name,"/",dow,sep=""))
        
        #Cycle through each speciated bin
        for(bin in 1:nrow(voc.wt.df))
        {
          voc.name <- voc.wt.df$Bin[bin]
          
          print(paste("Writing sum of ", voc.name, " on ", dow, "...", sep=""))
          
          voc.out <- cbind(coord.df, sum.voc.list[[bin]])
          
          if(out.rds.flag)
          {
            
            saveRDS(voc.out, file = paste(out.name,"_",voc.name,".rds",sep=""))
            
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
