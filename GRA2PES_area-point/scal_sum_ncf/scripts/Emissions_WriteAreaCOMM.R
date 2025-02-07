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
month.list   <- c(1,2,3,4,5,6,7,8,9,10,11,12) ###CHANGE
dow.list     <- c("weekdy","satdy","sundy")

base.dir     <- paste("/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_GRA2PES2021/scal_sum_ncf/base2017_rds/area/Month00",sep="") ###CHANGE
#temporal scaling factors, these don't have the later months in the year, clyu rearranged some files from charkins original files to match new sectoral definitions 
input.dir    <- "/wrk/charkins/emissions/GRA2PES/V7_NRT_scaling/AREA21_202404/input" ###CHANGE
output.dir   <- "/wrk/users/charkins/emissions/V7_GRA2PES/AREA21" ###CHANGE

#Set ArcGIS shapefile with domain
domain.dir   <- "/wrk/csd4/bmcdonald/COVID_Modeling/Domains"
domain.name  <- "nei04k_domain"

#Define sectors following CAMS
comm.list    <- c("COMM_Coal","COMM_Oil","COMM_NG","COMM_Wood")
res.list     <- c("RES_Coal","RES_Oil","RES_NG","RES_Wood")
indf.list    <- c("IND_Coal","IND_Oil","IND_NG","IND_Wood")
indp.list    <- c("CHEM","FOOD","METAL","PULP","MACHINE","OTH","MISC")
fog.list     <- c("OnG") #use Colin's file, no need to scale, sector name changed from AreaOG to AreaFOG
ognofog.list <- c("OffG","STORAGE","EVAPGAS") #Oil and Gas noFOG: this is a new sector as we decided to define these as oil and gas instead of fugitives
fug.list     <- c("CONST","DUST") #three old fugitives now defined as oil and gas
rail.list    <- c("RAIL") #we are not taking CMV from NEI17 anymore because of using CAMS shipping, so MOB is changed to RAIL
vcp.list     <- c("IndCoat","Degreasing","Inks","IndAdhesive","AgPesticide","TotlVCP") #use Colin's file, no need to scale
ag.list      <- c("CROP","AgBURN","LIVESTOCK") #link base year, no need to scale
waste.list   <- c("WASTE") #link base year, no need to scale

sector.list  <- c(comm.list)

out.name     <- "AreaCOMM"

scale.file   <- paste(out.name,"_monthly.csv",sep="")

#Variables to process
crit.list         <- c("CO2","CO","NH3","NOX","PM10-PRI","PM25-PRI","SO2","VOC")
hc.list           <- c("HC01","HC02","HC03","HC04","HC05","HC06","HC07","HC08","HC09","HC10",
                       "HC11","HC12","HC13","HC14","HC15","HC16","HC17","HC18","HC19","HC20",
                       "HC21","HC22","HC23","HC24","HC25","HC26","HC27","HC28","HC29","HC30",
                       "HC31","HC32","HC33","HC34","HC35","HC36","HC37","HC38","HC39","HC40",
                       "HC41","HC42","HC43","HC44","HC45","HC46","HC47","HC48","HC49","HC50")
pm.list           <- c("PM01","PM02","PM03","PM04","PM05","PM06","PM07","PM08","PM09","PM10",
                       "PM11","PM12","PM13","PM14","PM15","PM16","PM17","PM18","PM19")

var.list          <- c(crit.list,hc.list,pm.list)

#Set flags to turn on functions
read.domain.flag <- FALSE
out.rds.flag     <- TRUE
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
      
      exist.count <- 0
      
      #Cycle by sector
      for(s in 1:length(sector.list))
      {
        
        sector <- sector.list[s]
        
        #---------------------------------------
        # Sum emissions data across sectors
        #---------------------------------------
        setwd(paste(base.dir,"/",sector,"/",dow,sep=""))
        
        emis.file <- grep(paste(sector,"_",pollutant,".rds",sep=""),list.files(),value = TRUE)
        
        if(length(emis.file) > 0)
        {
          
          print(paste("Processing ",pollutant," for ",sector," in month ",month," on ",dow,"...",sep=""))
          
        }else{
          
          print(paste("Missing ",pollutant," for ",sector," in month ",month," on ",dow,"...",sep=""))
          
        }
        
        if(length(emis.file) > 0)
        {
          
          input.df  <- readRDS(emis.file)
          
          coord.df  <- domain.sp[,c(1:6)]
          
          start.col <- ncol(input.df)-24
          end.col   <- ncol(input.df)
          
          emis.df   <- input.df[,c(start.col:end.col)] #Exclude coordinate data
          
          #Read monthly scale factors
          scale.data  <- read.csv(paste(input.dir,"/",scale.file,sep=""))
          scale.col   <- grep(paste("^",sector,"$",sep=""), colnames(scale.data))
          scale.value <- scale.data[month,scale.col]
          
          scale.factors <- rep(scale.value,nrow(emis.df))
          
          #Initialize summation array
          if(exist.count == 0)
          {
            
            sum.emis.df <- array(0, c(nrow(emis.df),ncol(emis.df)))
            sum.emis.df <- sum.emis.df + emis.df * scale.factors[s]
            
            exist.count <- exist.count + 1
            
          }else{
            
            sum.emis.df <- sum.emis.df + emis.df * scale.factors[s]
            
          }
          
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
      
    } #End pollutant for loop
    
  } #End day-of-week for loop
  
} #End month for loop

print(proc.time()-ptm)
