#--------------------------------------------------------
# bash
# eval "$(/home/bmcdonald/anaconda3/bin/conda shell.bash hook)"
# conda activate r-gdal
#--------------------------------------------------------
rm(list=ls()) #Clear memory
ptm <- proc.time()
library(foreign, lib.loc="//home/bmcdonald/R/x86_64-redhat-linux-gnu-library/3.6")
#library(rgdal)
#library(sp)
#library(maptools)

#--------------------------------------------------------
# Set input variables
#--------------------------------------------------------
#Input directories and files
base.dir     <- "/wrk/d2/bmcdonald/NEI17/area/Month01"
domain.file  <- "/wrk/csd4/bmcdonald/COVID_Modeling/Domains/nei04k_domain.rds"
#base.dir    <- "//ozone10g/csd4_bmcdonald/NYC_Modeling/Emissions/n11_bySector"
#domain.file <- "//ozone10g/csd4_bmcdonald/NYC_Modeling/Domains/us04k_domain.dbf"

#Day-of-week subdirectories
dow.list     <- c("weekdy","satdy","sundy")

#Sectors to process
comm.list    <- c("COMM_Coal","COMM_Oil","COMM_NG","COMM_Wood")
res.list     <- c("RES_Coal","RES_Oil","RES_NG","RES_Wood")
indf.list    <- c("IND_Coal","IND_Oil","IND_NG","IND_Wood")
indp.list    <- c("CHEM","FOOD","METAL","PULP","MACHINE","OTH","MISC")
fug.list     <- c("OnG","OffG","EVAPGAS","STORAGE","CONST","DUST")
mob.list     <- c("CMV_DSL","CMV_OGV","RAIL")
vcp.list     <- c("IndCoat","Degreasing","Inks","IndAdhesive","AgPesticide","TotlVCP")
ag.list      <- c("CROP","AgBURN","LIVESTOCK")
waste.list   <- c("WASTE")

sector.list  <- c(comm.list)

#Variables to process
var.list    <- c("CO",    "NH3",  "NOX",  "SO2",  "VOC",
                 "HC01", "HC02", "HC03", "HC04", "HC05", "HC06", "HC07", "HC08", "HC09", "HC10",
                 "HC11", "HC12", "HC13", "HC14", "HC15", "HC16", "HC17", "HC18", "HC19", "HC20",
                 "HC21", "HC22", "HC23", "HC24", "HC25", "HC26", "HC27", "HC28", "HC29", "HC30",
                 "HC31", "HC32", "HC33", "HC34", "HC35", "HC36", "HC37", "HC38", "HC39", "HC40",
                 "HC41", "HC42", "HC43", "HC44", "HC45", "HC46", "HC47", "HC48", "HC49", "HC50",
                 "PM01", "PM02", "PM03", "PM04", "PM05", "PM06", "PM07", "PM08", "PM09", "PM10",
                 "PM11", "PM12", "PM13", "PM14", "PM15", "PM16", "PM17", "PM18", "PM19",
                 "PM10-PRI", "PM25-PRI")

#Emission column names
emis.header  <- c("dayav","HR01","HR02","HR03","HR04","HR05","HR06","HR07","HR08",
                   "HR09","HR10","HR11","HR12","HR13","HR14","HR15","HR16","HR17",
                   "HR18","HR19","HR20","HR21","HR22","HR23","HR24")

output.flag = TRUE

#--------------------------------------------------------
# Read Domain
#--------------------------------------------------------
print("Reading domain file...")

#domain <- read.dbf(domain.file, as.is = FALSE)
domain <- readRDS(domain.file)

#Cycle by sector
for(s in 1:length(sector.list))
{
  
  sector <- sector.list[s]
  
  #Cycle by variable
  for(v in 1:length(var.list))
  {
    pollutant <- var.list[v]
    
    #Cycle by day-of-week
    for(d in 1:length(dow.list))
    {
      dow <- dow.list[d]
      
      input.dir <- paste(base.dir,"/",sector,"/",dow,sep="")
      
      #---------------------------------------
      #Get emissions text file data
      #---------------------------------------
      setwd(input.dir)
      
      file.list <- grep(pattern=paste("/",pollutant,".gz",sep=""),list.files(input.dir,recursive=TRUE),value=TRUE)
      
      exist.count <- 0
      
      #Check if file exists
      if(length(file.list)>0)
      {
        
        print(paste("Processing ",pollutant," for ",sector," on ",dow,"...",sep=""))
        
        #Cycles through hours of the day, note: h=1 corresponds to dayav
        for(h in 1:25)
        {
          
          #Read text file
          input.file <- paste(input.dir,file.list[h],sep="/")
          
          print(paste("Processing ",emis.header[h],"...",sep=""))
          
          data <- read.table(input.file,header=FALSE, sep="")
          
          #Get dimensions of text file matrix
          n.cols <- ncol(data)
          n.rows <- nrow(data)
          
          #Initialize single column array
          temp.emis <- rep(0,n.cols*n.rows)
          
          #Convert text file matrix to single column array
          temp.emis <- as.vector(t(data))
          
          #Append emissions (by hour) to domain file
          temp.df <- data.frame(temp.emis)
          colnames(temp.df) <- emis.header[h]
          
          #Convert from short tons to metric tons
          if(substr(pollutant,1,2) != "HC")
          {
            temp.df <- temp.df * 0.907
          }
            
            emis.df <- cbind(domain,temp.df)
            
          }else{
            
            emis.df <- cbind(emis.df,temp.df)
            
          } #End if statement combining data frames
          
          exist.count <- exist.count + 1
          
        } #End for loop by hour
        
      }else{
        
        print(paste("Missing ",pollutant," for ",sector," on ",dow,"...",sep=""))
        
      } #End if file exists
      
      
      #--------------------------------------------------------
      # Output Data
      #--------------------------------------------------------
      if(output.flag & exist.count > 0)
      {
        #Output directories
        #output.dir  <- paste(base.dir,"/",sector,"/Month",sprintf("%02d",1),"/",dow,sep="")
        #out.file    <- paste(sector,"_",pollutant,".rds",sep="")
        output.dir  <- paste(base.dir,"/",sector,"/",dow,sep="")
        out.file    <- paste(sector,"_",pollutant,".rds",sep="")
        #out.file   <- paste(sector,"_",pollutant,".dbf",sep="")
        #out.file   <- paste(sector,"_",pollutant,".Rdata",sep="")
        #out.file   <- paste(sector,"_",pollutant,".csv",sep="")
        
        setwd(output.dir)
        
        saveRDS(emis.df, file = out.file)
        #save(emis.df, file=out.file)
        #write.dbf(emis.df,out.file)
        #write.csv(emis.df,out.file)
        
      } #End output if statement
      
    } #End for loop by day-of-week
    
  } #End for loop by variable
  
} #End for loop by sector

print(proc.time()-ptm)		