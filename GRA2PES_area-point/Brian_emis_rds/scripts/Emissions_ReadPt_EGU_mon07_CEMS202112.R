#--------------------------------------------------------
# bash
# eval "$(/home/bmcdonald/anaconda3/bin/conda shell.bash hook)"
# conda activate r_env
#--------------------------------------------------------

rm(list=ls()) #Clear memory
ptm <- proc.time()
library("ncdf4")
library("foreign")
library(rgdal)
library(sp)
#library(maptools)

#--------------------------------------------------------
# Set input variables
#--------------------------------------------------------
#category    <- "point" 
#base.dir    <- "/wrk/d2/bmcdonald/NEI17"
base.dir    <- "/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_GRA2PES2021/Stu_RELPT/output_RELPT_emis_CEMS202112"

#domain.file <- "/wrk/csd4/bmcdonald/NEI17/Relpointinfo_Criteria_n_VCP.csv"
domain.file <- "/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_GRA2PES2021/Stu_RELPT/output_RELPT_template_CEMS202112/Relpointinfo_Criteria_n_VCP.csv"

#Day-of-week subdirectories
dow.list    <- c("weekdy","satdy","sundy")

month       <- "Month07" ###CHANGE

#Sectors to process
egu.list     <- c("PtEGU_Coal","PtEGU_NG","PtEGU_Oil","PtEGU_BIO")
comm.list    <- c("PtCOMM_Coal","PtCOMM_NG","PtCOMM_Oil","PtCOMM_BIO")
indf.list    <- c("PtIND_Coal","PtIND_NG","PtIND_NG2","PtIND_Oil","PtIND_Oil2","PtIND_BIO")
indp.list    <- c("PtCHEM","PtFOOD","PtMETAL","PtREFINE","PtPULP","PtELECT","PtMOTOR","PtAPPAREL","PtPHOTO","PtDRUG","PtMISC","PtMISC2")
vcp.list     <- c("PtAdhesives","PtCoatings","PtDegreasing","PtInks")
fug.list     <- c("PtOnG","PtSTORAGE","PtEVAPGAS","PtCONST")
mob.list     <- c("PtAVIATION","PtRAIL")

sector.list <- c(egu.list)

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

#Emission column names
emis.header  <- c("dayav","HR01","HR02","HR03","HR04","HR05","HR06","HR07","HR08",
                  "HR09","HR10","HR11","HR12","HR13","HR14","HR15","HR16","HR17",
                  "HR18","HR19","HR20","HR21","HR22","HR23","HR24")

#Set ArcGIS shapefile with state boundaries
state.dir  <- "/wrk/csd4/bmcdonald/NYC_Modeling/Domains"
state.file <- "statesp020"

output.flag  <- TRUE

#---------------------------------------
# Get point file location information
#---------------------------------------
library(readr)
domain <- read_csv(domain.file)#,header=TRUE)

#Keep only lat, lon coordinates
domain <- domain[,c("Lon_deg","Lat_deg")]
colnames(domain) <- c("LON","LAT")

#Create point shapefile of stations
crswgs84 <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
site.sp  <- SpatialPointsDataFrame(coords = domain, data = domain,
                                   proj4string = crswgs84)

#Read state polygon shapefile
states.sp <- readOGR(state.dir,layer=state.file,verbose=TRUE)
proj4string(states.sp) <- crswgs84

#Intersect shapefiles (add state name to point file)
state.data = over(site.sp,states.sp)
state.data = data.frame(state.data$STATE,state.data$STATE_FIPS)
colnames(state.data) <- c("STATE","STATE_FIPS")

#Merge site location data with state information
domain <- cbind(site.sp@data,state.data)

#---------------------------------------
# Process Emission Text Files
#---------------------------------------
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

      #input.dir <- paste(base.dir,"/",category,"/",month,"/",sector,"/",dow,sep="")
      input.dir <- paste(base.dir,"/",month,"/",sector,"/",dow,sep="")      

      #---------------------------------------
      #Get emissions text file data
      #---------------------------------------
      setwd(input.dir)

      print(paste("Processing ",pollutant," for ",sector," on ",dow,"...",sep=""))

      file.list <- grep(pattern=paste("/",pollutant,".gz",sep=""),list.files(input.dir,recursive=TRUE),value=TRUE)
      
      exist.count <- 0
      
      if(length(file.list)>0)
      {
        
        #Cycles through hours of the day, note: h=1 corresponds to dayav
        for(h in 1:25)
        {
          
          print(paste("Processing ",emis.header[h],"...",sep=""))
          
          #Read text file
          input.file <- paste(input.dir,file.list[h],sep="/")
          
          data <- read.table(input.file,fill=TRUE,header=FALSE, sep="")
          
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
          
          #Remove NA values
          temp.df <- na.omit(temp.df)
          
          #Convert from short tons to metric tons
          if(substr(pollutant,1,2) != "HC")
          {
            temp.df <- temp.df * 0.907
          }
          
          if(h==1)
          {
            
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
        #output.dir  <- paste(base.dir,"/",category,"/",month,"/",sector,"/",dow,sep="")
        output.dir  <- paste(base.dir,"/",month,"/",sector,"/",dow,sep="")        
        setwd(output.dir)

        out.file <- paste(sector,"_",pollutant,"_",dow,".rds",sep="")
        saveRDS(emis.df, file = out.file)

        #out.file <- paste(sector,"_",pollutant,"_",dow,".dbf",sep="")
        #write.dbf(emis.df,out.file)

      } #End output if statement

    } #End for loop by day-of-week

  } #End for loop by variable

} #End for loop by sector

print(proc.time()-ptm)	
