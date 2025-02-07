rm(list=ls()) #Clear memory
ptm <- proc.time()
library(foreign)
#library(rgdal)
#library(sp)
#library(maptools)

#--------------------------------------------------------
# Set input variables
#--------------------------------------------------------
#Set ArcGIS shapefile with domain
domain.dir   <- "/wrk/csd4/bmcdonald/COVID_Modeling/Domains"
domain.name  <- "nei04k_domain"

#Set ArcGIS shapefile with on-road CO2 mapped
base.dir     <- "/wrk/users/charkins/emissions/FIVE/offroad_marGASnomonthly_marGASnoDOW_ffCO2/2021/Month04"
input.dir    <- "/wrk/charkins/FIVE/2019_2020_inputdir"
five.dir     <- paste(input.dir,"",sep="")
five.dir.wd  <- paste(base.dir,"/weekdy",sep="")
five.dir.sa  <- paste(base.dir,"/satdy",sep="")
five.dir.su  <- paste(base.dir,"/sundy",sep="")

nei.dir      <- "/wrk/charkins/FIVE/nei_nonroad_shp"
agDSL.shp    <- "Ag_Diesel_emis_tpd"
offDSL.shp   <- "DieselOffroad_emis_tpd"
offGAS2.shp  <- "Gas2StrkOffrd_emis_tpd"
offGAS4.shp  <- "Gas4StrkOffrd_emis_tpd"
marGAS2.shp  <- "Gas2StrkMarine_emis_tpd"
marGAS4.shp  <- "Gas4StrkMarine_emis_tpd"

five.name    <- "offroad2019.rds"

#Onroad fuel scaling factors
fuel.file    <- "offroad_fuel2019.csv"

#Fuel-based emission factors
pollutant.list <- c("CO2","ffCO2", "CO", "NOX", "NH3", "SO2", "VOC_EXH", "PM25_PRI_EXH","PM10_PRI_EXH")

pollutant.marGAS.ef.file   <- "offroad_marGASef2021_ffCO2.csv"    #run exhuast only
pollutant.twoGAS.ef.file   <- "offroad_twoGASef2021_ffCO2.csv"    #run exhuast only
pollutant.fourGAS.ef.file  <- "offroad_fourGASef2021_ffCO2.csv"   #run exhuast only
pollutant.nonagDSL.ef.file <- "offroad_nonagDSLef2021_ffCO2.csv"  #run exhuast only
pollutant.agDSL.ef.file    <- "offroad_agDSLef2021_ffCO2.csv"     #run exhuast only

#Notes:
#1. Exclude ptVCP use from industrial coatings, inks, and adhesives
#2. Account for 20% reduction in VOCs in coatings in AIM areas (based on ARB report)
#3. Printing ink shops in LA and NYC assumed to have maximum controls
#4. Account for CA, NY, LA, and NYC regulations on coatings
#5. Account for CA and NY regulations on adhesives & sealants

#Speciation files (number of speciation profiles must match number of sectors)
voc.wt.file    <- "offroad_wt.csv"
voc.mw.file    <- "offroad_mw.csv"
pm25.wt.file   <- "offroad_pm25.csv"

#Diurnal files (number of diurnal profiles must match number of sectors)
month               <- 4
monthly.file        <- "offroad_monthly_2021_noMarGasMonthlyScale.csv"
diurnal.file        <- "offroad_diurnal_noMarGASDOW.csv" #local time
dow.list            <- c("weekdy","satdy","sundy")
timezone.list  <- c("Atlantic","Eastern","Central","Mountain","Pacific","Alaska")
utc.list            <- c(-3,-4,-5,-6,-7,-8)

#Set output directory
output.dir  <- base.dir

#Set flags to turn on functions
read.domain.flag   <- FALSE
read.fuel.flag     <- FALSE
speciate.voc.flag  <- TRUE
speciate.pm25.flag <- TRUE
output.flag        <- TRUE

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

  print("Reading domain file...")
  
  #Read domain polygon shapefile
  domain.sp <- readRDS(paste(domain.name,".rds",sep=""))
  
}

#--------------------------------------------------------
# Estimate Offroad Fuel Data
#--------------------------------------------------------
if(read.fuel.flag)
{
  setwd(nei.dir)

  #------------------------------------
  # Read off-road NOx emissions (mt/d)
  #------------------------------------
  print("Reading agricultural diesel NOx emissions file...")
  agDSL.sp   <- read.dbf(paste(agDSL.shp,   ".dbf", sep="")) #Agricultural diesel 
  
  print("Reading total offroad diesel NOx emissions file...")
  offDSL.sp  <- read.dbf(paste(offDSL.shp,  ".dbf", sep="")) #Total offroad diesel
  
  print("Reading 2-stroke offroad gasoline NOx emissions file...")
  offGAS2.sp <- read.dbf(paste(offGAS2.shp, ".dbf", sep="")) #Two-stroke offroad gasoline
  
  print("Reading 4-stroke offroad gasoline NOx emissions file...")  
  offGAS4.sp <- read.dbf(paste(offGAS4.shp, ".dbf", sep="")) #Four-stroke offroad gasoline
  
  print("Reading 2-stroke marine gasoline NOx emissions file...")
  marGAS2.sp <- read.dbf(paste(marGAS2.shp, ".dbf", sep="")) #Two-stroke marine gasoline
  
  print("Reading 4-stroke marine gasoline NOx emissions file...")
  marGAS4.sp <- read.dbf(paste(marGAS4.shp, ".dbf", sep="")) #Four-stroke marine gasoline

  #Combine into one file
  nei.sp <- cbind(domain.sp,
                   marGAS2.sp$avgNOX + marGAS4.sp$avgNOX,
                   offGAS2.sp$avgNOX,
                   offGAS4.sp$avgNOX,
                   offDSL.sp$avgNOX - agDSL.sp$avgNOX,
                   agDSL.sp$avgNOX)
  
  colnames(nei.sp) <- c(colnames(domain.sp), "marGAS", "twoGAS", "fourGAS", "nonagDSL", "agDSL")
  
  #------------------------------------
  # Calculate Current Year Fuel Data
  #------------------------------------
  setwd(five.dir)

  #Get fuel sales by state (in kg/d)
  fuel.df <- read.csv(fuel.file, header = TRUE)

  #Initialize arrays
  domain.sp$Id <- paste(domain.sp$Row, domain.sp$Col, sep=",")
  nei.sp$Id    <- paste(nei.sp$Row, nei.sp$Col, sep=",")
 
  marGAS.index   <- ncol(domain.sp)+1
  twoGAS.index   <- ncol(domain.sp)+2
  fourGAS.index  <- ncol(domain.sp)+3
  nonagDSL.index <- ncol(domain.sp)+4
  agDSL.index    <- ncol(domain.sp)+5

  domain.sp$marGAS   <- 0
  domain.sp$twoGAS   <- 0
  domain.sp$fourGAS  <- 0
  domain.sp$nonagDSL <- 0
  domain.sp$agDSL    <- 0

  max.row <- max(domain.sp$Row) + 1
  max.col <- max(domain.sp$Col) + 1
  
  #Normalize NOx emissions by state and multiple with fuel sales(except marine)
  for(s in 1:nrow(fuel.df))
  {
    
    state.name <- fuel.df$Name[s]
    
    print(paste("Processing fuel sales for ",state.name,"...",sep=""))
    
    twoGAS.fuel   <- fuel.df[s,grep("^twoGAS$",colnames(fuel.df),value=TRUE)]   #in kg/d
    fourGAS.fuel  <- fuel.df[s,grep("^fourGAS$",colnames(fuel.df),value=TRUE)]  #in kg/d
    nonagDSL.fuel <- fuel.df[s,grep("^nonagDSL$",colnames(fuel.df),value=TRUE)] #in kg/d
    agDSL.fuel    <- fuel.df[s,grep("^agDSL$",colnames(fuel.df),value=TRUE)]    #in kg/d
    
    if(state.name != "US")
    {

      state.pos <- grep(paste("^",state.name,"$",sep=""), nei.sp$STATE_NAME)
            
      #Calculate state total NOx emissions (except marine)
      if(length(state.pos) > 0)
      {
        
        twoGAS.nox   <- sum(nei.sp[state.pos,c("twoGAS")])
        fourGAS.nox  <- sum(nei.sp[state.pos,c("fourGAS")])
        nonagDSL.nox <- sum(nei.sp[state.pos,c("nonagDSL")])
        agDSL.nox    <- sum(nei.sp[state.pos,c("agDSL")])
        
        #Calculate fuel sales in each grid cell by state (except marine)
        if(twoGAS.nox > 0){
          
          domain.sp[state.pos,c("twoGAS")]   <- twoGAS.fuel   * nei.sp[state.pos,c("twoGAS")]   / twoGAS.nox
          
        } #End twoGAS if statement
        
        if(fourGAS.nox > 0){
          
          domain.sp[state.pos,c("fourGAS")]  <- fourGAS.fuel  * nei.sp[state.pos,c("fourGAS")]  / fourGAS.nox       
          
        } #End fourGAS if statement
        
        if(nonagDSL.nox > 0){
          
          domain.sp[state.pos,c("nonagDSL")] <- nonagDSL.fuel * nei.sp[state.pos,c("nonagDSL")] / nonagDSL.nox        
          
        } #End nonagDsl if statement
        
        if(agDSL.nox > 0){
          
          domain.sp[state.pos,c("agDSL")]    <- agDSL.fuel    * nei.sp[state.pos,c("agDSL")]    / agDSL.nox        
          
        } #End agDsl if statement
        
      } #End statepos if statement

      
    }else{
      
      #Calculate marine fuel sales (normalize nationally)
      marGAS.fuel  <- fuel.df[grep("^US$",fuel.df$Name),grep(paste("^marGAS$",sep=""),colnames(fuel.df),value=TRUE)] #in kg/d
      marGAS.nox   <- sum(nei.sp[,marGAS.index])
      
      if(marGAS.nox > 0)
      {
        
        domain.sp[,marGAS.index] <- marGAS.fuel * nei.sp[,c("marGAS")] / marGAS.nox
        
      } #End marGAS if statement
      
    } #End state if statement
    
  } #End state for loop

  #------------------------------------
  # Save Updated FIVE Fuel Data
  #------------------------------------
  setwd(five.dir)

  saveRDS(domain.sp, file = five.name)

}else{

  print("Reading offroad fuel sales (in kg/d) emissions file...")

  #Read population polygon shapefile
  setwd(five.dir)
  domain.sp <- readRDS(five.name)

}

#--------------------------------------------------------
# Calculate On-Road Co-Emitted Pollutant Emissions
#--------------------------------------------------------
for(p in 1:length(pollutant.list))
{
  setwd(five.dir)

  pollutant.name <- pollutant.list[p]
  print(paste("Processing ",pollutant.name," emissions...",sep=""))

  #------------------------------------
  # Estimate Daily Emissions
  #------------------------------------
  marGAS.index   <- grep("^marGAS$",colnames(domain.sp))
  twoGAS.index   <- grep("^twoGAS$",colnames(domain.sp))
  fourGAS.index  <- grep("^fourGAS$",colnames(domain.sp))
  nonagDSL.index <- grep("^nonagDSL$",colnames(domain.sp))
  agDSL.index    <- grep("^agDSL$",colnames(domain.sp))  

  #Initialize sectoral emission dataframes
  marGAS.emis   <- data.frame(rep(0,nrow(domain.sp)))
  twoGAS.emis   <- data.frame(rep(0,nrow(domain.sp)))
  fourGAS.emis  <- data.frame(rep(0,nrow(domain.sp)))
  nonagDSL.emis <- data.frame(rep(0,nrow(domain.sp)))
  agDSL.emis    <- data.frame(rep(0,nrow(domain.sp)))

  colnames(marGAS.emis)   <- paste(colnames(domain.sp)[marGAS.index],"_",pollutant.name,sep="")
  colnames(twoGAS.emis)   <- paste(colnames(domain.sp)[twoGAS.index],"_",pollutant.name,sep="")
  colnames(fourGAS.emis)  <- paste(colnames(domain.sp)[fourGAS.index],"_",pollutant.name,sep="")
  colnames(nonagDSL.emis) <- paste(colnames(domain.sp)[nonagDSL.index],"_",pollutant.name,sep="")
  colnames(agDSL.emis)    <- paste(colnames(domain.sp)[agDSL.index],"_",pollutant.name,sep="")

  #Read emission factors
  marGAS.ef.df   <- read.csv(pollutant.marGAS.ef.file, header = TRUE)
  twoGAS.ef.df   <- read.csv(pollutant.twoGAS.ef.file, header = TRUE)
  fourGAS.ef.df  <- read.csv(pollutant.fourGAS.ef.file, header = TRUE)
  nonagDSL.ef.df <- read.csv(pollutant.nonagDSL.ef.file, header = TRUE)
  agDSL.ef.df    <- read.csv(pollutant.agDSL.ef.file, header = TRUE)

  # Read monthly PADD Scalings
  monthly.df <- read.csv(monthly.file, header = TRUE)
  
  #Multiply with on-road fuel EF by state
  for(s in 1:nrow(marGAS.ef.df))
  {
    state.name <- marGAS.ef.df$Name[s]
    
    marGAS.factor   <- marGAS.ef.df[s,grep(paste("^",pollutant.name,"$",sep=""),colnames(marGAS.ef.df),value=TRUE)]
    twoGAS.factor   <- twoGAS.ef.df[s,grep(paste("^",pollutant.name,"$",sep=""),colnames(twoGAS.ef.df),value=TRUE)]
    fourGAS.factor  <- fourGAS.ef.df[s,grep(paste("^",pollutant.name,"$",sep=""),colnames(fourGAS.ef.df),value=TRUE)]
    nonagDSL.factor <- nonagDSL.ef.df[s,grep(paste("^",pollutant.name,"$",sep=""),colnames(nonagDSL.ef.df),value=TRUE)]
    agDSL.factor    <- agDSL.ef.df[s,grep(paste("^",pollutant.name,"$",sep=""),colnames(agDSL.ef.df),value=TRUE)]

    if(state.name != "US")
    {
      state.pos <- grep(paste("^",state.name,"$",sep=""), domain.sp$STATE_NAME)
      padd.no <- round(mean(domain.sp[state.pos,c("PADD")]),0)
    }else{
      state.pos <- which(is.na(domain.sp$STATE_NAME))
      padd.no <- "US"
    }
    
    # Subset monthly padd scalings and monthly padd to state scalings
    monthly.sub <- subset(monthly.df, PADD == padd.no)
    
    #Monthly factors
    marGAS.month.factor   <- monthly.sub[month, grep("^marGAS$",  colnames(monthly.sub))]  #select by sector (month)
    twoGAS.month.factor   <- monthly.sub[month, grep("^twoGAS$",  colnames(monthly.sub))]  #select by sector (month)
    fourGAS.month.factor  <- monthly.sub[month, grep("^fourGAS$", colnames(monthly.sub))]  #select by sector (month)
    nonagDSL.month.factor <- monthly.sub[month, grep("^nonagDSL$", colnames(monthly.sub))] #select by sector (month)
    agDSL.month.factor    <- monthly.sub[month, grep("^agDSL$", colnames(monthly.sub))]    #select by sector (month)

    marGAS.emis[state.pos,]   <- domain.sp[state.pos,marGAS.index]   * marGAS.factor   / 10^6  * marGAS.month.factor   #Convert g/d to mt/d
    twoGAS.emis[state.pos,]   <- domain.sp[state.pos,twoGAS.index]   * twoGAS.factor   / 10^6  * twoGAS.month.factor   #Convert g/d to mt/d
    fourGAS.emis[state.pos,]  <- domain.sp[state.pos,fourGAS.index]  * fourGAS.factor  / 10^6  * fourGAS.month.factor  #Convert g/d to mt/d
    nonagDSL.emis[state.pos,] <- domain.sp[state.pos,nonagDSL.index] * nonagDSL.factor / 10^6  * nonagDSL.month.factor #Convert g/d to mt/d
    agDSL.emis[state.pos,]    <- domain.sp[state.pos,agDSL.index]    * agDSL.factor    / 10^6  * agDSL.month.factor    #Convert g/d to mt/d

  }

  #------------------------------------
  # Estimate Diurnal Emissions
  #------------------------------------
  
  diurnal.df <- read.csv(diurnal.file, header = TRUE)
  
  #Cycle by day of week
  for(d in 1:length(dow.list))
  {
    dow.name <- dow.list[d]

    #Initalize emissions data frame
    emis.df    <- data.frame(matrix(0, nrow = nrow(domain.sp), ncol = 25))
    colnames(emis.df) <- c("dayav","HR00","HR01","HR02","HR03","HR04","HR05",
                           "HR06","HR07","HR08","HR09","HR10","HR11","HR12",
                           "HR13","HR14","HR15","HR16","HR17","HR18","HR19",
                           "HR20","HR21","HR22","HR23")

    emis.marGAS.df    <- emis.df
    emis.twoGAS.df    <- emis.df
    emis.fourGAS.df   <- emis.df
    emis.nonagDSL.df  <- emis.df
    emis.agDSL.df     <- emis.df

    #Select diurnal profiles for day-of-week
    if(dow.name == "weekdy")
    {

      weekdy.df <- diurnal.df[grep("weekdy",diurnal.df$Day),]

      diurnal.dow.df <- rbind(weekdy.df, weekdy.df)

    }else if(dow.name == "satdy"){

      weekdy.df <- diurnal.df[grep("weekdy",diurnal.df$Day),]
      satdy.df <- diurnal.df[grep("satdy",diurnal.df$Day),]

      diurnal.dow.df <- rbind(weekdy.df, satdy.df)

    }else if(dow.name == "sundy"){

      satdy.df <- diurnal.df[grep("satdy",diurnal.df$Day),]
      sundy.df <- diurnal.df[grep("sundy",diurnal.df$Day),]

      diurnal.dow.df <- rbind(satdy.df, sundy.df)

    }

    #Diurnal factors
    marGAS.col    <- grep("^marGAS$",   colnames(diurnal.dow.df)) #select by sector (dow)
    twoGAS.col    <- grep("^twoGAS$",   colnames(diurnal.dow.df)) #select by sector (dow)
    fourGAS.col   <- grep("^fourGAS$",  colnames(diurnal.dow.df)) #select by sector (dow)
    nonagDSL.col  <- grep("^nonagDSL$", colnames(diurnal.dow.df)) #select by sector (dow)
    agDSL.col     <- grep("^agDSL$",    colnames(diurnal.dow.df)) #select by sector (dow)

    #Adjust for dow and month
    marGAS.emis.dow   <- marGAS.emis*0
    twoGAS.emis.dow   <- twoGAS.emis*0
    fourGAS.emis.dow  <- fourGAS.emis*0
    nonagDSL.emis.dow <- nonagDSL.emis*0
    agDSL.emis.dow    <- agDSL.emis*0

    emis.df[,c("dayav")] <- marGAS.emis.dow + twoGAS.emis.dow + fourGAS.emis.dow + nonagDSL.emis.dow + agDSL.emis.dow
    
    emis.marGAS.df[,c("dayav")]   <- marGAS.emis.dow
    emis.twoGAS.df[,c("dayav")]   <- twoGAS.emis.dow
    emis.fourGAS.df[,c("dayav")]  <- fourGAS.emis.dow
    emis.nonagDSL.df[,c("dayav")] <- nonagDSL.emis.dow
    emis.agDSL.df[,c("dayav")]    <- agDSL.emis.dow

    #Cycle by timezone
    for(t in 1:length(timezone.list))
    {

      tz.name   <- timezone.list[t]
      tz.offset <- -utc.list[t]

      print(paste("Performing temporal allocation for ",dow.name," in ",tz.name," timezone...",sep=""))

      #Create diurnal factor table that accounts for different timezones
      start.pos  <- (24 - tz.offset)
      end.pos    <- start.pos + 24 - 1
      #start.pos <- (24 - tz.offset) + 1
      #end.pos   <- start.pos + 24 - 1

      diurnal.utc.df <- diurnal.dow.df[c(start.pos:end.pos),]
      diurnal.utc.df[1:24,c("Day")]  <- rep(dow.name,24)
      diurnal.utc.df[1:24,c("Hour")] <- c(0:23)
    
      #diurnal.utc.df <- diurnal.dow.df[start.pos:end.pos, 2:ncol(diurnal.dow.df)]

      # for(h in 1:24)
      # {
      # 
      #  #Add hours to convert from local time to UTC time
      #  if(h+tz.offset <= 24)
      #  {
      # 
      #    diurnal.utc.df[h+tz.offset, 2:ncol(diurnal.df)] <- diurnal.dow.df[start.pos, 2:ncol(diurnal.dow.df)]
      # 
      #  }else{
      # 
      #    diurnal.utc.df[h+tz.offset-24, 2:ncol(diurnal.df)] <- diurnal.dow.df[h, 2:ncol(diurnal.df)]
      #  }
      # 
      # } #End hour for loop

      #Multiply daily average emissions to hourly emissions
      tz.indices       <- grep(tz.name, domain.sp$TZ)            #select by timezone

      #Add emissions
      if(length(tz.indices) > 0)
      {

        daily.marGAS     <- matrix(rep(marGAS.emis[tz.indices,],each=24), ncol=24, byrow=TRUE)
        diurnal.marGAS   <- matrix(rep(diurnal.utc.df[,marGAS.col],each=nrow(daily.marGAS)), nrow=nrow(daily.marGAS))
        
        daily.twoGAS     <- matrix(rep(twoGAS.emis[tz.indices,],each=24), ncol=24, byrow=TRUE)
        diurnal.twoGAS   <- matrix(rep(diurnal.utc.df[,twoGAS.col],each=nrow(daily.twoGAS)), nrow=nrow(daily.twoGAS))
        
        daily.fourGAS    <- matrix(rep(fourGAS.emis[tz.indices,],each=24), ncol=24, byrow=TRUE)
        diurnal.fourGAS  <- matrix(rep(diurnal.utc.df[,fourGAS.col],each=nrow(daily.fourGAS)), nrow=nrow(daily.fourGAS))
        
        daily.nonagDSL   <- matrix(rep(nonagDSL.emis[tz.indices,],each=24), ncol=24, byrow=TRUE)
        diurnal.nonagDSL <- matrix(rep(diurnal.utc.df[,nonagDSL.col],each=nrow(daily.nonagDSL)), nrow=nrow(daily.nonagDSL))
        
        daily.agDSL      <- matrix(rep(agDSL.emis[tz.indices,],each=24), ncol=24, byrow=TRUE)
        diurnal.agDSL    <- matrix(rep(diurnal.utc.df[,agDSL.col],each=nrow(daily.agDSL)), nrow=nrow(daily.agDSL))
        
        emis.df[tz.indices,2:25] <- emis.df[tz.indices,2:25] + (daily.marGAS   / 24 * diurnal.marGAS)
        emis.df[tz.indices,2:25] <- emis.df[tz.indices,2:25] + (daily.twoGAS   / 24 * diurnal.twoGAS)
        emis.df[tz.indices,2:25] <- emis.df[tz.indices,2:25] + (daily.fourGAS  / 24 * diurnal.fourGAS)
        emis.df[tz.indices,2:25] <- emis.df[tz.indices,2:25] + (daily.nonagDSL / 24 * diurnal.nonagDSL)
        emis.df[tz.indices,2:25] <- emis.df[tz.indices,2:25] + (daily.agDSL    / 24 * diurnal.agDSL)
        
        emis.marGAS.df[tz.indices,2:25]   <- emis.marGAS.df[tz.indices,2:25]   + (daily.marGAS   / 24 * diurnal.marGAS)
        emis.twoGAS.df[tz.indices,2:25]   <- emis.twoGAS.df[tz.indices,2:25]   + (daily.twoGAS   / 24 * diurnal.twoGAS)
        emis.fourGAS.df[tz.indices,2:25]  <- emis.fourGAS.df[tz.indices,2:25]  + (daily.fourGAS  / 24 * diurnal.fourGAS)
        emis.nonagDSL.df[tz.indices,2:25] <- emis.nonagDSL.df[tz.indices,2:25] + (daily.nonagDSL / 24 * diurnal.nonagDSL)
        emis.agDSL.df[tz.indices,2:25]    <- emis.agDSL.df[tz.indices,2:25]    + (daily.agDSL    / 24 * diurnal.agDSL)

      }

    } #End time zone for loop
    
    emis.df[,1]           <- rowSums(emis.df[,2:25])
    emis.marGAS.df[,1]    <- rowSums(emis.marGAS.df[,2:25])
    emis.twoGAS.df[,1]    <- rowSums(emis.twoGAS.df[,2:25])
    emis.fourGAS.df[,1]   <- rowSums(emis.fourGAS.df[,2:25])
    emis.nonagDSL.df[,1]  <- rowSums(emis.nonagDSL.df[,2:25])
    emis.agDSL.df[,1]     <- rowSums(emis.agDSL.df[,2:25])

    #------------------------------------
    # Output Pollutant Emissions
    #------------------------------------
    if(output.flag)
    {
      setwd(paste(output.dir,"/",dow.name,sep=""))

      print(paste("Writing output of ",pollutant.name," for ",dow.name,"...",sep=""))

      #Reorder HR00 to HR24
      emis24.df <- emis.df[,c("dayav")]
      emis24.df <- cbind(emis24.df,emis.df[,c(3:ncol(emis.df))])
      emis24.df <- cbind(emis24.df,emis.df[,c("HR00")])
      colnames(emis24.df)[1] <- c("dayav")
      colnames(emis24.df)[ncol(emis24.df)] <- c("HR24")

      out.file <- paste("Offroad","_",pollutant.name,".rds",sep="")
      out.df <- cbind(domain.sp[,c("Id","Row","Col","LON","LAT","TZ","URB_RUR")], emis24.df)

      #write.csv(out.df, out.file, row.names = FALSE)
      saveRDS(out.df, file=out.file)
    }
 
    #------------------------------------
    # VOC Speciation
    #------------------------------------
    if(speciate.voc.flag)
    {
      setwd(five.dir)

      voc.wt.df <- read.csv(voc.wt.file, header = TRUE)
      voc.mw.df <- read.csv(voc.mw.file, header = TRUE)

      #Cycle through each hydrocarbon bin
      for(bin in 1:nrow(voc.wt.df))
      {

        bin.name <- voc.wt.df$Bin[bin]

        #Obtain VOC weight fraction (in percent) and MW
        if(pollutant.name == "VOC_EXH")
        {

          print(paste("Processing VOC EXH for ",bin.name,"...",sep=""))

          marGAS.exh.col    <- grep("^marGAS$",   colnames(voc.wt.df))
          twoGAS.exh.col    <- grep("^twoGAS$",   colnames(voc.wt.df))
          fourGAS.exh.col   <- grep("^fourGAS$",  colnames(voc.wt.df))
          nonagDSL.exh.col  <- grep("^nonagDSL$", colnames(voc.wt.df))
          agDSL.exh.col     <- grep("^agDSL$",    colnames(voc.wt.df))

          marGAS.exh.voc.wt   <- voc.wt.df[bin,"marGAS"]
          marGAS.exh.voc.mw   <- voc.mw.df[bin,"marGAS"]
          twoGAS.exh.voc.wt   <- voc.wt.df[bin,"twoGAS"]
          twoGAS.exh.voc.mw   <- voc.mw.df[bin,"twoGAS"]
          fourGAS.exh.voc.wt  <- voc.wt.df[bin,"fourGAS"]
          fourGAS.exh.voc.mw  <- voc.mw.df[bin,"fourGAS"]
          nonagDSL.exh.voc.wt <- voc.wt.df[bin,"nonagDSL"]
          nonagDSL.exh.voc.mw <- voc.mw.df[bin,"nonagDSL"]
          agDSL.exh.voc.wt    <- voc.wt.df[bin,"agDSL"]
          agDSL.exh.voc.mw    <- voc.mw.df[bin,"agDSL"]

          #Estimate Marine Gasoline VOC emissions (in moles)
          if(is.na(marGAS.exh.voc.mw))
          {
            #Create blank array if MW does not exist
            marGAS.exh.voc.df <- data.frame(matrix(0, nrow=nrow(emis.marGAS.df), ncol=ncol(emis.marGAS.df)))
            colnames(marGAS.exh.voc.df) <- colnames(emis.marGAS.df)

          }else{

            #Save molar emissions if MW exists
            marGAS.exh.voc.df <- emis.marGAS.df * (marGAS.exh.voc.wt / 100) * 10^6 / marGAS.exh.voc.mw

          } #End marine gasoline exhaust if statement

          #Estimate Two-Stroke Gasoline VOC emissions (in moles)
          if(is.na(twoGAS.exh.voc.mw))
          {
            #Create blank array if MW does not exist
            twoGAS.exh.voc.df <- data.frame(matrix(0, nrow=nrow(emis.twoGAS.df), ncol=ncol(emis.twoGAS.df)))
            colnames(twoGAS.exh.voc.df) <- colnames(emis.twoGAS.df)
            
          }else{
            
            #Save molar emissions if MW exists
            twoGAS.exh.voc.df <- emis.twoGAS.df * (twoGAS.exh.voc.wt / 100) * 10^6 / twoGAS.exh.voc.mw
            
          } #End two-stroke gasoline exhaust if statement
          
          #Estimate Four-Stroke Gasoline VOC emissions (in moles)
          if(is.na(fourGAS.exh.voc.mw))
          {
            #Create blank array if MW does not exist
            fourGAS.exh.voc.df <- data.frame(matrix(0, nrow=nrow(emis.fourGAS.df), ncol=ncol(emis.fourGAS.df)))
            colnames(fourGAS.exh.voc.df) <- colnames(emis.fourGAS.df)
            
          }else{
            
            #Save molar emissions if MW exists
            fourGAS.exh.voc.df <- emis.fourGAS.df * (fourGAS.exh.voc.wt / 100) * 10^6 / fourGAS.exh.voc.mw
            
          } #End four-stroke gasoline exhaust if statement
          
          #Estimate non-Agriculture Diesel VOC emissions (in moles)
          if(is.na(nonagDSL.exh.voc.mw))
          {
            #Create blank array if MW does not exist
            nonagDSL.exh.voc.df <- data.frame(matrix(0, nrow=nrow(emis.nonagDSL.df), ncol=ncol(emis.nonagDSL.df)))
            colnames(nonagDSL.exh.voc.df) <- colnames(emis.nonagDSL.df)
            
          }else{
            
            #Save molar emissions if MW exists
            nonagDSL.exh.voc.df <- emis.nonagDSL.df * (nonagDSL.exh.voc.wt / 100) * 10^6 / nonagDSL.exh.voc.mw
            
          } #End non-agricultural diesel exhaust if statement
          
          #Estimate non-Agriculture Diesel VOC emissions (in moles)
          if(is.na(agDSL.exh.voc.mw))
          {
            #Create blank array if MW does not exist
            agDSL.exh.voc.df <- data.frame(matrix(0, nrow=nrow(emis.agDSL.df), ncol=ncol(emis.agDSL.df)))
            colnames(agDSL.exh.voc.df) <- colnames(emis.agDSL.df)
            
          }else{
            
            #Save molar emissions if MW exists
            agDSL.exh.voc.df <- emis.agDSL.df * (agDSL.exh.voc.wt / 100) * 10^6 / agDSL.exh.voc.mw
            
          } #End non-agricultural diesel exhaust if statement

          #------------------------------------
          # Output VOC EXH Emissions (in moles)
          #------------------------------------
          tot.exh.voc.df <- marGAS.exh.voc.df + twoGAS.exh.voc.df + fourGAS.exh.voc.df + nonagDSL.exh.voc.df + agDSL.exh.voc.df

          if(output.flag & pollutant.name == "VOC_EXH")
          {
            setwd(paste(output.dir,"/",dow.name,sep=""))

            #Reorder HR00 to HR24
            voc24.df <- tot.exh.voc.df[,c("dayav")]
            voc24.df <- cbind(voc24.df, tot.exh.voc.df[,c(3:ncol(tot.exh.voc.df))])
            voc24.df <- cbind(voc24.df, tot.exh.voc.df[,c("HR00")])
            colnames(voc24.df)[1] <- c("dayav")
            colnames(voc24.df)[ncol(voc24.df)] <- c("HR24")

            #voc.out <- cbind(domain.sp[,c("Id","Row","Col","LON","LAT","TZ","URB_RUR")], tot.exh.voc.df)
            voc.out  <- cbind(domain.sp[,c("Id","Row","Col","LON","LAT","TZ","URB_RUR")], voc24.df)
            voc.file <- paste("Offroad","_",bin.name,"_EXH.rds",sep="")

            #write.csv(voc.out, file = voc.file, row.names = FALSE)
            saveRDS(voc.out, file=voc.file)
            setwd(five.dir)

          } #End hydrocarbon output if statement

        } #End exh if statement

      } #End hydrocarbon bin for loop

    } #End VOC speciation if statement

    #------------------------------------
    # PM2.5 Speciation
    #------------------------------------
    if(speciate.pm25.flag)
    {
      setwd(five.dir)

      pm25.wt.df <- read.csv(pm25.wt.file, header = TRUE)

      #Cycle through each PM2.5 bin
      for(pm.bin in 1:nrow(pm25.wt.df))
      {
        pm.bin.name <- pm25.wt.df$Bin[pm.bin]

        #------------------------------------
        # Calculate Exhaust PM (in t/d)
        #------------------------------------
        if(pollutant.name == "PM25_PRI_EXH")
        {
          print(paste("Processing PM25 EXH for ",pm.bin.name,"...",sep=""))

          marGAS.exh.pm25.wt   <- pm25.wt.df[pm.bin,"marGAS"]
          twoGAS.exh.pm25.wt   <- pm25.wt.df[pm.bin,"twoGAS"]
          fourGAS.exh.pm25.wt  <- pm25.wt.df[pm.bin,"fourGAS"]
          nonagDSL.exh.pm25.wt <- pm25.wt.df[pm.bin,"nonagDSL"]
          agDSL.exh.pm25.wt    <- pm25.wt.df[pm.bin,"agDSL"]

          marGAS.exh.pm25.df   <- emis.marGAS.df   * (marGAS.exh.pm25.wt / 100)
          twoGAS.exh.pm25.df   <- emis.twoGAS.df   * (twoGAS.exh.pm25.wt / 100)
          fourGAS.exh.pm25.df  <- emis.fourGAS.df  * (fourGAS.exh.pm25.wt / 100)
          nonagDSL.exh.pm25.df <- emis.nonagDSL.df * (nonagDSL.exh.pm25.wt / 100)
          agDSL.exh.pm25.df    <- emis.agDSL.df    * (agDSL.exh.pm25.wt / 100)

          tot.exh.pm25.df <- marGAS.exh.pm25.df + twoGAS.exh.pm25.df + fourGAS.exh.pm25.df + nonagDSL.exh.pm25.df + agDSL.exh.pm25.df

          if(output.flag & pollutant.name == "PM25_PRI_EXH")
          {
            setwd(paste(output.dir,"/",dow.name,sep=""))

            #Reorder HR00 to HR24
            pm24.df <- tot.exh.pm25.df[,c("dayav")]
            pm24.df <- cbind(pm24.df, tot.exh.pm25.df[,c(3:ncol(tot.exh.pm25.df))])
            pm24.df <- cbind(pm24.df, tot.exh.pm25.df[,c("HR00")])
            colnames(pm24.df)[1] <- c("dayav")
            colnames(pm24.df)[ncol(pm24.df)] <- c("HR24")

            #pm25.exh.out <- cbind(domain.sp[,c("Id","Row","Col","LON","LAT","TZ","URB_RUR")], tot.exh.pm25.df)
            pm25.exh.out  <- cbind(domain.sp[,c("Id","Row","Col","LON","LAT","TZ","URB_RUR")], pm24.df)
            pm25.exh.file <- paste("Offroad","_",pm.bin.name,"_EXH.rds",sep="")

            #write.csv(pm25.exh.out, file = pm25.exh.file, row.names = FALSE)
            saveRDS(pm25.exh.out, file=pm25.exh.file)
            setwd(five.dir)

          } #End pm2.5 output if statement
        } #End pm2.5 exhaust if statement

      } #End pm2.5 bin for loop

    } #End speciate pm2.5 if statement

  } #End day-of-week for loop

} #End pollutant for loop

file.rename(paste(five.dir.wd,'/','Offroad_PM25_PRI_EXH.rds',sep=''), paste(five.dir.wd,'/','Offroad_PM25-PRI_EXH.rds',sep=''))
file.rename(paste(five.dir.sa,'/','Offroad_PM25_PRI_EXH.rds',sep=''), paste(five.dir.sa,'/','Offroad_PM25-PRI_EXH.rds',sep=''))
file.rename(paste(five.dir.su,'/','Offroad_PM25_PRI_EXH.rds',sep=''), paste(five.dir.su,'/','Offroad_PM25-PRI_EXH.rds',sep=''))

file.rename(paste(five.dir.wd,'/','Offroad_PM10_PRI_EXH.rds',sep=''), paste(five.dir.wd,'/','Offroad_PM10-PRI_EXH.rds',sep=''))
file.rename(paste(five.dir.sa,'/','Offroad_PM10_PRI_EXH.rds',sep=''), paste(five.dir.sa,'/','Offroad_PM10-PRI_EXH.rds',sep=''))
file.rename(paste(five.dir.su,'/','Offroad_PM10_PRI_EXH.rds',sep=''), paste(five.dir.su,'/','Offroad_PM10-PRI_EXH.rds',sep=''))

print(paste("Finished processing all pollutants.",sep=""))
print(proc.time()-ptm)
