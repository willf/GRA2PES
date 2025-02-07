rm(list=ls()) #Clear memory
ptm <- proc.time()
library(foreign)
#library(rgdal)
#library(sp)
#library(maptools)

#--------------------------------------------------------
# Set input variables
#--------------------------------------------------------
#Select year, month, and day-of-week
year         <- 2021
month.list   <- c(1,2,3,4,5,6,7,8,9,10,11,12)
dow.list     <- c("weekdy","satdy","sundy")

#Set directories
input.dir <- "/wrk/charkins/emissions/GRA2PES/V7_NRT_scaling/VCP21_202404/input"
output.dir   <- "/wrk/users/charkins/emissions/V7_GRA2PES/VCP21"


out.name     <- "AreaVCP"

#Set ArcGIS shapefile with domain
domain.dir   <- "/wrk/csd4/bmcdonald/COVID_Modeling/Domains"
domain.name  <- "nei04k_domain"

#Set ArcGIS shapefile with population data
#pop.dir      <- "/wrk/csd4/bmcdonald/COVID_Modeling/Domains"
pop.dir      <- "/wrk/charkins/Domains"
pop.shp      <- "Tract_2010Census_DP1_WGS_us04k_DISSOLVE"
pop.name     <- paste("pop",year,".rds",sep="")

#Population & GDP scaling factors
scale.file   <- "pop_factors.csv"
gdp.file     <- paste(out.name,"_McD_monthly.csv",sep="")

#VCP use and emission factors
sector.list  <- c("Pesticides","ConsCoatings","IndCoatings","Inks","ConsAdhesives","IndAdhesives","Personal","Cleaning")
sector.ind   <- c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE)

voc.use      <- c(4.3, 21.5, 0.74*23.3, 0.67*9.5, 11.6, 0.79*23.5, 27.9, 158)/1000 #kg/person/d (US)

voc.ef.us    <- c(370, 1.0*450*1.00, 410*1.00, 200*1.00, 410*1.00, 410*1.00, 400, 34) #g/kg (US)
voc.ef.ca    <- c(370, 0.8*450*0.65, 410*0.65, 200*1.00, 410*0.50, 410*0.70, 400, 34) #g/kg (CA)
voc.ef.la    <- c(370, 0.8*450*0.45, 410*0.45, 200*0.33, 410*0.50, 410*0.70, 400, 34) #g/kg (LA)
voc.ef.ny    <- c(370, 0.8*450*0.80, 410*0.80, 200*1.00, 410*0.60, 410*0.90, 400, 34) #g/kg (NY)
voc.ef.nj    <- c(370, 0.8*450*0.80, 410*0.80, 200*1.00, 410*0.60, 410*0.90, 400, 34) #g/kg (NJ)
voc.ef.ct    <- c(370, 0.8*450*0.80, 410*0.80, 200*1.00, 410*0.60, 410*0.90, 400, 34) #g/kg (CT)
voc.ef.nyc   <- c(370, 0.8*450*0.75, 410*0.75, 200*0.33, 410*0.60, 410*0.90, 400, 34) #g/kg (NYC)

#Notes:
#1. Exclude ptVCP use from industrial coatings, inks, and adhesives
#2. Account for 20% reduction in VOCs in coatings in AIM areas (based on ARB report)
#3. Printing ink shops in LA and NYC assumed to have maximum controls
#4. Account for CA, NY, LA, and NYC regulations on coatings
#5. Account for CA and NY regulations on adhesives & sealants

#Speciation files (number of speciation profiles must match number of sectors)
voc.wt.file  <- "areaVCP_wt.csv"
voc.mw.file  <- "areaVCP_mw.csv"

#Diurnal files (number of diurnal profiles must match number of sectors)
diurnal.file   <- "areaVCP_diurnal.csv" #local time
timezone.list  <- c("Atlantic","Eastern","Central","Mountain","Pacific","Alaska")
utc.list       <- c(-3,-4,-5,-6,-7,-8)

#Set flags to turn on functions
read.domain.flag <- FALSE
read.pop.flag    <- TRUE
speciate.flag    <- TRUE
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
# Process Population Data
#--------------------------------------------------------
if(read.pop.flag)
{
  setwd(pop.dir)

  print("Reading population file...")
  
  #------------------------------------
  # Read 2010 Population Data
  #------------------------------------
  pop.sp  <- read.dbf(paste(pop.shp,".dbf",sep=""))
  
  #------------------------------------
  # Calculate 2018 Population Data
  #------------------------------------
  #Get 2010 to 2018 scaling factors
  setwd(input.dir)
  
  scale.df <- read.csv(scale.file, header = TRUE)

  #Add 2010 population data to domain
  setwd(pop.dir)
  
  domain.sp$Id <- paste(domain.sp$Row, domain.sp$Col, sep=",")
  pop.sp$Id    <- paste(pop.sp$Row, pop.sp$Col, sep=",")
  
  domain.sp$POP10  <- 0
  domain.sp$POPAdj <- 0
  domain.sp$MANF   <- 0
  
  max.row <- max(domain.sp$Row) + 1
  max.col <- max(domain.sp$Col) + 1
  
  #Cycle through each role of dissolved population layer
  for(r in 1:nrow(pop.sp))
  {
    
    if(r %% 100 == 0)
    {
      print(paste("Processing grid ",r," of ",nrow(pop.sp),"...",sep=""))    
    }
    
    #Match grid cell in domain with population grid cell
    curr.row  <- pop.sp$Row[r]
    curr.col  <- pop.sp$Col[r]
    curr.cell <- curr.row * max.col + (curr.col + 1)
    
    #Add population data to domain
    domain.sp$POP10[curr.cell] <- pop.sp$POP10[r]
    
  }
  
  #Scale to 2018 population by state
  for(i in 1:nrow(scale.df))
  {
    state.name  <- scale.df$Name[i]
    pop.factor  <- scale.df$Factor[i]
    manf.factor <- scale.df$MANF[i]
    
    if(state.name != "US")
    {
      
      state.pos <- grep(state.name, domain.sp$STATE_NAME)
      
    }else{
      
      state.pos <- which(is.na(domain.sp$STATE_NAME))
      
    }
    
    domain.sp$POPAdj[state.pos] <- domain.sp$POP10[state.pos] * pop.factor
    domain.sp$MANF[state.pos]   <- manf.factor #Adjustment for manufacturing GDP by state
  }
  
  #------------------------------------
  # Save Updated Population Data
  #------------------------------------
  setwd(pop.dir)
  
  saveRDS(domain.sp, file = pop.name)  
   
}else{
  
  #Read population polygon shapefile
  setwd(pop.dir)
  domain.sp <- readRDS(pop.name)
  
}

#--------------------------------------------------------
# Generate non-Agriculture VCP Emissions
#--------------------------------------------------------
#Cycle by month
for(m in 1:length(month.list))
{
  
  month <- month.list[m]

  #Cycle by day-of-week
  for(d in 1:length(dow.list))
  {
    
    dow <- dow.list[d]
    
    setwd(input.dir)
    
    us.indices  <- grep("",    domain.sp$URBAN)
    ca.indices  <- grep("CA",  domain.sp$URBAN)
    la.indices  <- grep("LA",  domain.sp$URBAN)
    ny.indices  <- grep("NY",  domain.sp$URBAN)
    nj.indices  <- grep("NJ",  domain.sp$URBAN)
    ct.indices  <- grep("CT",  domain.sp$URBAN)
    nyc.indices <- grep("NYC", domain.sp$URBAN)
    
    voc.wt.df   <- read.csv(voc.wt.file, header = TRUE)
    voc.mw.df   <- read.csv(voc.mw.file, header = TRUE)
    diurnal.df  <- read.csv(diurnal.file, header = TRUE) 
    
    emis.df     <- data.frame(matrix(0, nrow = nrow(domain.sp), ncol = 25))
    colnames(emis.df) <- c("dayav","HR01","HR02","HR03","HR04","HR05","HR06",
                           "HR07","HR08","HR09","HR10","HR11","HR12",
                           "HR13","HR14","HR15","HR16","HR17","HR18",
                           "HR19","HR20","HR21","HR22","HR23","HR24")
    
    sum.emis.df <- emis.df
    
    #Cycle through non-agricultural sectors
    for(s in 1:length(sector.list))
    {
      sector   <- sector.list[s]
      ind.flag <- sector.ind[s]
      
      #Scale by monthly activity data
      gdp.data  <- read.csv(paste(input.dir,"/",gdp.file,sep=""))
      gdp.col   <- grep(paste("^",sector,"$",sep=""), colnames(gdp.data))
      gdp.value <- gdp.data[month,gdp.col]
      
      gdp.factors <- rep(gdp.value,nrow(emis.df))
      
      #------------------------------------
      # Estimate Daily Emissions
      #------------------------------------
      if(!ind.flag)
      {
        
        #Calculate US emissions first, then regional emissions
        #convert from g/d to metric tons/d
        daily.emis              <- (domain.sp$POPAdj              * voc.use[s] * gdp.factors * voc.ef.us[s])  / 10^6  
        daily.emis[ca.indices]  <- (domain.sp$POPAdj[ca.indices]  * voc.use[s] * gdp.factors * voc.ef.ca[s])  / 10^6  
        daily.emis[la.indices]  <- (domain.sp$POPAdj[la.indices]  * voc.use[s] * gdp.factors * voc.ef.la[s])  / 10^6
        daily.emis[ny.indices]  <- (domain.sp$POPAdj[ny.indices]  * voc.use[s] * gdp.factors * voc.ef.ny[s])  / 10^6
        daily.emis[nj.indices]  <- (domain.sp$POPAdj[nj.indices]  * voc.use[s] * gdp.factors * voc.ef.nj[s])  / 10^6
        daily.emis[ct.indices]  <- (domain.sp$POPAdj[ct.indices]  * voc.use[s] * gdp.factors * voc.ef.ct[s])  / 10^6
        daily.emis[nyc.indices] <- (domain.sp$POPAdj[nyc.indices] * voc.use[s] * gdp.factors * voc.ef.nyc[s]) / 10^6
        
      }else{
        
        #If industrial source, adjust by manufacturing GDP factor
        daily.emis              <- (domain.sp$POPAdj              * domain.sp$MANF              * voc.use[s] * gdp.factors * voc.ef.us[s])  / 10^6 
        daily.emis[ca.indices]  <- (domain.sp$POPAdj[ca.indices]  * domain.sp$MANF[ca.indices]  * voc.use[s] * gdp.factors * voc.ef.ca[s])  / 10^6 
        daily.emis[la.indices]  <- (domain.sp$POPAdj[la.indices]  * domain.sp$MANF[la.indices]  * voc.use[s] * gdp.factors * voc.ef.la[s])  / 10^6
        daily.emis[ny.indices]  <- (domain.sp$POPAdj[ny.indices]  * domain.sp$MANF[ny.indices]  * voc.use[s] * gdp.factors * voc.ef.ny[s])  / 10^6
        daily.emis[nj.indices]  <- (domain.sp$POPAdj[nj.indices]  * domain.sp$MANF[nj.indices]  * voc.use[s] * gdp.factors * voc.ef.nj[s])  / 10^6
        daily.emis[ct.indices]  <- (domain.sp$POPAdj[ct.indices]  * domain.sp$MANF[ct.indices]  * voc.use[s] * gdp.factors * voc.ef.ct[s])  / 10^6
        daily.emis[nyc.indices] <- (domain.sp$POPAdj[nyc.indices] * domain.sp$MANF[nyc.indices] * voc.use[s] * gdp.factors * voc.ef.nyc[s]) / 10^6
        
      }
      
      emis.df$dayav     <- daily.emis
      sum.emis.df$dayav <- sum.emis.df$dayav + daily.emis
      
      #------------------------------------
      # Temporal Allocation
      #------------------------------------
      #Cycle by timezone
      for(t in 1:length(timezone.list))
      {
        
        tz.name   <- timezone.list[t]
        tz.offset <- -utc.list[t]
        
        print(paste("Performing temporal allocation for ",sector," in ",tz.name," timezone...",sep=""))
        
        #Create diurnal factor table that accounts for different timezones
        diurnal.utc.df <- diurnal.df
        
        for(h in 1:24)
        {
          #Add hours to convert from local time to UTC time
          if(h+tz.offset <= 24)
          {
            
            diurnal.utc.df[h+tz.offset, 2:ncol(diurnal.df)] <- diurnal.df[h, 2:ncol(diurnal.df)]
            
          }else{
            
            diurnal.utc.df[h+tz.offset-24, 2:ncol(diurnal.df)] <- diurnal.df[h, 2:ncol(diurnal.df)]
          }
          
        } #End hour for loop
        
        #Multiply daily average emissions to hourly emissions
        tz.indices     <- grep(tz.name, domain.sp$TZ)            #select by timezone
        sector.col     <- grep(sector, colnames(diurnal.utc.df)) #select by sector
        diurnal.utc    <- diurnal.utc.df[,sector.col]            #select by sector
        
        daily.matrix   <- matrix(rep(daily.emis[tz.indices],each=24), ncol=24, byrow=TRUE)
        diurnal.matrix <- matrix(rep(diurnal.utc,each=nrow(daily.matrix)), nrow=nrow(daily.matrix))
        
        emis.df[tz.indices,2:25] <- (daily.matrix * diurnal.matrix)
        sum.emis.df[tz.indices,2:25] <- sum.emis.df[tz.indices,2:25] + (daily.matrix * diurnal.matrix)
        
      } #End time zone for loop
      
      #------------------------------------
      # Speciate VOCs
      #------------------------------------
      if(speciate.flag){
        
        #Cycle through each hydrocarbon bin
        for(bin in 1:nrow(voc.wt.df))
        {
          
          bin.name <- voc.wt.df$Bin[bin]
          
          print(paste("Processing ",bin.name," for month ",month," ",dow," ",sector,"...",sep=""))
          
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
            sector.voc.df <- emis.df * (sector.voc.wt / 100) * 10^6 / sector.voc.mw
            
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
      
    } #End sector for loop
    
    #---------------------------------------
    # Write Total VOC Output
    #---------------------------------------
    month.dir <- paste(output.dir,"/","Month",sprintf("%02d",month),sep="")
    
    if(output.flag)
    {
      setwd(month.dir)
      
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
      
      #Save coordinate variables (ID, Row, Col, LON, LAT, TZ)
      coord.df <- domain.sp[,1:6]
      
      #Output total mass emissions (in metric tons)
      print(paste("Writing output for sum of VOC...", sep=""))
      
      emis.out <- cbind(coord.df, sum.emis.df)
      
      if(out.rds.flag)
      {
        
        saveRDS(emis.out, file = paste(out.name,"_VOC.rds",sep=""))
        
      }else{
      
        write.csv(emis.out, file = paste(out.name,"_VOC.csv",sep=""),
                  row.names = FALSE)
        
      }
      
    } #End output flag if statement
    
    #---------------------------------------
    # Write Speciated VOC Emissions Output
    #---------------------------------------    
    #Output speciated VOC emissions (in moles)
    if(output.flag & speciate.flag)
    {
      setwd(paste(month.dir,"/",out.name,"/",dow,sep=""))
      
      #Save coordinate variables (ID, Row, Col, LON, LAT, TZ)
      coord.df <- domain.sp[,1:6]
      
      #Cycle through each speciated bin
      for(bin in 1:nrow(voc.wt.df))
      {
        voc.name <- voc.wt.df$Bin[bin]
        
        print(paste("Writing output for sum of ", voc.name,"...", sep=""))
        
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
    
  } #End dow for loop
  
} #End month for loop
