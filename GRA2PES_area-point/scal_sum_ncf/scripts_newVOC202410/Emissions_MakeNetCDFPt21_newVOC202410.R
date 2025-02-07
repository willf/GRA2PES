#--------------------------------------------------------
# bash
# eval "$(/home/bmcdonald/anaconda3/bin/conda shell.bash hook)"
# conda activate r-gdal
#--------------------------------------------------------
rm(list=ls()) #Clear memory
ptm <- proc.time()

#--------------------------------------------------------
# Emissions_Make_wrfchemi.R (created 8/16/2020)
#
# Creates wrfchemi file from onroad, nonroad, area .csv
# files and given wrfinput file. Also maps hydrocarbon
# bins in emission file to WRF-Chem VOC speciation.
#
# Input file units are in metric tons/d, expect for VOCs,
# which are in moles/d.
#
#--------------------------------------------------------
library(ncdf4)
library(foreign)

#--------------------------------------------------------
# Set input variables
#--------------------------------------------------------
#Select month and day-of-week
dow.list     <- c("weekdy","satdy","sundy")
month.list   <- c(1,2,3,4,5,6,7,8,9,10,11,12)
year         <- 2021

base.dir     <- paste("/wrk/users/charkins/emissions/V7_GRA2PES/POINT21",sep="")
out.dir      <- paste("/wrk/users/charkins/emissions/V7_GRA2PES/POINT21_ncf",sep="")

sector.list  <- c("TotlPoint_newVCPVOC202410")

#Set template NetCDF file
domain.dir   <- "/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_GRA2PES2021/Stu_RELPT/output_RELPT_template_CEMS2021"
domain.file  <- "Relpointinfo_Criteria_n_VCP.csv"

#Variables to process
crit.list    <- c("CO2","CO","NH3","NOX","PM10-PRI","PM25-PRI","SO2","VOC")
hc.list      <- c("HC01","HC02","HC03","HC04","HC05","HC06","HC07","HC08","HC09","HC10",
                  "HC11","HC12","HC13","HC14","HC15","HC16","HC17","HC18","HC19","HC20",
                  "HC21","HC22","HC23","HC24","HC25","HC26","HC27","HC28","HC29","HC30",
                  "HC31","HC32","HC33","HC34","HC35","HC36","HC37","HC38","HC39","HC40",
                  "HC41","HC42","HC43","HC44","HC45","HC46","HC47","HC48","HC49","HC50",
                  "HC51","HC52","HC53","HC54","HC55","HC56","HC57","HC58","HC59","HC60",
                  "HC61","HC62","HC63","HC64","HC65","HC66","HC67","HC68","HC69","HC70",
				  "HC71","HC72","HC73","HC74","HC75","HC76","HC77","HC78","HC79","HC80",
				  "HC81","HC82","HC83","HC84")
pm.list      <- c("PM01","PM02","PM03","PM04","PM05","PM06","PM07","PM08","PM09","PM10",
                  "PM11","PM12","PM13","PM14","PM15","PM16","PM17","PM18","PM19")

var.list     <- c(crit.list, hc.list, pm.list)

#Set flags
in.rds.flag  <- TRUE

#--------------------------------------------------------
# Read Relpoint File
#--------------------------------------------------------
setwd(domain.dir)

relpt <- read.csv(domain.file, header=TRUE)

n.pts <- nrow(relpt)

#--------------------------------------------------------
# Process Emissions
#--------------------------------------------------------
#Cycle by month
for(m in 1:length(month.list))
{

  month <- month.list[m]
  
  setwd(paste(base.dir,"/Month",sprintf("%02d",month),sep=""))
  
  #Cycle by sector
  for(s in 1:length(sector.list))
  {
    
    sector <- sector.list[s]
    
    setwd(paste(base.dir,"/Month",sprintf("%02d",month),"/",sector,sep=""))
    
    #Cycle by day-of-week
    for(d in 1:length(dow.list))
    {
      
      dow <- dow.list[d]
      
      setwd(paste(base.dir,"/Month",sprintf("%02d",month),"/",sector,"/",dow,sep=""))
      
      #--------------------------------------------------------
      # Create NetCDF with Stack Parameters
      #--------------------------------------------------------
      dim.row  <- ncdim_def("ROW",  units="", vals=c(1:n.pts), unlim=FALSE, create_dimvar=FALSE)
      dim.time <- ncdim_def("Time", units="", vals=c(1:12),    unlim=FALSE, create_dimvar=FALSE)
      #londim  <- ncdim_def("lon",  units="degrees_east",  as.double(relpt$Lon_deg))
      #latdim  <- ncdim_def("lat",  units="degrees_north", as.double(relpt$Lat_deg))
      #timedim <- ncdim_def("time", units="", c(1:12))
      
      ITYPE     <- ncvar_def("ITYPE",     units = "",     dim=list(dim.row), 0, prec = "float")
      STKht     <- ncvar_def("STKht",     units = "m",    dim=list(dim.row), 0, prec = "float")
      STKdiam   <- ncvar_def("STKdiam",   units = "m",    dim=list(dim.row), 0, prec = "float")
      STKtemp   <- ncvar_def("STKtemp",   units = "K",    dim=list(dim.row), 0, prec = "float")
      STKve     <- ncvar_def("STKve",     units = "m/s",  dim=list(dim.row), 0, prec = "float")
      STKflw    <- ncvar_def("STKflw",    units = "m3/s", dim=list(dim.row), 0, prec = "float")
      FUGht     <- ncvar_def("FUGht",     units = "m",    dim=list(dim.row), 0, prec = "float")
      XLONG     <- ncvar_def("XLONG", units = "degree_east", dim=list(dim.row), 0, prec = "float")
      XLAT      <- ncvar_def("XLAT", units = "degree_north", dim=list(dim.row), 0, prec = "float")
      
      wrf00z.file <- paste(out.dir,"/Month",sprintf("%02d",month),"/",sector,"/",dow,"/",sector,"_00to12Z.nc",sep="")
      wrf12z.file <- paste(out.dir,"/Month",sprintf("%02d",month),"/",sector,"/",dow,"/",sector,"_12to24Z.nc",sep="")
      
      #Delete file if exists
      if(file.exists(wrf00z.file))
      {
        
        print("Deleting existing 00to12Z file...")
        file.remove(wrf00z.file)
        
      }
      
      if(file.exists(wrf12z.file))
      {
        
        print("Deleting existing 12to24Z file...")
        file.remove(wrf12z.file)
        
      }
      
      wrf00z.nc <- nc_create(wrf00z.file, list(ITYPE, STKht, STKdiam, STKtemp, STKve, STKflw, FUGht, XLONG, XLAT), force_v4=TRUE)

      ncvar_put(wrf00z.nc, varid = "ITYPE", vals = relpt$ITYPE)
      ncatt_put(wrf00z.nc, varid = "ITYPE", attname = "FieldType",    attval = 104)
      ncatt_put(wrf00z.nc, varid = "ITYPE", attname = "MemoryOrder",  attval = "X")
      ncatt_put(wrf00z.nc, varid = "ITYPE", attname = "description",  attval = "Stack Type")
      ncatt_put(wrf00z.nc, varid = "ITYPE", attname = "units",        attval = "unitless")
      ncatt_put(wrf00z.nc, varid = "ITYPE", attname = "coordinates",  attval = "XLONG XLAT")
      
      ncvar_put(wrf00z.nc, varid = "STKht", vals = relpt$STKht_m)
      ncatt_put(wrf00z.nc, varid = "STKht", attname = "FieldType",    attval = 104)
      ncatt_put(wrf00z.nc, varid = "STKht", attname = "MemoryOrder",  attval = "X")
      ncatt_put(wrf00z.nc, varid = "STKht", attname = "description",  attval = "Stack Height")
      ncatt_put(wrf00z.nc, varid = "STKht", attname = "units",        attval = "m")
      ncatt_put(wrf00z.nc, varid = "STKht", attname = "coordinates",  attval = "XLONG XLAT")
      
      ncvar_put(wrf00z.nc, varid = "STKdiam", vals = relpt$STKdiam_m)
      ncatt_put(wrf00z.nc, varid = "STKdiam", attname = "FieldType",    attval = 104)
      ncatt_put(wrf00z.nc, varid = "STKdiam", attname = "MemoryOrder",  attval = "X")
      ncatt_put(wrf00z.nc, varid = "STKdiam", attname = "description",  attval = "Stack Diameter")
      ncatt_put(wrf00z.nc, varid = "STKdiam", attname = "units",        attval = "m")
      ncatt_put(wrf00z.nc, varid = "STKdiam", attname = "coordinates",  attval = "XLONG XLAT")

      ncvar_put(wrf00z.nc, varid = "STKtemp", vals = relpt$STKtemp_K)
      ncatt_put(wrf00z.nc, varid = "STKtemp", attname = "FieldType",    attval = 104)
      ncatt_put(wrf00z.nc, varid = "STKtemp", attname = "MemoryOrder",  attval = "X")
      ncatt_put(wrf00z.nc, varid = "STKtemp", attname = "description",  attval = "Stack Temperature")
      ncatt_put(wrf00z.nc, varid = "STKtemp", attname = "units",        attval = "K")
      ncatt_put(wrf00z.nc, varid = "STKtemp", attname = "coordinates",  attval = "XLONG XLAT")
            
      ncvar_put(wrf00z.nc, varid = "STKve", vals = relpt$STKve_mps)
      ncatt_put(wrf00z.nc, varid = "STKve", attname = "FieldType",    attval = 104)
      ncatt_put(wrf00z.nc, varid = "STKve", attname = "MemoryOrder",  attval = "XY")
      ncatt_put(wrf00z.nc, varid = "STKve", attname = "description",  attval = "Stack Velocity")
      ncatt_put(wrf00z.nc, varid = "STKve", attname = "units",        attval = "m/s")
      ncatt_put(wrf00z.nc, varid = "STKve", attname = "coordinates",  attval = "XLONG XLAT")
      
      ncvar_put(wrf00z.nc, varid = "STKflw", vals = relpt$STKflw_m3ps)
      ncatt_put(wrf00z.nc, varid = "STKflw", attname = "FieldType",    attval = 104)
      ncatt_put(wrf00z.nc, varid = "STKflw", attname = "MemoryOrder",  attval = "XY")
      ncatt_put(wrf00z.nc, varid = "STKflw", attname = "description",  attval = "Stack Flow")
      ncatt_put(wrf00z.nc, varid = "STKflw", attname = "units",        attval = "m3/s")
      ncatt_put(wrf00z.nc, varid = "STKflw", attname = "coordinates",  attval = "XLONG XLAT")
      
      ncvar_put(wrf00z.nc, varid = "FUGht", vals = relpt$FUGht_m)
      ncatt_put(wrf00z.nc, varid = "FUGht", attname = "FieldType",    attval = 104)
      ncatt_put(wrf00z.nc, varid = "FUGht", attname = "MemoryOrder",  attval = "XY")
      ncatt_put(wrf00z.nc, varid = "FUGht", attname = "description",  attval = "Fugitive Height")
      ncatt_put(wrf00z.nc, varid = "FUGht", attname = "units",        attval = "m")
      ncatt_put(wrf00z.nc, varid = "FUGht", attname = "coordinates",  attval = "XLONG XLAT")
      
      ncvar_put(wrf00z.nc, varid = "XLONG", vals = relpt$Lon_deg)
      ncatt_put(wrf00z.nc, varid = "XLONG", attname = "FieldType",    attval = 104)
      ncatt_put(wrf00z.nc, varid = "XLONG", attname = "MemoryOrder",  attval = "XY")
      ncatt_put(wrf00z.nc, varid = "XLONG", attname = "description",  attval = "LONGITUDE, WEST IS NEGATIVE")
      ncatt_put(wrf00z.nc, varid = "XLONG", attname = "units",        attval = "degree_east")
      ncatt_put(wrf00z.nc, varid = "XLONG", attname = "coordinates",  attval = "XLONG XLAT")
      
      ncvar_put(wrf00z.nc, varid = "XLAT", vals = relpt$Lat_deg)
      ncatt_put(wrf00z.nc, varid = "XLAT", attname = "FieldType",    attval = 104)
      ncatt_put(wrf00z.nc, varid = "XLAT", attname = "MemoryOrder",  attval = "XY")
      ncatt_put(wrf00z.nc, varid = "XLAT", attname = "description",  attval = "LATITUDE, SOUTH IS NEGATIVE")
      ncatt_put(wrf00z.nc, varid = "XLAT", attname = "units",        attval = "degree_north")
      ncatt_put(wrf00z.nc, varid = "XLAT", attname = "coordinates",  attval = "XLONG XLAT")
      
      nc_close(wrf00z.nc)
      
      file.copy(from=wrf00z.file, to=wrf12z.file)
      
      wrf00z.nc <- nc_open(wrf00z.file, write = TRUE)
      wrf12z.nc <- nc_open(wrf12z.file, write = TRUE)
      
      #Cycle by variable
      for(v in 1:length(var.list))
      {

        pollutant <- var.list[v]
        
        print(paste("Processing ",pollutant," for ",sector," in Month ",month," on ",dow,"...",sep=""))
        
        #--------------------------------------------------------
        # Read Input Files
        #--------------------------------------------------------
        if(in.rds.flag)
        {
          emis.file <- grep(paste(sector,"_",pollutant,".rds",sep=""),list.files(),value = TRUE)
		  if (length(emis.file)>0){
			emis.df <- readRDS(emis.file)
			hr.cols <- grep("HR",colnames(emis.df))
			emis00z.df <- emis.df[,hr.cols[1:12]]
			emis12z.df <- emis.df[,hr.cols[13:24]]
		  }else{
			print(paste("WARNING: For sector: ",sector," Pollutant: ",pollutant,", no .rds file exists. Setting values to zero. ",sep=""))
			emis.file <- grep(pattern = ".rds", list.files(), value = TRUE)
			emis.df <- readRDS(emis.file[1])
			hr.cols <- grep("HR",colnames(emis.df))
			emis00z.df <- emis.df[,hr.cols[1:12]]*0
			emis12z.df <- emis.df[,hr.cols[13:24]]*0
		  }
		  
        }else{

          emis.file <- grep(paste(sector,"_",pollutant,".csv",sep=""),list.files(),value = TRUE)
		  if (length(emis.file)>0){
			emis.df <- read.csv(emis.file, header=TRUE)
			hr.cols <- grep("HR",colnames(emis.df))
			emis00z.df <- emis.df[,hr.cols[1:12]]
			emis12z.df <- emis.df[,hr.cols[13:24]]
		  }else{
			print(paste("WARNING: For sector: ",sector," Pollutant: ",pollutant,", no .rds file exists. Setting values to zero. ",sep=""))
			emis.file <- grep(pattern = ".csv", list.files(), value = TRUE)
			emis.df <- read.csv(emis.file[1], header=TRUE)
			hr.cols <- grep("HR",colnames(emis.df))
			emis00z.df <- emis.df[,hr.cols[1:12]]*0
			emis12z.df <- emis.df[,hr.cols[13:24]]*0
		  }
          
        }
        
        #--------------------------------------------------------
        # Add Variable to NetCDF
        #--------------------------------------------------------
        if(substr(pollutant,1,2)=="HC")
        {
          
          unit.name = "mole hr^-1"
          
        }else{
          
          unit.name = "metric_Ton hr^-1"
          
        }
        
        #Initialize pollutant variable
        temp.var <- ncvar_def(pollutant, units = unit.name, dim=list(dim.row, dim.time), prec = "float")
        
        wrf00z.nc <- ncvar_add(wrf00z.nc, v = temp.var)
        ncatt_put(wrf00z.nc, varid = pollutant, attname = "FieldType",    attval = 104)
        ncatt_put(wrf00z.nc, varid = pollutant, attname = "MemoryOrder",  attval = "XY")
        ncatt_put(wrf00z.nc, varid = pollutant, attname = "description",  attval = pollutant)
        ncatt_put(wrf00z.nc, varid = pollutant, attname = "units",        attval = unit.name)
        ncatt_put(wrf00z.nc, varid = pollutant, attname = "coordinates",  attval = "XLONG XLAT")
        
        wrf12z.nc <- ncvar_add(wrf12z.nc, v = temp.var)
        ncatt_put(wrf12z.nc, varid = pollutant, attname = "FieldType",    attval = 104)
        ncatt_put(wrf12z.nc, varid = pollutant, attname = "MemoryOrder",  attval = "XY")
        ncatt_put(wrf12z.nc, varid = pollutant, attname = "description",  attval = pollutant)
        ncatt_put(wrf12z.nc, varid = pollutant, attname = "units",        attval = unit.name)
        ncatt_put(wrf12z.nc, varid = pollutant, attname = "coordinates",  attval = "XLONG XLAT")
        
        #Cycle by hour
        for(h in 1:24)
        {
          
          if(h <= 12)
          {
            
            ncvar_put(wrf00z.nc, varid = pollutant, vals = emis00z.df[,h], start=c(1,h), count=c(n.pts,1))
            
          }else{
            
            ncvar_put(wrf12z.nc, varid = pollutant, vals = emis12z.df[,h-12], start=c(1,h-12), count=c(n.pts,1))
                        
          }
          
        }

      } #End loop by variable
      
      #Re-do timestamp
      Times  <- ncvar_def("Times", units = "", dim=list(dim.time), prec = "char")
        
      wrf00z.nc <- ncvar_add(wrf00z.nc, v = Times)
      wrf12z.nc <- ncvar_add(wrf12z.nc, v = Times)
      
      date.str   <- rep(paste(year,"-",sprintf("%02d",month),"-15",sep=""),12)
      hour00.str <- rep(paste("_",sprintf("%02d",0:11),":00:00",sep=""),12)
      hour12.str <- rep(paste("_",sprintf("%02d",12:23),":00:00",sep=""),12)
      
      date00.str <- paste(date.str, hour00.str, sep="")
      date12.str <- paste(date.str, hour12.str, sep="")
      
      ncvar_put(wrf00z.nc, varid = "Times", vals = date00.str)
      ncvar_put(wrf12z.nc, varid = "Times", vals = date12.str)
      
      nc_close(wrf00z.nc)
      nc_close(wrf12z.nc)
            
    } #End loop by day-of-week
    
  } #End loop by sector

} #End loop by month

