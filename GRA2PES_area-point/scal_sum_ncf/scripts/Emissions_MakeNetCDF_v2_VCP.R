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

base.dir     <- paste("/wrk/users/charkins/emissions/V7_GRA2PES/VCP",substr(year,3,4),sep="")
out.dir      <- paste("/wrk/users/charkins/emissions/V7_GRA2PES/VCP",substr(year,3,4),"_ncf",sep="")


sector.list  <- c("TotlAreaVCP")

#Set template NetCDF file
# domain.dir   <- "/wrk/csd4/bmcdonald/COVID_Modeling/CSD_Emissions/conus4k"
domain.dir   <- "/wrk/d2/charkins/templates/with_additional_VCPs/"
domain.files <- c("template_new_00to12Z.nc","template_new_12to24Z.nc")

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
# var.list     <- c(crit.list)

#--------------------------------------------------------
# Read wrfinput File
#--------------------------------------------------------
setwd(domain.dir)

wrf00z.nc <- nc_open(domain.files[1])
wrf12z.nc <- nc_open(domain.files[2])

#Get dimensions
n.col     <- max(wrf00z.nc$dim$west_east$vals)
n.row     <- max(wrf00z.nc$dim$south_north$vals)

#Get attributes
glob00z.attr <- ncatt_get(wrf00z.nc, varid=0)
time00z.attr <- ncatt_get(wrf00z.nc, varid="Times")

glob12z.attr <- ncatt_get(wrf12z.nc, varid=0)
time12z.attr <- ncatt_get(wrf12z.nc, varid="Times")

# nc_close(wrf00z.nc)
# nc_close(wrf12z.nc)

#--------------------------------------------------------
# Processing Emissions
#--------------------------------------------------------
#Cycle by month
for(m in 1:length(month.list))
{
  
  month <- month.list[m]
  
  #Cycle by sector
  for(s in 1:length(sector.list))
  {
    
    sector <- sector.list[s]
    
    #Cycle by day-of-week
    for(d in 1:length(dow.list))
    {
      
      dow <- dow.list[d]
      
      setwd(paste(base.dir,"/Month",sprintf("%02d",month),"/",sector,"/",dow,sep=""))
      
      print(paste("Copying 00to12Z file for ",sector," on month ",month," ",dow,"...",sep=""))
      
      out00z.file <- paste(sector,"_00to12Z.nc",sep="")
      out12z.file <- paste(sector,"_12to24Z.nc",sep="")
      
      ncks_command = paste("ncks --overwrite -x -v ",paste(var.list,collapse=",")," ",
                           paste(domain.dir,"/",domain.files[1],sep="")," ",
                           paste(out.dir,"/Month",sprintf("%02d",month),"/",sector,"/",dow,"/",out00z.file,sep=""),sep="")
      # file.copy(from = paste(domain.dir,"/",domain.files[1],sep=""), 
      #           to = paste(out.dir,"/Month",sprintf("%02d",month),"/",sector,"/",dow,"/",out00z.file,sep=""), 
      #           overwrite = TRUE)
      
      system(ncks_command,wait=TRUE)
      
      print(paste("Copying 12to24Z file for ",sector," on month ",month," ",dow,"...",sep=""))
      
      ncks_command = paste("ncks --overwrite -x -v ",paste(var.list,collapse=",")," ",
                           paste(domain.dir,"/",domain.files[2],sep="")," ",
                           paste(out.dir,"/Month",sprintf("%02d",month),"/",sector,"/",dow,"/",out12z.file,sep=""),sep="")
      # file.copy(from = paste(domain.dir,"/",domain.files[1],sep=""), 
      #           to = paste(out.dir,"/Month",sprintf("%02d",month),"/",sector,"/",dow,"/",out12z.file,sep=""), 
      #           overwrite = TRUE)
      system(ncks_command,wait=TRUE)
      
      out00z.nc <- nc_open(paste(out.dir,"/Month",sprintf("%02d",month),"/",sector,"/",dow,"/",out00z.file,sep=""), write = TRUE)
      out12z.nc <- nc_open(paste(out.dir,"/Month",sprintf("%02d",month),"/",sector,"/",dow,"/",out12z.file,sep=""), write = TRUE)
      
      dim.lon   <- ncdim_def("west_east", units="", vals=out00z.nc$dim$west_east$vals, unlim=FALSE, create_dimvar=FALSE)
      dim.lat   <- ncdim_def("south_north", units="", vals=out00z.nc$dim$south_north$vals, unlim=FALSE, create_dimvar=FALSE)
      dim.time  <- ncdim_def("Time", units="", vals=c(1:12), unlim=FALSE, create_dimvar=FALSE)
      
      #Cycle by variables
      for(v in 1:length(var.list))
      {

        pollutant <- var.list[v]
        
        print(paste("Processing ",sector," ",pollutant," emissions for month ",month," on ",dow,"...",sep=""))

        #--------------------------------------------------------
        # Input Emissions to wrfchemi
        #--------------------------------------------------------
        var.file <- grep(pattern = paste("_",pollutant,".rds",sep=""), list.files(), value = TRUE)
        
        # if the file exists open it, otherwise find the first file available and multiply values by zero
        if (length(var.file)>0){
          input.df <- readRDS(var.file)
          hr.cols <- grep("HR",colnames(input.df))
          var00z.df <- input.df[,hr.cols[1:12]]
          var12z.df <- input.df[,hr.cols[13:24]]
        } else {
          print(paste("WARNING: For sector: ",sector," Pollutant: ",pollutant,", no .rds file exists. Setting values to zero. ",sep=""))
          var.file <- grep(pattern = ".rds", list.files(), value = TRUE)
          input.df <- readRDS(var.file[1])
          hr.cols <- grep("HR",colnames(input.df))
          var00z.df <- input.df[,hr.cols[1:12]]*0
          var12z.df <- input.df[,hr.cols[13:24]]*0
        }
        
        
        
        
        
        #### Pull the original attributes from the variable in the file
        attrs= ncatt_get(wrf00z.nc,varid=pollutant)
        FieldType = attrs$FieldType
        MemoryOrder = attrs$MemoryOrder
        description = attrs$description
        units = attrs$units
        stagger = attrs$stagger
        grid_mapping = attrs$grid_mapping
        coordinates = attrs$coordinates
        
        temp.var <- ncvar_def(pollutant, units = "", dim=list(dim.lon,dim.lat,dim.time),chunksizes=c(400,400,1),compression=1, prec = "float")
        out00z.nc <- ncvar_add(out00z.nc, v = temp.var)
        out12z.nc <- ncvar_add(out12z.nc, v = temp.var)

        for(h in 1:24)
        {
           

          if(h <= 12)
          {
            mat = array(var00z.df[,h],dim=c(n.col,n.row))
            # empty[input.df['Col'],input.df['Row']] = var12z.df[,h]
            ncvar_put(out00z.nc, varid = pollutant, vals = mat, start=c(1,1,h), count=c(n.col,n.row,1))

          }else{
            mat = array(var12z.df[,c(h-12)],dim=c(n.col,n.row))
            ncvar_put(out12z.nc, varid = pollutant, vals = mat, start=c(1,1,h-12), count=c(n.col,n.row,1))

          }

        } #End hour loop
        
        # Add the variable attributes back to the new variable
        ncatt_put(out00z.nc, varid = pollutant, attname = "FieldType", attval = FieldType)
        ncatt_put(out00z.nc, varid = pollutant, attname = "MemoryOrder", attval = MemoryOrder)
        ncatt_put(out00z.nc, varid = pollutant, attname = "description", attval = description)
        ncatt_put(out00z.nc, varid = pollutant, attname = "units", attval = units)
        ncatt_put(out00z.nc, varid = pollutant, attname = "stagger", attval =   stagger)
        ncatt_put(out00z.nc, varid = pollutant, attname = "grid_mapping", attval = grid_mapping)
        ncatt_put(out00z.nc, varid = pollutant, attname = "coordinates", attval = coordinates)
        
        ncatt_put(out12z.nc, varid = pollutant, attname = "FieldType", attval = FieldType)
        ncatt_put(out12z.nc, varid = pollutant, attname = "MemoryOrder", attval = MemoryOrder)
        ncatt_put(out12z.nc, varid = pollutant, attname = "description", attval = description)
        ncatt_put(out12z.nc, varid = pollutant, attname = "units", attval = units)
        ncatt_put(out12z.nc, varid = pollutant, attname = "stagger", attval =   stagger)
        ncatt_put(out12z.nc, varid = pollutant, attname = "grid_mapping", attval = grid_mapping)
        ncatt_put(out12z.nc, varid = pollutant, attname = "coordinates", attval = coordinates)
            
        #End row loop

      } #End variable for loop
      
      #Re-do timestamp
      times00z <- ncvar_get(out00z.nc, varid = "Times")
      times12z <- ncvar_get(out12z.nc, varid = "Times")
      
      date.str   <- rep(paste(year,"-",sprintf("%02d",month),"-15",sep=""),length(times00z))
      hour00.str <- rep(paste("_",sprintf("%02d",0:11),":00:00",sep=""),length(times00z))
      hour12.str <- rep(paste("_",sprintf("%02d",12:23),":00:00",sep=""),length(times00z))
      
      date00.str <- paste(date.str, hour00.str, sep="")
      date12.str <- paste(date.str, hour12.str, sep="")
      
      ncvar_put(out00z.nc, varid = "Times", vals = date00.str)
      ncvar_put(out12z.nc, varid = "Times", vals = date12.str)
      
      nc_close(out00z.nc)
      nc_close(out12z.nc)
      
    } #End day-of-week loop
    
  } #End sector loop
  
} #End month loop

nc_close(wrf00z.nc)
nc_close(wrf12z.nc)


