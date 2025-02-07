
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
 dow.list     <- c("weekdy","satdy","sundy")
 month.list   <- c(1,7)
 
 base.dir     <- "/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_GRA2PES2021/Stu_RELPT/output_RELPT_emis_CEMS202102"
 out.dir      <- "/wrk/csd4/clyu/GHG_CO2/Improving_inventory/V7_GRA2PES2021/Stu_RELPT/output_RELPT_emis_CEMS202102"
 
 #Define sectors following CAMS
 egu.list     <- c("PtEGU_Coal","PtEGU_NG","PtEGU_Oil","PtEGU_BIO")
 comm.list    <- c("PtCOMM_Coal","PtCOMM_NG","PtCOMM_Oil","PtCOMM_BIO")
 indf.list    <- c("PtIND_Coal","PtIND_NG","PtIND_NG2","PtIND_Oil","PtIND_Oil2","PtIND_BIO")
 indp.list    <- c("PtCHEM","PtFOOD","PtMETAL","PtREFINE","PtPULP","PtELECT","PtMOTOR","PtAPPAREL","PtPHOTO","PtDRUG","PtMISC","PtMISC2")
 vcp.list     <- c("PtAdhesives","PtCoatings","PtDegreasing","PtInks")
 fug.list     <- c("PtOnG","PtSTORAGE","PtEVAPGAS","PtCONST")
 mob.list     <- c("PtAVIATION","PtRAIL")

 sector.list  <- c(egu.list) 

 #Variables to process
 crit.list    <- c("CO2","CO","NH3","NOX","PM10-PRI","PM25-PRI","SO2","VOC")
 hc.list      <- c("HC01","HC02","HC03","HC04","HC05","HC06","HC07","HC08","HC09","HC10",
                   "HC11","HC12","HC13","HC14","HC15","HC16","HC17","HC18","HC19","HC20",
                   "HC21","HC22","HC23","HC24","HC25","HC26","HC27","HC28","HC29","HC30",
                   "HC31","HC32","HC33","HC34","HC35","HC36","HC37","HC38","HC39","HC40",
                   "HC41","HC42","HC43","HC44","HC45","HC46","HC47","HC48","HC49","HC50")
 pm.list      <- c("PM01","PM02","PM03","PM04","PM05","PM06","PM07","PM08","PM09","PM10",
                   "PM11","PM12","PM13","PM14","PM15","PM16","PM17","PM18","PM19")
 
 var.list     <- c(crit.list, hc.list, pm.list)
 
 #Set flags to turn on functions
 output.flag     <- TRUE
 
 #--------------------------------------------------------
 # Summarize Emissions
 #--------------------------------------------------------
 #Cycle by sector
 for(s in 1:length(sector.list))
 {
 
   sector <- sector.list[s]
 
   #Initialize output
   weekdy.array <- rep(0, length(var.list))
   satdy.array  <- rep(0, length(var.list))
   sundy.array  <- rep(0, length(var.list))
 
   out.df <- data.frame(var.list, weekdy.array, satdy.array, sundy.array)
   colnames(out.df) <- c("Species","weekdy","satdy","sundy")
 
   #Cycle by variable
   for(v in 1:length(var.list))
   {
     pollutant <- var.list[v]
 
     #Cycle by day-of-week
     for(d in 1:length(dow.list))
     {
 
       dow <- dow.list[d]
       exist.count <- 0
       
       #Cycle by months
       for(m in 1:length(month.list))
       {
         
         month <- month.list[m]
         
         setwd(paste(base.dir,"/Month",sprintf("%02d",month),"/",sector,"/",dow,sep=""))
         
         print(paste("Reading ",pollutant," for ",sector," in month ",sprintf("%02d",month)," on ",dow,"...",sep=""))
         
         emis.file <- grep(paste(sector,"_",pollutant,"_",dow,".rds",sep=""),list.files(),value = TRUE)
         
         if(length(emis.file) > 0)
         {
           emis.df <- readRDS(emis.file)
           
           hr.cols <- grep("HR",colnames(emis.df))
           day.col <- grep("dayav",colnames(emis.df))
           
           if(m == 1)
           {
             
             sum.df <- emis.df[,c(day.col,hr.cols)]
             
           }else{
             
             sum.df <- sum.df + emis.df[,c(day.col,hr.cols)]
             
           }
           
           exist.count <- exist.count + 1
         }
         
       } #End month for loop
       
       if(exist.count > 0)
       {
         avg.df <- sum.df / length(month.list)
         
         info.col <- day.col - 1
         
         out.df <- cbind(emis.df[,c(1:info.col)], avg.df)
         
         #---------------------------------------
         # Output Avg
         #---------------------------------------
         if(output.flag)
         {
           
           month.dir <- paste(out.dir,"/Month00",sep="")
           
           #Set Month Directory
           if(dir.exists(month.dir))
           {
             
             setwd(month.dir)
             
           }else{
             
             dir.create(month.dir)
             setwd(month.dir)
             
           }
           
           #Set Sector Directory
           if(dir.exists(sector))
           {
             
             setwd(sector)
             
           }else{
             
             dir.create(sector)
             setwd(sector)
             
           }
           
           #Set Day-of-Week Directory
           if(dir.exists(dow))
           {
             
             setwd(dow)
             
           }else{
             
             dir.create(dow)
             setwd(dow)
             
           }
           
           print(paste("Writing ",pollutant," for ",sector," on ",dow,"...",sep=""))
           saveRDS(out.df, file=paste(sector,"_",pollutant,"_",dow,".rds",sep=""))
           
         } #End output flag
         
       } #End exist flag
 
     } #End day-of-week loop
 
   } #End pollutant for loop
 
 } #End sector for loop
