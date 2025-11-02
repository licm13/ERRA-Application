

# this script tests out ERRA on Roanoke River hourly records from USGS, 
# with precipitation records from Roanoke and Blacksburg airports

# this is the analysis behind Figure 4


setwd("~/ERRA introduction") # working directory -- change as needed!

rm(list=ls())  #clear environment

library(data.table)
library(dplyr)


source("ERRA_v1.0.R")  # ensemble rainfall runoff analysis script 



# read in data
dat <- fread("Roanoke_hourly_Q_and_airport_P_2006-2023.txt", header=TRUE, sep="\t", na.strings=c("NA",".","","#N/A"))    #this is the test data file


# $ date_time_UTC          : chr  "1990/10/01 05:00" "1990/10/01 06:00" "1990/10/01 07:00" "1990/10/01 08:00" ...
# $ date_time_EST          : chr  "1990/10/01 00:00" "1990/10/01 01:00" "1990/10/01 02:00" "1990/10/01 03:00" ...
# $ Blacksburg_METAR_P_mm.h: num  NA NA NA NA NA NA NA NA NA NA ...
# $ Roanoke_METAR_P_mm.h   : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Roanoke_Q_mm.h         : num  0.00758 0.00758 0.0079 0.00779 0.00747 ...




# Here we calculate the RRD using the average of Roanoke and Blacksburg precip
# This is the basis for Figure 4a

# average the two precipitation records
p <- 0.5*(dat$Roanoke_METAR_P_mm.h + dat$Blacksburg_METAR_P_mm.h) 

# discharge measured at Roanoke USGS gauge
q <- dat$Roanoke_Q_mm.h   


# now call ERRA
zz <- ERRA(p=p, q=q, m=168)

# now save the RRD and the peak statistics
fileID <- "Roanoke_avgP_Fig4a"
with(zz, {
  fwrite(RRD, paste0(fileID, "_RRD_", options, ".txt"), sep="\t") 
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
})




# Here we calculate the RRD using Roanoke precip only
# This is used in Figure 4b

fileID <- "RoanokeR_RoanokeP_Fig4b"

# use Roanoke precipitation
p <- dat$Roanoke_METAR_P_mm.h

# now call ERRA, and save the RRD and the peak statistics
zz <- ERRA(p=p, q=q, m=168)
with(zz, {
  fwrite(RRD, paste0(fileID, "_RRD_", options, ".txt"), sep="\t") 
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
})






# Here we calculate the RRD using Blacksburg precip only
# This is used in Figure 4b

fileID <- "RoanokeR_BlacksburgP_Fig4b"

# use Blacksburg precipitation
p <- dat$Blacksburg_METAR_P_mm.h

# now call ERRA, and save the RRD and the peak statistics
zz <- ERRA(p=p, q=q, m=168)
with(zz, {
  fwrite(RRD, paste0(fileID, "_RRD_", options, ".txt"), sep="\t") 
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
})









# Here we solve jointly for the RRDs of Roanoke and Blacksburg precipitation
# via joint deconvolution and demixing to separate their effects on streamflow.
# This is the basis of Figure 4c

fileID <- "RoanokeR_TwoP_Fig4c"

# precipitation input is now a matrix that includes both Roanoke and Blacksburg precipitation
p <- cbind(pR=0.5*dat$Roanoke_METAR_P_mm.h, pB=0.5*dat$Blacksburg_METAR_P_mm.h)        

# now call ERRA, and save the RRD and the peak statistics
zz <- ERRA(p=p, q=q, m=168)
with(zz, {
  fwrite(RRD, paste0(fileID, "_RRD_", options, ".txt"), sep="\t") 
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
})


