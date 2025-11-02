

# this script estimates RRDs and weighted average RRDs using hourly records from five MOPEX sites


setwd("~/ERRA introduction") # working direcory -- change as needed!

rm(list=ls())  #clear environment

library(data.table)
library(dplyr)


source("ERRA_v1.0.R")  # ensemble rainfall runoff analysis script






dat <- fread("MOPEX AnacostiaR hourly.txt", header=TRUE, sep="\t", na.strings=c("NA",".","","#N/A"))    #this is the test data file


# $ date_time: chr  "1990/10/02 00:00:00" "1990/10/02 01:00:00" "1990/10/02 02:00:00" "1990/10/02 03:00:00" ...
# $ p        : num  0 0 0 0 0 0 0 0 0 0 ...
# $ q        : num  0.0156 0.0156 0.0151 0.0156 0.0146 ...
# $ PET      : num  0.11 0.11 0.11 0.11 0.11 ...
# $ yr       : int  1990 1990 1990 1990 1990 1990 1990 1990 1990 1990 ...
# $ mm       : int  10 10 10 10 10 10 10 10 10 10 ...
# $ dd       : int  2 2 2 2 2 2 2 2 2 2 ...
# $ PET_mm/d : num  2.65 2.65 2.65 2.65 2.65 ...
# $ year     : num  1991 1991 1991 1991 1991 ...
# $ Tmax     : num  22.5 22.5 22.5 22.5 22.5 ...
# $ Tmin     : num  9.1 9.1 9.1 9.1 9.1 9.1 9.1 9.1 9.1 9.1 ...



fileID <- "Anacostia_hourly_Fig10ab"

p <- dat$p        
q <- dat$q

zz <- ERRA(p=p, q=q, m=100)       # run ERRA without nonlinear analysis
ww <- ERRA(p=p, q=q, m=100, xknots=c(6,20), xknot_type="even")       # run ERRA with nonlinear analysis

with(ww, {                                          # write the output to files
  fwrite(cbind(zz$RRD, wtd_avg_RRD[,-1]), paste0(fileID, "_RRDs_", options, ".txt"), sep="\t")  # combine RRD and weighted_average_RRD
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
})





# Now re-do this analysis for the other sites

dat <- fread("MOPEX NantahalaR hourly.txt", header=TRUE, sep="\t", na.strings=c("NA",".","","#N/A"))    #this is the test data file
fileID <- "Nantahala_hourly_Fig10cd"

p <- dat$p        
q <- dat$q

zz <- ERRA(p=p, q=q, m=100)       # run ERRA without nonlinear analysis
ww <- ERRA(p=p, q=q, m=100, xknots=c(6,20), xknot_type="even")       # run ERRA with nonlinear analysis

with(ww, {                                          # write the output to files
  fwrite(cbind(zz$RRD, wtd_avg_RRD[,-1]), paste0(fileID, "_RRDs_", options, ".txt"), sep="\t")  # combine RRD and weighted_average_RRD
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
})



dat <- fread("MOPEX SFNewR hourly.txt", header=TRUE, sep="\t", na.strings=c("NA",".","","#N/A"))    #this is the test data file
fileID <- "SFNew_hourly_Fig10ef"

p <- dat$p        
q <- dat$q

zz <- ERRA(p=p, q=q, m=100, agg=3)       # run ERRA without nonlinear analysis
ww <- ERRA(p=p, q=q, m=100, agg=3, xknots=c(6,20), xknot_type="even")       # run ERRA with nonlinear analysis

with(ww, {                                          # write the output to files
  fwrite(cbind(zz$RRD, wtd_avg_RRD[,-1]), paste0(fileID, "_RRDs_", options, ".txt"), sep="\t")  # combine RRD and weighted_average_RRD
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
})



dat <- fread("MOPEX FrenchBroadR hourly.txt", header=TRUE, sep="\t", na.strings=c("NA",".","","#N/A"))    #this is the test data file
fileID <- "FrenchBroad_hourly_Fig10gh"

p <- dat$p        
q <- dat$q

zz <- ERRA(p=p, q=q, m=100, agg=6)       # run ERRA without nonlinear analysis
ww <- ERRA(p=p, q=q, m=100, agg=6, xknots=c(6,20), xknot_type="even")       # run ERRA with nonlinear analysis

with(ww, {                                          # write the output to files
  fwrite(cbind(zz$RRD, wtd_avg_RRD[,-1]), paste0(fileID, "_RRDs_", options, ".txt"), sep="\t")  # combine RRD and weighted_average_RRD
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
})






dat <- fread("MOPEX ClinchR hourly.txt", header=TRUE, sep="\t", na.strings=c("NA",".","","#N/A"))    #this is the test data file
fileID <- "Clinch_hourly_Fig10ij"

p <- dat$p        
q <- dat$q

zz <- ERRA(p=p, q=q, m=100, agg=6)       # run ERRA without nonlinear analysis
ww <- ERRA(p=p, q=q, m=100, agg=6, xknots=c(6,20), xknot_type="even")       # run ERRA with nonlinear analysis

with(ww, {                                          # write the output to files
  fwrite(cbind(zz$RRD, wtd_avg_RRD[,-1]), paste0(fileID, "_RRDs_", options, ".txt"), sep="\t")  # combine RRD and weighted_average_RRD
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
})


