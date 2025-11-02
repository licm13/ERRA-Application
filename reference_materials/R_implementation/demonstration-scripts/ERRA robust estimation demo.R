

setwd("~/ERRA introduction")   # working directory -- change as needed!

rm(list=ls())  #clear environment

library(data.table)
library(dplyr)



source("ERRA_v1.0.R")  # ensemble rainfall runoff analysis script



dat <- fread("Plynlimon hourly data for ERRA.txt", header=TRUE, sep="\t", na.strings=c("NA",".","","#N/A"))    #this is the data file


# Classes 'data.table' and 'data.frame':	326545 obs. of  22 variables:
#   $ hour_ending                            : chr  "01/10/1973 00:00" "01/10/1973 01:00" "01/10/1973 02:00" "01/10/1973 03:00" ...
# $ year                                   : num  1974 1974 1974 1974 1974 ...
# $ yr                                     : int  1973 1973 1973 1973 1973 1973 1973 1973 1973 1973 ...
# $ mo                                     : int  10 10 10 10 10 10 10 10 10 10 ...
# $ day                                    : int  1 1 1 1 1 1 1 1 1 1 ...
# $ hour                                   : int  0 1 2 3 4 5 6 7 8 9 ...
# $ Upper_Hafren_avgQ_mm.hr                : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Hafren_avgQ_mm.hr                      : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Tanllwyth_avgQ_mm.hr                   : num  NA NA NA NA NA NA NA NA NA 0.144 ...
# $ Upper_Hore_avgQ_mm.hr                  : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Hore_avgQ_mm.hr                        : num  0.237 0.236 0.232 0.229 0.225 ...
# $ Iago_avgQ_mm.hr                        : num  NA NA NA NA NA ...
# $ Gwy_avgQ_mm.hr                         : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Cyff_avgQ_mm.hr                        : num  NA NA NA NA NA NA NA NA NA NA ...
# $ Cefn_Brwyn_avgQ_mm.hr                  : num  0.242 0.241 0.237 0.232 0.228 ...
# $ Severn_avgQ_mm.hr                      : num  0.235 0.234 0.23 0.226 0.222 ...
# $ P_Upper_Hafren_mm.hr                   : num  NA NA NA NA NA NA NA NA NA NA ...
# $ P_Hafren_Tanllwyth_Hore_UpperHore_mm.hr: num  NA NA NA NA NA NA NA NA NA NA ...
# $ P_Iago_Gwy_mm.hr                       : num  NA NA NA NA NA NA NA NA NA NA ...
# $ P_Cyff_mm.hr                           : num  NA NA NA NA NA NA NA NA NA NA ...
# $ P_Cefn_Brwyn_mm.hr                     : num  NA NA NA NA NA NA NA NA NA NA ...
# $ AllSitesAWS_VPD                        : num  NA NA NA NA NA NA NA NA NA NA ...# Classes 'data.table' and 'data.frame':	326545 obs. of  97 variables:






# first let's estimate a simple RRD, without accounting for nonlinear dependence on precipitation intensity

fileID <- "Robust est demo Hafren"
p <- dat$P_Hafren_Tanllwyth_Hore_UpperHore_mm.hr
q <- dat$Hafren_avgQ_mm.hr             

zz <- ERRA(p=p, q=q, m=100, robust=FALSE)        # first without robust estimation

with(zz, {                                          # write the output to files
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
  fwrite(RRD, paste0(fileID, "_RRD_", options, ".txt"), sep="\t") 
})


zz <- ERRA(p=p, q=q, m=100, robust=TRUE)        # now WITH robust estimation

with(zz, {                                          # write the output to files
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
  fwrite(RRD, paste0(fileID, "_RRD_", options, ".txt"), sep="\t") 
})



# now let's do that again, but this time we'll explicitly account for the nonlinear effects of precipitation intensity

fileID <- "Robust est demo Hafren"
p <- dat$P_Hafren_Tanllwyth_Hore_UpperHore_mm.hr
q <- dat$Hafren_avgQ_mm.hr             

# first non-robust
zz <- ERRA(p=p, q=q, m=100, xknots=c(3, 6, 9, 12), xknot_type="values", robust=FALSE)        # first without robust estimation

with(zz, {                                          # write the output to files
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
  fwrite(wtd_avg_RRD, paste0(fileID, "_avgRRD_", options, ".txt"), sep="\t") 
  fwrite(NRF, paste0(fileID, "_NRF_", options, ".txt"), sep="\t") 
})

# now robust
zz <- ERRA(p=p, q=q, m=100, xknots=c(3, 6, 9, 12), xknot_type="values", robust=TRUE)        # now WITH robust estimation

with(zz, {                                          # write the output to files
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
  fwrite(wtd_avg_RRD, paste0(fileID, "_avgRRD_", options, ".txt"), sep="\t") 
  fwrite(NRF, paste0(fileID, "_NRF_", options, ".txt"), sep="\t") 
})






