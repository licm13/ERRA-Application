

setwd("~/ERRA introduction") # working directory -- change as needed!

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






# For Figure 16, analyze long-tail recession curves with broken-stick methods


fileID <- "Hafren_Fig16"
p <- dat$P_Hafren_Tanllwyth_Hore_UpperHore_mm.hr        # use combined precip
q <- dat$Hafren_avgQ_mm.hr                                   # Hafren streamflow


# Figure 16a,b,c: 1000 lags, non-robust
zz <- ERRA(p=p, q=q, xknots=c(5, 20), m=1000, h=3, xknot_type="even", robust=FALSE)
with(zz, fwrite(wtd_avg_RRD, paste0(fileID, "_avgRRD_", options, ".txt"), sep="\t") )

# Figure 16a,b,c: 1000 lags, robust
zz <- ERRA(p=p, q=q, xknots=c(5, 20), m=1000, h=3, xknot_type="even", robust=TRUE)
with(zz, fwrite(wtd_avg_RRD, paste0(fileID, "_avgRRD_", options, ".txt"), sep="\t") )


# Figure 16a,b,c: broken-stick, non-robust
zz <- ERRA(p=p, q=q, xknots=c(5, 20), m=1000, nk=40, h=3, xknot_type="even", robust=FALSE)
with(zz, fwrite(wtd_avg_RRD, paste0(fileID, "_avgRRD_", options, ".txt"), sep="\t") )

# Figure 16a,b,c: broken-stick, robust
zz <- ERRA(p=p, q=q, xknots=c(5, 20), m=1000, nk=40, h=3, xknot_type="even", robust=TRUE)
with(zz, fwrite(wtd_avg_RRD, paste0(fileID, "_avgRRD_", options, ".txt"), sep="\t") )


# Figure 16d: Hore stream
fileID <- "Hore_Fig16d"
p <- dat$P_Hafren_Tanllwyth_Hore_UpperHore_mm.hr        # use combined precip
q <- dat$Hore_avgQ_mm.hr                                  
zz <- ERRA(p=p, q=q, xknots=c(5, 20), m=1000, nk=40, h=3, xknot_type="even", robust=TRUE)
with(zz, fwrite(wtd_avg_RRD, paste0(fileID, "_avgRRD_", options, ".txt"), sep="\t") )


# Figure 16d: Tanllwyth stream
fileID <- "Tanllwyth_Fig16d"
p <- dat$P_Hafren_Tanllwyth_Hore_UpperHore_mm.hr        # use combined precip
q <- dat$Tanllwyth_avgQ_mm.hr                                
zz <- ERRA(p=p, q=q, xknots=c(5, 20), m=1000, nk=40, h=3, xknot_type="even", robust=TRUE)
with(zz, fwrite(wtd_avg_RRD, paste0(fileID, "_avgRRD_", options, ".txt"), sep="\t") )


# Figure 16d: Severn River
fileID <- "Severn_Fig16d"
p <- dat$P_Hafren_Tanllwyth_Hore_UpperHore_mm.hr        # use combined precip
q <- dat$Severn_avgQ_mm.hr                                 
zz <- ERRA(p=p, q=q, xknots=c(5, 20), m=1000, nk=40, h=3, xknot_type="even", robust=TRUE)
with(zz, fwrite(wtd_avg_RRD, paste0(fileID, "_avgRRD_", options, ".txt"), sep="\t") )


# Figure 16d: Gwy stream
fileID <- "Gwy_Fig16d"
p <- dat$P_Iago_Gwy_mm.hr        # use combined precip
q <- dat$Gwy_avgQ_mm.hr                              
zz <- ERRA(p=p, q=q, xknots=c(5, 20), m=1000, nk=40, h=3, xknot_type="even", robust=TRUE)
with(zz, fwrite(wtd_avg_RRD, paste0(fileID, "_avgRRD_", options, ".txt"), sep="\t") )


# Figure 16d: Cyff stream
fileID <- "Cyff_Fig16d"
p <- dat$P_Cyff_mm.hr       
q <- dat$Cyff_avgQ_mm.hr                         
zz <- ERRA(p=p, q=q, xknots=c(5, 20), m=1000, nk=40, h=3, xknot_type="even", robust=TRUE)
with(zz, fwrite(wtd_avg_RRD, paste0(fileID, "_avgRRD_", options, ".txt"), sep="\t") )


# Figure 16d: Wye River at Cefn Brwyn
fileID <- "Cefn_Brwyn_Fig16d"
p <- dat$P_Cefn_Brwyn_mm.hr        
q <- dat$Cefn_Brwyn_avgQ_mm.hr                 
zz <- ERRA(p=p, q=q, xknots=c(5, 20), m=1000, nk=40, h=3, xknot_type="even", robust=TRUE)
with(zz, fwrite(wtd_avg_RRD, paste0(fileID, "_avgRRD_", options, ".txt"), sep="\t") )





