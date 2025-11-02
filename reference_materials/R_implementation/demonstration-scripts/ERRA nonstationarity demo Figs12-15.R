

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








fileID <- "Hafren_Fig12a"
p <- dat$P_Hafren_Tanllwyth_Hore_UpperHore_mm.hr
q <- dat$Hafren_avgQ_mm.hr             


# construct splitting parameter list based on 1-hour lagged discharge
breakpts <- c(0.1, 0.2, 0.5, 1)
antQparams <- list(crit = list(q) ,                  # must be LIST of criterion variables; each must be of same length as P
                   crit_label = c("antQ") ,                  # must be vector of same length as crit, with string labels for each criterion
                   crit_lag = 1 ,                       # integer vector of same length as crit;
                   breakpts = list(breakpts) ,         # LIST of vectors of breakpoints (each of length >=1) for each criterion
                   pct_breakpts = c(FALSE) ,             # boolean vector of same length as crit: FALSE means breakpoints are raw values, TRUE means breakpoints are percentiles
                   thresh = 0.0 ,                       # numeric vector of same length as crit
                   by_bin = c(FALSE)                    # boolean vector of same length as crit
) # end split_params


# For figure 12a, split by antQ but don't account for differences in precipitation intensity
# (note this over-estimates the effects of antQ)
zz <- ERRA(p=p, q=q, m=100, split_params=antQparams)

with(zz, {                                          # write the output to files
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
  fwrite(RRD, paste0(fileID, "_RRD_", options, ".txt"), sep="\t") 
})


fileID <- "Hafren_Fig12bcd_13"

# For figure 12b,c,d and figure 13, split by antQ *and* take nonlinear effects of precipitation intensity into account
zz <- ERRA(p=p, q=q, xknots=c(3, 6, 9, 12), xknot_type="values", m=100, split_params=antQparams)

with(zz, {                                          # write the output to files
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
  fwrite(wtd_avg_RRD, paste0(fileID, "_avgRRD_", options, ".txt"), sep="\t") 
  fwrite(NRF, paste0(fileID, "_NRF_", options, ".txt"), sep="\t") 
})






# For figure 14 (a cautionary tale), split by VPD alone

fileID <- "Hafren_Fig14"

p <- dat$P_Hafren_Tanllwyth_Hore_UpperHore_mm.hr
q <- dat$Hafren_avgQ_mm.hr             
vpd <- dat$AllSitesAWS_VPD

breakpts <- c(0.0315, 0.0576, 0.0931, 0.185) # quintiles of VPD distribution
antVPDparams <- list(crit = list(vpd) ,                  # must be LIST of criterion variables; each must be of same length as P
                     crit_label = c("antVPD") ,                  # must be vector of same length as crit, with string labels for each criterion
                     crit_lag = 1 ,                       # integer vector of same length as crit;
                     breakpts = list(breakpts) ,         # LIST of vectors of breakpoints (each of length >=1) for each criterion
                     pct_breakpts = c(FALSE) ,             # boolean vector of same length as crit: FALSE means breakpoints are raw values, TRUE means breakpoints are percentiles
                     thresh = 0.0 ,                       # numeric vector of same length as crit
                     by_bin = c(FALSE)                    # boolean vector of same length as crit
) # end antVPDparams


zz <- ERRA(p=p, q=q, xknots=c(2, 4, 6, 10), xknot_type="values", m=100, split_params=antVPDparams)

with(zz, {                                          # write the output to files
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
  fwrite(NRF, paste0(fileID, "_NRF_", options, ".txt"), sep="\t") 
})




# For figure 15, split by antecedent VPD and antecedent Q

fileID <- "Hafren_Fig15"
p <- dat$P_Hafren_Tanllwyth_Hore_UpperHore_mm.hr
q <- dat$Hafren_avgQ_mm.hr
vpd <- dat$AllSitesAWS_VPD


qbreakpts <- c(0.1, 0.2, 0.5, 1)
vbreakpts <- c(0.0932) #60% of VPD distribution
antVPD_Qparams <- list(crit = list(vpd, q) ,                  # must be LIST of criterion variables; each must be of same length as P
                   crit_label = c("antVPD", "antQ") ,                  # must be vector of same length as crit, with string labels for each criterion
                   crit_lag = c(1, 1) ,                       # integer vector of same length as crit;
                   breakpts = list(vbreakpts, qbreakpts) ,         # LIST of vectors of breakpoints (each of length >=1) for each criterion
                   pct_breakpts = c(FALSE, FALSE) ,             # boolean vector of same length as crit: FALSE means breakpoints are raw values, TRUE means breakpoints are percentiles
                   thresh = c(-999, 0.0) ,                       # numeric vector of same length as crit
                   by_bin = c(FALSE, FALSE)                    # boolean vector of same length as crit
) # end antVPD_Qparams

zz <- ERRA(p=p, q=q, xknots=c(2, 4, 6, 10), xknot_type="values", m=100, split_params=antVPD_Qparams)

with(zz, {                                          # write the output to files
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
  fwrite(NRF, paste0(fileID, "_NRF_", options, ".txt"), sep="\t") 
})


























