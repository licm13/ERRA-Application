

# this script tests out ERRA on MOPEX SF New River hourly records


setwd("~/ERRA introduction")  #

rm(list=ls())  #clear environment

library(data.table)
library(dplyr)


source("ERRA_v1.0.R")  # ensemble rainfall runoff analysis script






dat <- fread("MOPEX SacoR hourly.txt", header=TRUE, sep="\t", na.strings=c("NA",".","","#N/A"))    #this is the test data file


# $ p        : num  0.07 1.54 0.04 0.04 0.66 1.02 0.29 1.47 0 0 ...
# $ PET      : num  0.0578 0.0578 0.0578 0.0578 0.0578 ...
# $ q        : num  0.0397 0.0397 0.0397 0.0397 0.0397 ...
# $ yr       : int  1987 1987 1987 1987 1987 1987 1987 1987 1987 1987 ...
# $ mm       : int  10 10 10 10 10 10 10 10 10 10 ...
# $ dd       : int  21 21 21 21 21 21 21 21 21 21 ...
# $ PET_mm/d : num  1.39 1.39 1.39 1.39 1.39 ...
# $ Tmax     : num  10.9 10.9 10.9 10.9 10.9 ...
# $ Tmin     : num  1.44 1.44 1.44 1.44 1.44 ...





# first, a simple RRD for the snow-free months
fileID <- "SacoR_snowfree"

p <- dat$p        
q <- dat$q

snowfree <- ifelse((dat$mm>4)&(dat$mm<11), 1, 0)    # define snowfree months

zz <- ERRA(p=p, q=q, m=100, Qfilter=snowfree)       # run ERRA, excluding snowy months

with(zz, {                                          # write the output to files
  fwrite(RRD, paste0(fileID, "_RRD_", options, ".txt"), sep="\t") 
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
})





# now, a nonlinear analysis with xknot_type="even" to get a first impression of what we have

fileID <- "SacoR_snowfree"

p <- dat$p        
q <- dat$q

snowfree <- ifelse((dat$mm>4)&(dat$mm<11), 1, 0)    # define snowfree months

# run ERRA with xknot_type="even"
zz <- ERRA(p=p, q=q, xknots=c(6,40), xknot_type="even", m=100, Qfilter=snowfree, show_top_xknot=TRUE)

zz$knot_peakstats$knot   # inspect the resulting outputs (in particular the knot values that were found)



fileID <- "SacoR_snowfree Figs7-9"

# now re-run ERRA with fixed xknot values close to those found above 
# (this gives us the numerical values behind Figures 7, 8, and 9 of K2024)

zz <- ERRA(p=p, q=q, xknots=c(1.5, 3, 4.5, 6, 7.5, 9), xknot_type="values", m=100, Qfilter=snowfree, show_top_xknot=TRUE)

with(zz, {                                          # write the output to files
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
  fwrite(NRF, paste0(fileID, "_NRF_", options, ".txt"), sep="\t") 
})


fileID <- "SacoR_snowfree Figs8-9"

# now aggregate time steps to 3 hours (for Figures 8 and 9 of K2024)
# note we need different xknots because aggregation shortens the tail of the precipitation intensity distribution

zz <- ERRA(p=p, q=q, xknots=c(1, 2, 3, 4, 5, 6), xknot_type="values", agg=3, m=34, Qfilter=snowfree, show_top_xknot=TRUE)

with(zz, {                                          # write the output to files
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
  fwrite(NRF, paste0(fileID, "_NRF_", options, ".txt"), sep="\t") 
})




# now aggregate time steps to 6 hours (for Figures 8 and 9 of K2024)
# note we need different xknots because aggregation shortens the tail of the precipitation intensity distribution

zz <- ERRA(p=p, q=q, xknots=c(0.7, 1.4, 2.1, 2.8, 3.5, 4.2), xknot_type="values", agg=6, m=17, Qfilter=snowfree, show_top_xknot=TRUE)

with(zz, {                                          # write the output to files
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
  fwrite(NRF, paste0(fileID, "_NRF_", options, ".txt"), sep="\t") 
})





# now aggregate time steps to 12 hours (for Figures 8 and 9 of K2024)
# note we need different xknots because aggregation shortens the tail of the precipitation intensity distribution

zz <- ERRA(p=p, q=q, xknots=c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0), xknot_type="values", agg=12, m=9, Qfilter=snowfree, show_top_xknot=TRUE)

with(zz, {                                          # write the output to files
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
  fwrite(NRF, paste0(fileID, "_NRF_", options, ".txt"), sep="\t") 
})




