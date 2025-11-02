

setwd("~/ERRA introduction")   # working directory -- change as needed!

rm(list=ls())  #clear environment

library(data.table)
library(dplyr)



source("ERRA_v1.0.R")  # ensemble rainfall runoff analysis script





dat <- fread("MOPEX NantahalaR hourly.txt", header=TRUE, sep="\t", na.strings=c("NA",".","","#N/A"))    #this is the data file


# $ date_time: chr  "1985-10-08 23:00:00" "1985-10-09 00:00:00" "1985-10-09 01:00:00" "1985-10-09 02:00:00" ...
# $ p        : num  0 0 0 0 0 0 0 0 0 0 ...
# $ q        : num  0.0563 0.0563 0.0563 0.0563 0.0563 ...
# $ PET      : num  0.0719 0.0706 0.0706 0.0706 0.0706 ...
# $ yr       : int  1985 1985 1985 1985 1985 1985 1985 1985 1985 1985 ...
# $ mm       : int  10 10 10 10 10 10 10 10 10 10 ...
# $ dd       : int  8 9 9 9 9 9 9 9 9 9 ...
# $ PET_mm/d : num  1.73 1.69 1.69 1.69 1.69 ...
# $ year     : num  1986 1986 1986 1986 1986 ...





fileID <- "Baseline filtering demo Nantahala"
p <- dat$p
q <- dat$q             


# construct splitting parameter list based on 1-hour lagged discharge
breakpts <- c(0.1, 0.2, 0.4)
antQparams <- list(crit = list(q) ,                  # must be LIST of criterion variables; each must be of same length as P
                   crit_label = c("antQ") ,                  # must be vector of same length as crit, with string labels for each criterion
                   crit_lag = 1 ,                       # integer vector of same length as crit;
                   breakpts = list(breakpts) ,         # LIST of vectors of breakpoints (each of length >=1) for each criterion
                   pct_breakpts = c(FALSE) ,             # boolean vector of same length as crit: FALSE means breakpoints are raw values, TRUE means breakpoints are percentiles
                   thresh = 0.0 ,                       # numeric vector of same length as crit
                   by_bin = c(FALSE)                    # boolean vector of same length as crit
) # end split_params


# without baseline filtering
zz <- ERRA(p=p, q=q, xknots=c(4, 20), xknot_type="even", m=200, h=3, split_params=antQparams)

with(zz, {                                          # write the output to files
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
  fwrite(wtd_avg_RRD, paste0(fileID, "_avgRRD_", options, ".txt"), sep="\t") 
  fwrite(NRF, paste0(fileID, "_NRF_", options, ".txt"), sep="\t") 
  fwrite(Qcomp, paste0(fileID, "_Qcomp_", options, ".txt"), sep="\t") 
})


with(zz$Qcomp, {
  NSE <- 1-var(Qresidual, na.rm=TRUE)/var(Q[!is.na(Qresidual)])
  print(paste0("NSE ", NSE))
})


# with baseline filtering
zz <- ERRA(p=p, q=q, xknots=c(4, 20), xknot_type="even", m=200, h=3, fq=0.1, split_params=antQparams)


with(zz, {                                          # write the output to files
  fwrite(peakstats, paste0(fileID, "_peakstats_", options, ".txt"), sep="\t") 
  fwrite(wtd_avg_RRD, paste0(fileID, "_avgRRD_", options, ".txt"), sep="\t") 
  fwrite(NRF, paste0(fileID, "_NRF_", options, ".txt"), sep="\t") 
  fwrite(Qcomp, paste0(fileID, "_Qcomp_", options, ".txt"), sep="\t") 
})


with(zz$Qcomp, {
  NSE <- 1-var(Qresidual, na.rm=TRUE)/var(Q[!is.na(Qresidual)])
  print(paste0("NSE ", NSE))
})








