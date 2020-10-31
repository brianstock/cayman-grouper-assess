# Catch curve analysis of Pickle Bank lengths
# Estimate natural mortality

# Use Chapman-Robson method implemented in FSA / fishR package
# recommended by Smith 2012, with first age group = Peak Plus 1
#   alternatively weighted regression with no truncation of older ages
# convert lengths to ages using overall pars from m2, hierarchical growth model with island-specific k
#     Estimate Std. Error    95% LCI    95% UCI
# Z  0.2762924 0.05440519  0.1696602  0.3829246

# M = 0.276

library(TMB)
library(tidyverse)
library(FSA)

dir <- here()
setwd(dir)

# ----------------------------------------------
# 1. get growth pars from m2 (k varies by island)
# use global/mean pars for Pickle Bank
mods <- readRDS("results/growth_REisland.rds")
meanpars <- summary(mods[[2]]$sdrep)[rownames(summary(mods[[2]]$sdrep))=="theta_global",1]
names(meanpars) <- c("linf","vbk","t0")
CVlen <- summary(mods[[2]]$sdrep)[rownames(summary(mods[[2]]$sdrep))=="CV_L",1]

## ----------------------------------------------------------------------------
## 2. Get length data (all islands, both catch 1987-1995 and lasers 2004-2017)
## ----------------------------------------------------------------------------
dat <- read.csv("data/CaymansNG_DOE_biological_data.csv", header=T)
dat_PB <- dat %>% filter(Island=="PickleBank")
dat_PB$Age <- -1*log(1-dat_PB$Length_cm/meanpars[1])/meanpars[2] + meanpars[3]
amax=30
dat_PB$Age[dat_PB$Length_cm > .99*meanpars[1]] <- amax-0.5
# If fish bigger than linf, assign age = 29.5
# plot(dat_PB$Age, dat_PB$Length_cm,xlim=c(0,30),ylim=c(0,85))

# --------------------------------------------------------
# Use Chapman-Robson method implemented in FSA / fishR package
# recommended by Smith 2012, with first age group = Peak Plus 1
#   alternatively weighted regression with no truncation of older ages
dat_PB$age_bin <- cut(dat_PB$Age, breaks=seq(5.5,amax+0.5,by=1))
levels(dat_PB$age_bin) <- 6:amax
dat_PB$age_bin <- as.numeric(as.character(dat_PB$age_bin))
df <- dat_PB %>% group_by(age_bin) %>% summarize(ct=n())
peak = df$age_bin[df$ct==max(df$ct)]
plot(log(ct)~age_bin,data=df, xlab="Age (yrs)",ylab="Log Catch",pch=19)
thcr <- chapmanRobson(ct~age_bin,data=df,ages2use=(peak+1):amax)
cbind(summary(thcr),confint(thcr))
plot(thcr)
# Hierarchical growth model -- overall pars
#     Estimate Std. Error    95% LCI    95% UCI
# S 75.7692308 2.66244201 70.5509403 80.9875212
# Z  0.2762924 0.05440519  0.1696602  0.3829246

# ----------------------------------------------------
# https://academic.oup.com/icesjms/article/72/1/82/2804320
# Then et al. (2015) 
# Tested M estimators, recommends:
#   1. updated Hoenig estimator
#   2. updated Pauly estimator w/o temp

# 1. updated Hoenig
tmax <- 29
M.Hoenig <- 4.899 * tmax^-0.916
M.Hoenig
# 0.224

# 2. updated Pauly
M.Pauly.new <- 4.118 * meanpars[2]^0.73 * meanpars[1]^-0.333
M.Pauly.new
# 0.245

# old Pauly
temp = 27
M.Pauly <- exp(-0.0066 - 0.279*log(meanpars[1]) + 0.6543*log(meanpars[2]) + 0.4634*log(temp))
M.Pauly
# 0.397
