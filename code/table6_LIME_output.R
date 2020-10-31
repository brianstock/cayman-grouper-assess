# Table 6
# Parameter estimates for 4 models
#   growth
#   recruitment spike 2011
#   depletion

library(LIME)
library(tidyverse)

dir <- here()
setwd(dir)

# note table 6 is in different order than in Fig 9
mods <- list(readRDS("results/LIME_integrated_francis/output/output_CV0.096_sigI0.175_sigC0.2_francis_2.rds"),
            readRDS("results/LIME_integrated_dirichlet/output/output_CV0.096_integrated.rds"),	
            readRDS("results/LIME_fixed_growth_francis/output/output_CV0.096_francis_2.rds"),
            readRDS("results/LIME_fixed_growth_dirichlet/output/output_CV0.096.rds"))
yrs <- mods[[1]]$years

# 2011 multiple of average recruitment with 95% CI
R2011 <- lapply(mods, function(x) matrix(summary(x$Sdreport)[which(rownames(summary(x$Sdreport))=="R2011"),], nrow=1))
sapply(R2011, function(x) round(x[,1],1))
do.call(rbind, lapply(R2011, function(x) round(read_sdreport(x, log=FALSE),1)))

# 5.3 7.8 4.4 5.6

# [1,]  0.6  9.9
# [2,]  3.1 12.5
# [3,]  1.6  7.1
# [4,]  3.6  7.6

# sigma_R
sigR <- lapply(mods, function(x) matrix(summary(x$Sdreport)[which(rownames(summary(x$Sdreport))=="log_sigma_R"),], nrow=1))
sigR.cis <- lapply(sigR, function(x) data.frame(mean=x[,1], low=x[,1]-1.96*x[,2], upp=x[,1]+1.96*x[,2]))
print(exp(do.call(rbind, sigR.cis)), digits=3)
#    mean   low  upp
# 1 0.829 0.568 1.21
# 2 0.932 0.664 1.31
# 3 0.866 0.598 1.25
# 4 0.926 0.658 1.30

# year of SSB min = 2009
yrs[which(res$Report$SB_t==min(res$Report$SB_t))]

# year of #s min = 2008
yrs[which(res$Report$SN_t==min(res$Report$SN_t))]

# % increase SSB 2008-2019 = 414
tail(res$Report$SB_t,1)/min(res$Report$SB_t)*100

# % increase #s 2008-2019 = 341
tail(res$Report$SN_t,1)/min(res$Report$SN_t)*100

# min #: 1776 in 2008
min(res$Report$SN_t)
# 2019 #: 6057
tail(res$Report$SN_t,1)

# Biomass
ssb <- lapply(mods, function(x) matrix(summary(x$Sdreport)[which(rownames(summary(x$Sdreport))=="lSB_t"),], nrow=length(yrs)))
ssb.cis <- lapply(ssb, function(x) data.frame(mean=exp(x[,1]), low=exp(x[,1]-1.96*x[,2]), upp=exp(x[,1]+1.96*x[,2])))

# final year ssbletion = SB/SB0 2019
do.call(rbind, lapply(ssb.cis, function(x) tail(x,1)))
#         mean      low      upp
# 21  18635.35 14208.77 24441.00
# 211 18856.08 15579.38 22821.94
# 212 19619.98 16210.82 23746.10
# 213 18517.33 15472.81 22160.91

# min ssb 2009
do.call(rbind, lapply(ssb.cis, function(x) x[which(x$mean == min(x$mean)), ]))
#         mean      low      upp
# 11  4821.770 3966.285 5861.774
# 10  5793.374 4824.738 6956.478
# 111 4983.578 4132.921 6009.323
# 112 5820.423 4979.744 6803.027

# --------------------------------------------------------
# Depletion (SB/SB0)
dep <- lapply(mods, function(x) matrix(summary(x$Sdreport)[which(rownames(summary(x$Sdreport))=="lD_t"),], nrow=length(yrs)))
dep.cis <- lapply(dep, function(x) data.frame(mean=exp(x[,1]), low=exp(x[,1]-1.96*x[,2]), upp=exp(x[,1]+1.96*x[,2])))

# final year depletion = SB/SB0 2019
# print(do.call(rbind, lapply(dep.cis, function(x) tail(x,1))), digits=3)
print(do.call(rbind, lapply(dep.cis, function(x) tail(x,1))), digits=2)
#      mean   low  upp
# 21  0.899 0.647 1.25
# 211 0.957 0.731 1.25
# 212 1.061 0.802 1.40
# 213 1.045 0.801 1.36

# 21  0.90 0.65 1.2
# 211 0.96 0.73 1.3
# 212 1.06 0.80 1.4
# 213 1.04 0.80 1.4

# min depletion = SB/SB0 2009
print(do.call(rbind, lapply(dep.cis, function(x) x[which(x$mean == min(x$mean)), ])), digits=2)
#     mean  low  upp
# 11  0.23 0.17 0.32
# 10  0.29 0.23 0.38
# 111 0.27 0.20 0.36
# 112 0.33 0.26 0.42

# depletion in 2002
print(do.call(rbind, lapply(dep.cis, function(x) x[4,])), digits=2)
#    mean  low  upp
# 4  0.46 0.39 0.54
# 41 0.55 0.50 0.61
# 42 0.53 0.46 0.61
# 43 0.59 0.55 0.63

# % SSB removed in 2001-2
do.call(rbind, lapply(dep.cis, function(x) 1-x[4,1]))
# [1,] 0.5366119
# [2,] 0.4475096
# [3,] 0.4701871
# [4,] 0.4074884

# --------------------------------------------------------------
# Uncertainty in Linf
LCpars <- lapply(mods, function(x) tail(summary(x$Sdreport)[rownames(summary(x$Sdreport))=="growthpars",1:2],3))
pars.cis <- lapply(LCpars, function(x) data.frame(mean=x[,1], low=x[,1]-1.96*x[,2], upp=x[,1]+1.96*x[,2]))
print(pars.cis, digits=3)
#     mean    low    upp
# 1 81.193 78.056 84.330
# 2  0.141  0.126  0.156
# 3 -0.802 -0.951 -0.654

#     mean    low    upp
# 1 79.343 77.255 81.430
# 2  0.146  0.135  0.158
# 3 -0.778 -0.897 -0.659

S50 <- lapply(mods, function(x) matrix(summary(x$Sdreport)[which(rownames(summary(x$Sdreport))=="S50_f"),], nrow=1))
S50.cis <- lapply(S50, function(x) data.frame(mean=x[,1], low=x[,1]-1.96*x[,2], upp=x[,1]+1.96*x[,2]))
print(do.call(rbind, S50.cis), digits=3)
# 1 61.8 59.2 64.4
# 2 63.8 62.4 65.3
# 3 61.9 59.4 64.4
# 4 64.0 62.7 65.3

S95 <- lapply(mods, function(x) matrix(summary(x$Sdreport)[which(rownames(summary(x$Sdreport))=="S95_f"),], nrow=1))
S95.cis <- lapply(S95, function(x) data.frame(mean=x[,1], low=x[,1]-1.96*x[,2], upp=x[,1]+1.96*x[,2]))
print(do.call(rbind, S95.cis), digits=3)
# 1 66.0 59.4 72.6
# 2 69.8 67.3 72.3
# 3 66.1 59.9 72.4
# 4 70.1 67.8 72.3
