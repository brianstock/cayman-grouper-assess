# Collect results from best LIME model (integrated, francis weighting)
#   Figure 8: 6 panel LIME output 
#   Figure S6: fit to length composition data
#   Figure S8b: likelihood profile Linf
#   Figure S7: residuals index and mean length

# Parameter estimates (in text)
#   growth
#   recruitment spike 2011
#   SSB
#   numbers
#   depletion

dir <- here()
setwd(dir)
plots_dir <- file.path(dir, "plots")
source("code/plot_LIME_output.R")
source("code/plot_stdresiduals.R")
source("code/plot_LCfits.R")

library(LIME)
library(ggsidekick)
library(tidyverse)
library(TMB)

# read in fit from LIME integrated francis weighting (fig8a_run_LIME_integrated_francis.R)
#   CV_L = 0.96, sigI = 0.175, sigC = 0.2
res <- readRDS("results/LIME_integrated_francis/output/output_CV0.096_sigI0.175_sigC0.2_francis_2.rds")
yrs <- res$years

# selectivity
S50 <- matrix(summary(res$Sdreport)[which(rownames(summary(res$Sdreport))=="S50_f"),], nrow=1)
read_sdreport(S50, log=FALSE)
res$Report$S50_f
# 59.24563 64.42461
# 61.83512

S95 <- matrix(summary(res$Sdreport)[which(rownames(summary(res$Sdreport))=="S95_f"),], nrow=1)
read_sdreport(S95, log=FALSE)
res$Report$S95_f
# [1] 59.40353 72.60420
# [1] 66.00386

# sigmaR 
logsigR <- matrix(summary(res$Sdreport)[which(rownames(summary(res$Sdreport))=="log_sigma_R"),], nrow=1)
read_sdreport(logsigR, log=TRUE)
exp(res$Report$log_sigma_R)
# [1] 0.5683289 1.2102877
# [1] 0.8293621

# average recruitment (R0)
R0 <- matrix(summary(res$Sdreport)[which(rownames(summary(res$Sdreport))=="R0"),], nrow=1)
read_sdreport(R0, log=FALSE)
R0[1,1]
# [1] 5321.410 8473.918
# [1] 6897.664

# 2011 multiple of average recruitment with 95% CI
R2011 <- matrix(summary(res$Sdreport)[which(rownames(summary(res$Sdreport))=="R2011"),], nrow=1)
round(R2011[,1],1)
round(read_sdreport(R2011, log=FALSE),1)
# [1] 5.3
# [1] 0.6 9.9

# year of SSB min = 2009
yrs[which(res$Report$SB_t==min(res$Report$SB_t))]

# year of #s min = 2008
yrs[which(res$Report$SN_t==min(res$Report$SN_t))]

# % increase SSB 2008-2019 = 395
tail(res$Report$SB_t,1)/min(res$Report$SB_t)*100

# % increase #s 2008-2019 = 334
tail(res$Report$SN_t,1)/min(res$Report$SN_t)*100

# min #: 1757 in 2008
min(res$Report$SN_t)
# 2019 #: 5873
tail(res$Report$SN_t,1)

est_rdev <- est_F <- rep(1,length(yrs))
est_rdev[which(is.na(res$Inputs$Map$Nu_input))] = 0
est_F[which(is.na(res$Inputs$Map$log_F_ft))] = 0

# Figure 8: 6-panel output from best LIME integrated model w/ Francis weighting
cairo_pdf(filename=file.path(plots_dir,"fig8_LIME_best.pdf"), width=10, height=5.5)
plot_output(Inputs=res$Inputs, Report=res$Report, Sdreport=res$Sdreport, true_years=yrs,
	plot=c("ML","Rec","Ind","SB","Fish","Selex"), set_ylim=list("Ind" = c(0,8500), "Rec" = c(0,40000)), relative=FALSE, est=list(Rec=est_rdev, Fish=est_F))
dev.off()

# Figure S6: fit to length composition data
png(file.path(plots_dir,"figS7_LIME_LCfit.png"),width=10,height=7,units="in",res=300)
g <- plot_LCfits(Inputs=res$Inputs, Report=res$Report, n=TRUE, year_labels=yrs)
print(g + theme_sleek() + theme(legend.position='none'))
dev.off()

# --------------------------------------------------------------
# Uncertainty in Linf
LCpars <- tail(summary(res$Sdreport)[rownames(summary(res$Sdreport))=="growthpars",1:2],3)
rownames(LCpars) <- c("linf","vbk","t0")
#        Estimate  Std. Error
# linf 81.0894061 1.612818671
# vbk   0.1413890 0.007777532
# t0   -0.8012961 0.075443251

# Little Cayman specific growth pars when fit outside assessment
readRDS("results/growth_REisland_LCpars.rds")
#        linf         vbk          t0       CVlen 
# 80.23134337  0.14025094 -0.83162114  0.09233058 

# Overall/mean growth pars
growthmu <- summary(res$Sdreport)[rownames(summary(res$Sdreport))=="growth_mu",1:2]
growthmu[2,1] <- exp(growthmu[2,1])
rownames(growthmu) <- c("linf","vbk","t0")
#        Estimate Std. Error
# linf 81.0894061 1.61281867
# vbk   0.1528427 0.06112610
# t0   -0.8012961 0.07544325

# Overall growth pars when fit outside assessment
readRDS("results/growth_REisland_globalpars.rds")
#        linf         vbk          t0       CVlen 
# 80.23134337  0.15450038 -0.83162114  0.09233058 

# ----------------------------------------------------
# Figure S7: standardized residuals for index and mean length
stdres <- plot_stdresiduals(res)
png(file.path(plots_dir,paste0("figS7_residuals_LIME_integrated_francis.png")),width=6,height=3.5,units="in",res=300)
print(stdres[[1]])
dev.off()

# --------------------------------------------------------
# Figure S8b: likelihood profile over Linf for best model
dyn.load("code/LIME_integrated_REisland.so")
prof <- tmbprofile(res$obj, which(names(res$obj$par)=="growth_mu")[1])
saveRDS(prof, "results/LIME_integrated_francis/linf_profile.rds")
png("plots/profile_Linf_integrated_francis.png", width=7, height=7, units='in', res=300)
plot(prof[,1], prof[,2], type='l', xlab="Linf", ylab="Negative log-likelihood",bty='l')
i.min <- which(prof[,2] == min(prof[,2], na.rm=T))
points(prof[i.min,1], prof[i.min,2], bg='black',pch=21)
abline(h=(prof[i.min,2]+1.92), lty=3) # 1.92 NLL units up = 95% CI
dev.off()

