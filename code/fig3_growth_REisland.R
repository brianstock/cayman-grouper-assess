# Little Cayman Nassau Grouper length assessment
# Fit growth outside assessment
#   by island, sharing parameters as in Helser and Lai (2004) p.401-2
#   uses growth_REisland.cpp
# Saves intermediate results
#   growth_REisland.rds
# Makes Figure 3

library(tidyverse)
library(TMB)

dir <- here()
setwd(dir)
plots_dir <- file.path(dir, "plots")
res_dir <- file.path(dir, "results")
compile("code/growth_REisland.cpp","-O0 -g")
dyn.load("code/growth_REisland.so")

# ----------------------------------------------------------------------
# Get data
# convert lengths in inches to cm
dat <- read.csv("data/NG_pre2001_biological_data_2002.csv", header=T)
dat$Length_cm <- dat$Length
dat$Length_cm[dat$Length_units == "in"] <- dat$Length[dat$Length_units == "in"] * 2.54

x <- read.csv("data/NG_LC_lengths_2003_2019.csv")
x <- x %>% filter(usable!="n")
# x %>% group_by(year) %>% summarize(n=n(),median=median(length),mean=mean(length),range=max(length)-min(length),low=quantile(length,.25),high=quantile(length,.75))

dat.colnames <- names(dat)
dat_LC_lasers <- data.frame(matrix(ncol = length(dat.colnames), nrow = dim(x)[1]))
colnames(dat_LC_lasers) <- dat.colnames
dat_LC_lasers$Length_cm = x$length/10 # conver mm to cm
dat_LC_lasers$Island = "Little" # all from Little Cayman
dat_LC_lasers$Year = x$year
dat_LC_lasers$Date = x$date

# add new factor, whether data was from lasers or catch
dat$type <- "catch"
dat_LC_lasers$type <- "lasers"

# combine the two data frames
dat[] <- lapply(dat, as.character) # convert all columns to character
dat_LC_lasers[] <- lapply(dat_LC_lasers, as.character)
dat_all <- rbind(dat, dat_LC_lasers)

# reconvert columns to appropriate data classes
dat_all$Length_cm <- as.numeric(dat_all$Length_cm)
dat_all$Year <- factor(dat_all$Year, levels=1987:2019)
dat_all$type <- factor(dat_all$type)
dat_all$Date <- as.Date(dat_all$Date)
dat_all$Age <- as.numeric(dat_all$Age)
dat_all$Weight <- as.numeric(dat_all$Weight)
dat_all$Weight_kg <- ifelse(dat_all$Weight_units=="kg", dat_all$Weight, dat_all$Weight*0.453592) # convert lbs to kg

# note: vonB parameters sensitive to age for the few young fish (0-2 yr-old non-mature fish caught off SPAG)
#   best assigned age for the 0-2 yo = ring count + 0.5, as they were caught in July and Aug, 6 months after spawning
#   all others were caught on SPAG in Jan-Feb, so add 1 year to ring count (see 1988 Table 1 note)
dat_all$Age <- ifelse(dat_all$Age > 2, dat_all$Age + 1, dat_all$Age + 0.5)

# add 2012 YOY (assumed 1 year old, used in growth curve but not length analysis)
x <- c(8, 16.4, 15, 15, 23.25, 16.8, 14.6, 19.9, 12.3, 15.1, 14.7, 15.7, 18.6, 17.4, 22, 23.1, 22, 20, 19.1, 15, 19.9, 20.5, 20.1, 15, 15)
dat_yoy <- data.frame(matrix(NA,ncol = length(colnames(dat_all)), nrow = length(x)))
colnames(dat_yoy) <- colnames(dat_all)
dat_yoy$Length_cm = x
dat_yoy$Year = 2012
dat_yoy$Age = 1
dat_yoy$Island = "Little"
dat_yoy$type <- "catch"
dat_all <- rbind(dat_all, dat_yoy)
dat_all$Island[dat_all$Island=="Grand_12MileBanks"] = "Grand"
dat_all$Year[dat_all$Island=="PickleBank"] <- 2000

dat_growth <- filter(dat_all, !is.na(dat_all$Age) & !is.na(dat_all$Length_cm)) # for growth curve only, leave in the fish < 2 years old
dat_growth_fit <- select(dat_growth, Age, Length_cm) %>% as.matrix
island <- sapply(dat_growth$Island, match, c("Brac", "Grand","Little"))
n.islands <- length(unique(island))

dat_catch <- filter(dat_all, type=="catch")
dat_catch <- dat_catch %>% select(Year, Island, Age, Length_cm, Weight_kg)

# -----------------------------------------------------------
# set up TMB model
Data <- list(dat_growth=dat_growth_fit, n_islands=n.islands, island=island)

# initial values for mean growth pars
vbk = 0.171
linf = 78.2
t0 = -0.608
CVlen = 0.098
Parameters <- list(mu=c(linf, log(vbk), t0),
				log_sigma_mu=c(rep(log(0.1),n.islands)),
				rho_untrans=rep(0,n.islands),
				log_CV_L=log(CVlen),
				emat=matrix(0,nrow=3,ncol=n.islands))
Map <- list()
mods <- list()

# ----------------------------------------------------
# m1: ---
tmp <- Parameters$emat
tmp[] <- NA
ind.notNA <- which(!is.na(tmp))
tmp[ind.notNA] <- 1:length(ind.notNA)
Map$emat = factor(tmp)
tmp <- Parameters$log_sigma_mu
tmp[] <- NA
ind.notNA <- which(!is.na(tmp))
tmp[ind.notNA] <- 1:length(ind.notNA)
Map$log_sigma_mu = factor(tmp)
tmp <- Parameters$rho_untrans
tmp[] <- NA
ind.notNA <- which(!is.na(tmp))
tmp[ind.notNA] <- 1:length(ind.notNA)
Map$rho_untrans = factor(tmp)
Random <- c()

tmblist <- list(Parameters=Parameters, Data=Data, Random=Random, Map=Map)
mods[[1]] <- TMB::MakeADFun(data=tmblist$Data, parameters=tmblist$Parameters, map=tmblist$Map, DLL="growth_REisland")
mods[[1]]$opt <- TMBhelper::fit_tmb(obj=mods[[1]], newtonsteps=3, getsd=FALSE)
mods[[1]]$rep <- mods[[1]]$report()
mods[[1]]$sdrep <- sdreport(mods[[1]])
summary(mods[[1]]$sdrep)

# -----------------------------------------------------
# m2: k
Map <- list()
tmp <- Parameters$emat
tmp[c(1,3),] <- NA
ind.notNA <- which(!is.na(tmp))
tmp[ind.notNA] <- 1:length(ind.notNA)
Map$emat = factor(tmp)
tmp <- Parameters$log_sigma_mu
tmp[c(1,3)] <- NA
ind.notNA <- which(!is.na(tmp))
tmp[ind.notNA] <- 1:length(ind.notNA)
Map$log_sigma_mu = factor(tmp)
tmp <- Parameters$rho_untrans
tmp[] <- NA
ind.notNA <- which(!is.na(tmp))
tmp[ind.notNA] <- 1:length(ind.notNA)
Map$rho_untrans = factor(tmp)
Random <- c("emat")

tmblist <- list(Parameters=Parameters, Data=Data, Random=Random, Map=Map)
mods[[2]] <- TMB::MakeADFun(data=tmblist$Data, parameters=tmblist$Parameters, random=tmblist$Random, map=tmblist$Map, DLL="growth_REisland")
mods[[2]]$opt <- TMBhelper::fit_tmb(obj=mods[[2]], newtonsteps=3, getsd=FALSE)
mods[[2]]$rep <- mods[[2]]$report()
mods[[2]]$sdrep <- sdreport(mods[[2]])
summary(mods[[2]]$sdrep)

# ------------------------------------------------------
# m3: linf
Map <- list()
tmp <- Parameters$emat
tmp[c(2,3),] <- NA
ind.notNA <- which(!is.na(tmp))
tmp[ind.notNA] <- 1:length(ind.notNA)
Map$emat = factor(tmp)
tmp <- Parameters$log_sigma_mu
tmp[c(2,3)] <- NA
ind.notNA <- which(!is.na(tmp))
tmp[ind.notNA] <- 1:length(ind.notNA)
Map$log_sigma_mu = factor(tmp)
tmp <- Parameters$rho_untrans
tmp[] <- NA
ind.notNA <- which(!is.na(tmp))
tmp[ind.notNA] <- 1:length(ind.notNA)
Map$rho_untrans = factor(tmp)
Random <- c("emat")

tmblist <- list(Parameters=Parameters, Data=Data, Random=Random, Map=Map)
mods[[3]] <- TMB::MakeADFun(data=tmblist$Data, parameters=tmblist$Parameters, random=tmblist$Random, map=tmblist$Map, DLL="growth_REisland")
mods[[3]]$opt <- TMBhelper::fit_tmb(obj=mods[[3]], newtonsteps=3, getsd=FALSE)
mods[[3]]$rep <- mods[[3]]$report()
mods[[3]]$sdrep <- sdreport(mods[[3]])
summary(mods[[3]]$sdrep)

# ------------------------------------------------------
# m4: t0
Map <- list()
tmp <- Parameters$emat
tmp[c(1,2),] <- NA
ind.notNA <- which(!is.na(tmp))
tmp[ind.notNA] <- 1:length(ind.notNA)
Map$emat = factor(tmp)
tmp <- Parameters$log_sigma_mu
tmp[c(1,2)] <- NA
ind.notNA <- which(!is.na(tmp))
tmp[ind.notNA] <- 1:length(ind.notNA)
Map$log_sigma_mu = factor(tmp)
tmp <- Parameters$rho_untrans
tmp[] <- NA
ind.notNA <- which(!is.na(tmp))
tmp[ind.notNA] <- 1:length(ind.notNA)
Map$rho_untrans = factor(tmp)
Random <- c("emat")

tmblist <- list(Parameters=Parameters, Data=Data, Random=Random, Map=Map)
mods[[4]] <- TMB::MakeADFun(data=tmblist$Data, parameters=tmblist$Parameters, random=tmblist$Random, map=tmblist$Map, DLL="growth_REisland")
mods[[4]]$opt <- TMBhelper::fit_tmb(obj=mods[[4]], newtonsteps=3, getsd=FALSE)
mods[[4]]$rep <- mods[[4]]$report()
mods[[4]]$sdrep <- sdreport(mods[[4]])
summary(mods[[4]]$sdrep)

# -----------------------------------------------------
# m5: linf + k
Map <- list()
tmp <- Parameters$emat
tmp[3,] <- NA
ind.notNA <- which(!is.na(tmp))
tmp[ind.notNA] <- 1:length(ind.notNA)
Map$emat = factor(tmp)
tmp <- Parameters$log_sigma_mu
tmp[3] <- NA
ind.notNA <- which(!is.na(tmp))
tmp[ind.notNA] <- 1:length(ind.notNA)
Map$log_sigma_mu = factor(tmp)
tmp <- Parameters$rho_untrans
tmp[2:3] <- NA
ind.notNA <- which(!is.na(tmp))
tmp[ind.notNA] <- 1:length(ind.notNA)
Map$rho_untrans = factor(tmp)
Random <- c("emat")

tmblist <- list(Parameters=Parameters, Data=Data, Random=Random, Map=Map)
mods[[5]] <- TMB::MakeADFun(data=tmblist$Data, parameters=tmblist$Parameters, random=tmblist$Random, map=tmblist$Map, DLL="growth_REisland")
mods[[5]]$opt <- TMBhelper::fit_tmb(obj=mods[[5]], newtonsteps=3, getsd=FALSE)
mods[[5]]$rep <- mods[[5]]$report()
mods[[5]]$sdrep <- sdreport(mods[[5]])
summary(mods[[5]]$sdrep)

# -----------------------------------------------------
# m6: linf + t0
Map <- list()
tmp <- Parameters$emat
tmp[2,] <- NA
ind.notNA <- which(!is.na(tmp))
tmp[ind.notNA] <- 1:length(ind.notNA)
Map$emat = factor(tmp)
tmp <- Parameters$log_sigma_mu
tmp[2] <- NA
ind.notNA <- which(!is.na(tmp))
tmp[ind.notNA] <- 1:length(ind.notNA)
Map$log_sigma_mu = factor(tmp)
tmp <- Parameters$rho_untrans
tmp[c(1,3)] <- NA
ind.notNA <- which(!is.na(tmp))
tmp[ind.notNA] <- 1:length(ind.notNA)
Map$rho_untrans = factor(tmp)
Random <- c("emat")

tmblist <- list(Parameters=Parameters, Data=Data, Random=Random, Map=Map)
mods[[6]] <- TMB::MakeADFun(data=tmblist$Data, parameters=tmblist$Parameters, random=tmblist$Random, map=tmblist$Map, DLL="growth_REisland")
mods[[6]]$opt <- TMBhelper::fit_tmb(obj=mods[[6]], newtonsteps=3, getsd=FALSE)
mods[[6]]$rep <- mods[[6]]$report()
mods[[6]]$sdrep <- sdreport(mods[[6]])
summary(mods[[6]]$sdrep)

# -----------------------------------------------------
# m7: k + t0
Map <- list()
tmp <- Parameters$emat
tmp[1,] <- NA
ind.notNA <- which(!is.na(tmp))
tmp[ind.notNA] <- 1:length(ind.notNA)
Map$emat = factor(tmp)
tmp <- Parameters$log_sigma_mu
tmp[1] <- NA
ind.notNA <- which(!is.na(tmp))
tmp[ind.notNA] <- 1:length(ind.notNA)
Map$log_sigma_mu = factor(tmp)
tmp <- Parameters$rho_untrans
tmp[1:2] <- NA
ind.notNA <- which(!is.na(tmp))
tmp[ind.notNA] <- 1:length(ind.notNA)
Map$rho_untrans = factor(tmp)
Random <- c("emat")

tmblist <- list(Parameters=Parameters, Data=Data, Random=Random, Map=Map)
mods[[7]] <- TMB::MakeADFun(data=tmblist$Data, parameters=tmblist$Parameters, random=tmblist$Random, map=tmblist$Map, DLL="growth_REisland")
mods[[7]]$opt <- TMBhelper::fit_tmb(obj=mods[[7]], newtonsteps=3, getsd=FALSE)
mods[[7]]$rep <- mods[[7]]$report()
mods[[7]]$sdrep <- sdreport(mods[[7]])
summary(mods[[7]]$sdrep)

# -----------------------------------------------------
# m8: linf, k, t0 (all)
Map <- list()
Random <- c("emat")

tmblist <- list(Parameters=Parameters, Data=Data, Random=Random, Map=Map)
mods[[8]] <- TMB::MakeADFun(data=tmblist$Data, parameters=tmblist$Parameters, random=tmblist$Random, map=tmblist$Map, DLL="growth_REisland")
mods[[8]]$opt <- TMBhelper::fit_tmb(obj=mods[[8]], newtonsteps=3, getsd=FALSE)
mods[[8]]$rep <- mods[[8]]$report()
mods[[8]]$sdrep <- sdreport(mods[[8]])
summary(mods[[8]]$sdrep)

# -------------------------------------------------------
# Tables 1 and 2
saveRDS(mods, file.path(res_dir,"growth_REisland.rds"))
if(is.null(names(mods))) names(mods) <- paste0("m",1:length(mods))

# num fixed effect parameters
npar <- sapply(mods, function(x) length(x$par))
npar

# converged (pos def Hessian?)
pdHess <- sapply(mods[1:7], function(x) x$sdrep$pdHess)
 #   m1    m2    m3    m4    m5    m6    m7 
 # TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE

aic <- round(sapply(mods[1:7], function(x) 2*(x$opt$obj + length(x$opt$par))), 1)
#     m1     m2     m3     m4     m5     m6     m7 
# 3112.8 3078.4 3082.2 3087.9 3090.7 3086.4 3088.7 

min.aic <- min(aic)
daic <- aic-min.aic
#   m1   m2   m3   m4   m5   m6   m7 
# 34.4  0.0  3.8  9.5 12.3  8.0 10.3

# Little Cayman specific growth pars from m2
LCpars <- tail(summary(mods[[2]]$sdrep)[rownames(summary(mods[[2]]$sdrep))=="theta",],3)
LCpars <- rbind(LCpars, summary(mods[[2]]$sdrep)[rownames(summary(mods[[2]]$sdrep))=="CV_L",])
rownames(LCpars) <- c("linf","vbk","t0","CVlen")
LCpars
#          Estimate  Std. Error
# linf  80.23134337 1.760988083
# vbk    0.14025094 0.008003665
# t0    -0.83162114 0.077529320
# CVlen  0.09233058 0.002947301
LCpars.est <- LCpars[,1]
saveRDS(LCpars.est, file.path(res_dir,"growth_REisland_LCpars.rds"))

# overall/global growth pars from m2
globalpars <- summary(mods[[2]]$sdrep)[rownames(summary(mods[[2]]$sdrep))=="theta_global",]
globalpars <- rbind(globalpars, summary(mods[[2]]$sdrep)[rownames(summary(mods[[2]]$sdrep))=="CV_L",])
rownames(globalpars) <- c("linf","vbk","t0","CVlen")
globalpars <- as.data.frame(globalpars)
globalpars$low = apply(globalpars, 1, function(x) x[1]-1.96*x[2])
globalpars$high = apply(globalpars, 1, function(x) x[1]+1.96*x[2])
#          Estimate  Std. Error         low        high
# linf  80.23134337 1.760988083 76.77980673 83.68288001
# vbk    0.15450038 0.010570166  0.13378285  0.17521790
# t0    -0.83162114 0.077529320 -0.98357860 -0.67966367
# CVlen  0.09233058 0.002947301  0.08655387  0.09810729
globalpars.est <- globalpars[,1]
saveRDS(globalpars.est, file.path(res_dir,"growth_REisland_globalpars.rds"))

# Island-specific growth pars
theta <- summary(mods[[2]]$sdrep)[rownames(summary(mods[[2]]$sdrep))=="theta",]
theta <- as.data.frame(theta)
theta$low = apply(theta, 1, function(x) x[1]-1.96*x[2])
theta$high = apply(theta, 1, function(x) x[1]+1.96*x[2])
round(theta[c(2,5,8),],3)
#         Estimate Std. Error   low  high
# theta.1    0.160      0.009 0.143 0.178  Brac
# theta.4    0.164      0.009 0.146 0.182  Grand
# theta.7    0.140      0.008 0.125 0.156  Little

# ---------------------------------------------------------------------
# Fig. 3
ages <- seq(0,30,0.5)
pred_global <- globalpars[1,1]*(1-exp(-globalpars[2,1]*(ages-globalpars[3,1])))
pred_LC <- theta[7,1]*(1-exp(-theta[8,1]*(ages-theta[9,1])))
pred_CB <- theta[1,1]*(1-exp(-theta[2,1]*(ages-theta[3,1])))
pred_GC <- theta[4,1]*(1-exp(-theta[5,1]*(ages-theta[6,1])))
cols <- RColorBrewer::brewer.pal(3,"Dark2") # Little, Brac, Grand
cols2 <- c("black",cols)

cairo_pdf(filename="plots/fig3_growth.pdf", width=7, height=7)
par(mar=c(4.1,4.1,0.5,0.5))
plot(0, type='n', xlim=c(0,30),ylim=c(0,95), xlab="Age", ylab="Total length (cm)",bty='l')
polygon(y=c(pred_global-1.96*pred_global*globalpars[4,1], rev(pred_global+1.96*pred_global*globalpars[4,1])), x=c(ages, rev(ages)), col=alpha("black",0.08), border=NA)
points(dat_growth$Age[island==2], dat_growth$Length_cm[island==2], col=alpha(cols2[4],0.4), pch=19)
points(dat_growth$Age[island==1], dat_growth$Length_cm[island==1], col=alpha(cols2[3],0.4), pch=19)
points(dat_growth$Age[island==3], dat_growth$Length_cm[island==3], col=alpha(cols2[2],0.4), pch=19)
lines(ages, pred_GC, lwd=4, col=alpha(cols2[4],0.7), lty=1)
lines(ages, pred_CB, lwd=4, col=alpha(cols2[3],0.7), lty=1)
lines(ages, pred_LC, lwd=4, col=alpha(cols2[2],0.7), lty=1)
lines(ages, pred_global, lwd=3, col=cols2[1], lty=3)
legend(x=-1,y=98, legend=c("All islands","Little Cayman","Cayman Brac","Grand Cayman"), 
	col=cols2, lty=c(3,1,1,1), lwd=c(2,4,4,4), box.lty=0, bg="transparent")
dev.off()

# pdfjam --outfile fig3_growth_resize.pdf --papersize '{85mm,85mm}' fig3_growth.pdf
