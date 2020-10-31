# Little Cayman Nassau Grouper length assessment
#   LIME-integrated model (fit growth in model)
#   Dirichlet-multinomial likelihood for length composition data

dir <- here()
setwd(dir)
res_dir <- file.path(dir, "results","LIME_integrated_dirichlet")
out_dir <- file.path(res_dir,"output")
plots_dir <- file.path(res_dir,"LIME_plots")
fits_dir <- file.path(res_dir,"LCfits")
residuals_dir <- file.path(res_dir,"residuals")
if(!file.exists(res_dir)) dir.create(res_dir, showWarnings=FALSE)
if(!file.exists(plots_dir)) dir.create(plots_dir, showWarnings=FALSE)
if(!file.exists(fits_dir)) dir.create(fits_dir, showWarnings=FALSE)
if(!file.exists(residuals_dir)) dir.create(residuals_dir, showWarnings=FALSE)
if(!file.exists(out_dir)) dir.create(out_dir, showWarnings=FALSE)

library(viridis)
library(ggsidekick) # devtools::install_github("seananderson/ggsidekick")
library(LIME)
library(tidyverse)
library(TMB)
compile("code/LIME_integrated_REisland.cpp","-O0 -g")
dyn.load("code/LIME_integrated_REisland.so")

# helper functions, overwrite some LIME functions with modifications
source("code/plot_LIME_output.R")
source("code/plot_stdresiduals.R")
source("code/plot_LCfits.R")
source("code/create_lh_list.R")
source("code/fit_tmb.R")

# --------------------------------------------------------
# load data
# convert lengths in inches to cm
dat <- read.csv("data/NG_pre2001_biological_data_2002.csv", header=T)
dat$Length_cm <- dat$Length
dat$Length_cm[dat$Length_units == "in"] <- dat$Length[dat$Length_units == "in"] * 2.54

x <- read.csv("data/NG_LC_lengths_2003_2019.csv")
x <- x %>% filter(usable!="n")

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
#   best assigned age = ring count + 0.5, as they were caught in July and Aug, 6 months after spawning
#   all others were caught on SPAG in Jan-Feb, so add 1 year to ring count (see 1988 Table 1 note)
dat_all$Age <- ifelse(dat_all$Age > 2, dat_all$Age + 1, dat_all$Age + 0.5)

# add 2012 YOY (1 year olds, used in growth curve but not length analysis)
x <- c(8, 16.4, 15, 15, 23.25, 16.8, 14.6, 19.9, 12.3, 15.1, 14.7, 15.7, 18.6, 17.4, 22, 23.1, 22, 20, 19.1, 15, 19.9, 20.5, 20.1, 15, 15)
dat_yoy <- data.frame(matrix(NA,ncol = length(colnames(dat_all)), nrow = length(x)))
colnames(dat_yoy) <- colnames(dat_all)
dat_yoy$Length_cm = x
dat_yoy$Year = 2012
dat_yoy$Age = 1
dat_yoy$Island = "Little"
dat_all <- rbind(dat_all, dat_yoy)

dat_growth <- filter(dat_all, !is.na(dat_all$Age) & !is.na(dat_all$Length_cm)) # for growth curve, leave in the 1-2 yr-old fish
dat_growth$Island[dat_growth$Island=="Grand_12MileBanks"] = "Grand"
dat_growth_fit <- select(dat_growth, Age, Length_cm) %>% as.matrix
island <- sapply(dat_growth$Island, match, c("Brac", "Grand","Little"))
n.islands <- length(unique(island))

dat_all <- filter(dat_all, Age > 2 | is.na(Age)) # for length analysis, remove 1-2 yr-old fish (not caught on SPAG, summer instead of winter)
dat_all$Year[dat_all$Island=="PickleBank"] <- 2000
dat_LC <- filter(dat_all, Island=="Little")

# set up length-frequency data
years <- as.numeric(levels(dat_LC$Year)) # Years go from 1987 to 2018
l.min = 0
l.max = 100
binwidth = 1
LC_LF <- dat_LC %>% 
			mutate(bin=cut(Length_cm, breaks=seq(l.min+binwidth/2,l.max+binwidth/2,binwidth), labels=seq(l.min+binwidth, l.max, binwidth))) %>%
		   	group_by(Year) %>%
		   	count(bin) %>%
		   	tidyr::complete(bin) %>% 
		   	tidyr::spread(bin,n)
LC_LF[is.na(LC_LF)] <- 0
yr.labs <- unlist(lapply(LC_LF[1], as.character))
LC_LF <- as.matrix(LC_LF[,-1])
rownames(LC_LF) <- yr.labs

fit.wl <- lm(log(Weight_kg) ~ log(Length_cm), data=dat_all)

# save parameters from length - weight curve
lwa = exp(fit.wl$coefficients[1]) # 3.726e-06
lwb = fit.wl$coefficients[2] # 3.385
# initialize selectivity = maturity (all fish caught/measured on SPAG)
L50 = (2.24 + 1.11 * 425)/10 # 425 and 500 mm SL from Sadovy & Eklund p.17
L95 = (2.24 + 1.11 * 500)/10 # SL-TL equation from Sadovy & Eklund p.31 (Cuba, Claro 1990, largest sample size and closest to Caymans)
S50 = M50 = L50
S95 = M95 = L95

# --------------------------------------------------------
# use MLE growth pars as starting values in assessment
res <- readRDS("results/growth_REisland_globalpars.rds")
vbk = res[2]
linf = res[1]
t0 = res[3]

# get index of abundance
# as in Waterhouse et al. 2020 except modified so that yearly estimates are independent
LC.ind <- read.csv("data/NG_LC_index_indep.csv")
# 2007 is missing, add in as -99
r = 2
LC.ind <- rbind(LC.ind[1:r,], c(2007, -99, -99, -99), LC.ind[-(1:r),])
med <- LC.ind$median
low <- LC.ind$q0.025
upp <- LC.ind$q0.975

sigI_vec <- c(.1, .125, .15, .175, .2)
sigC_vec <- c(.05,.1,.2)

# other parameters
sigmaF <- 0.3
sigmaR <- 1

# M estimate from catch curve analysis of Pickle Bank "unfished" length distribution 
#   code/fig7a_estimate_M_PickleBank.R
M = 0.276
yrs <- 1999:2019
n.yrs <- length(yrs)

# CVlen = 0.096 # best from likelihood profile of CV_L
CVlen_vec <- seq(from=0.0915, to=0.1005, by=0.0005)

# for(cvi in 1:length(CVlen_vec)){ # use loop to do likelihood profile over CV_L
cvi=10
	output <- list()
	lh <- create_lh_list(vbk=vbk,
					 linf=linf,
					 t0=t0,
					 lwa=lwa,
					 lwb=lwb,
					 S50=L50, 
					 S95=L95, 
					 selex_input="length",
					 M50=L50,
					 M95=L95,
					 maturity_input = "length",
					 M=M,
					 binwidth=binwidth,		 
					 CVlen=CVlen_vec[cvi],
					 AgeMax=NULL,
					 SigmaR=sigmaR,
					 SigmaF=sigmaF,
					 SigmaC=sigmaC,
					 SigmaI=sigmaI,
					 R0=6000,	
					 qcoef=1,
					 start_ages=0,
					 rho=0.43,
					 nseasons=1,
					 nfleets=1)

	# Use 2002 FSA catch, but should be in model year 2001 (before spawning #s 2002)
	LF <- rbind(matrix(rep(0,2*(l.max-l.min)), nrow=2), # 1999-2000
				LC_LF[as.numeric(rownames(LC_LF)) == 2002,], # 2001 slot
				matrix(rep(0,1*(l.max-l.min)), nrow=1), # 2002
				LC_LF[as.numeric(rownames(LC_LF)) > 2002,]) # 2003-2019
	dimnames(LF) <- list(yrs,colnames(LC_LF))

	index.video <- matrix(nrow=1,ncol=length(yrs),
		data=c(-99,-99,-99,-99,-99,-99, round(med), -99)) 
	colnames(index.video) <- yrs

	catch <- matrix(nrow=1,ncol=length(yrs),
		data=c(-99,-99,2000,1934,-99,-99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99, -99)) 
	colnames(catch) <- yrs

	data_lasers <- list(years=yrs,
						LF=LF,
						I_ft=index.video,
						C_ft=catch)		
	inputs_all <- create_inputs(lh=lh, input_data=data_lasers)

	true_years <- yrs

	rdev <- rep(0,n.yrs)
	est_F_ft <- matrix(1,nrow=1,ncol=n.yrs) # est F in all years except 1999 and 2000 (fix at 0)
	fixF <- c(1,2,5,6,15,16,21)
	est_F_ft[1,fixF] = 0

	f_startval_ft <- matrix(0.05,nrow=1,ncol=n.yrs)
	f_startval_ft[1,3:4] <- 1 # high F expected in 2001-2002
	f_startval_ft[1,fixF] = 0

	# ==========================================================================
	# Prep LIME
	C_type = 1 # catch is in numbers
	est_more = FALSE
	fix_more = c("log_q_f")
	rdev_startval_t = rdev
	F_up=10 # default value	

	Sdreport <- NA
	ParList <- NA  
	df <- NULL
	input <- inputs_all

	data_avail <- "Index_Catch_LC"
	Fpen=1
	SigRpen=1
	SigRprior=c(0.737,0.3)
	LFdist=1 # Dirichlet-Multinomial (default)
	est_selex_f=TRUE
	vals_selex_ft=-1

	est_rdev_t = rep(1,n.yrs)
	est_rdev_t[(n.yrs-3):n.yrs] = 0 # fix last 4 rdevs at 0 (not estimable)

	newtonsteps=3
	S50_up=NULL
	est_totalF=FALSE
	prop_f=1

	if(all(vals_selex_ft < 0)) input$vals_selex_ft <- matrix(-1, nrow=input$nfleets, ncol=length(input$highs))
	if(any(vals_selex_ft >= 0)) input$vals_selex_ft <- vals_selex_ft
	if(all(prop_f==1)) prop_f_inp <- rep(1/input$nfleets, input$nfleets)

	tmblist <- with(input, {
		S_yrs_inp <- 1:Nyears
		selex_type_f <- rep(1,nfleets)
		for(i in 1:nfleets) if(any(vals_selex_ft[i,] > 0)) selex_type_f[i] <- 0
	    if(all(is.null(input$neff_ft))){
	        if(is.vector(LF[,,1])==FALSE) n_inp <- t(sapply(1:nfleets, function(x) rowSums(LF[,,x])))
	        if(is.vector(LF[,,1])) n_inp <- t(sapply(1:nfleets, function(x) sum(LF[,,x])))
	    }
	    if(all(is.null(input$neff_ft)==FALSE)) n_inp <- input$neff_ft
		if(all(est_rdev_t == TRUE)){ est_rdev_t_inp <- rep(1, Nyears) } else { est_rdev_t_inp <- est_rdev_t }

		ML <- apply(LF[,,1], 1, function(x) sum(x*mids)/sum(x))
		haveL <- rep(1,length(ML))
		haveL[is.nan(ML)] = 0
		ML[!haveL] = 0
		Data <- list("n_t"=dim(LF)[1],
	                 "n_lb"=dim(LF)[2],
	                 "n_fl"=dim(LF)[3],
	                 "n_a"=length(ages),
	                 "LF_tlf"=LF,
	                 "n_LF_ft"=n_inp,
	                 "I_ft"=I_ft,
	                 "C_ft"=C_ft,
	                 "C_type"=C_type,
	                 "ML_t"=ML,
	                 "ML_years"=haveL,		                 
	                 "ages"=ages,
	                 "match_ages"=seq(min(ages), max(ages), by=1),
	                 "lw_pars"=c(lwa, lwb),
	                 "Mat_l"=Mat_l,
	                 "M"=M,
	                 "h"=h,
	                 "lbhighs"=highs,
	                 "lbmids"=mids,
	                 "Fpen"=Fpen,
	                 "SigRpen"=SigRpen,
	                 "SigRprior"=SigRprior,
	                 "selex_type_f"=selex_type_f,
	                 "vals_selex_ft"=vals_selex_ft,
	                 "LFdist"=LFdist,
	                 "S_yrs"=S_yrs_inp,
	                 "n_s"=nseasons,
	                 "n_y"=Nyears,
	                 "mirror_theta"=0,
	                 "mirror_q"=0,
	                 "est_totalF"=ifelse(est_totalF==TRUE,1,0),
	                 "prop_f"=prop_f,
	                 "dat_growth"=dat_growth_fit,
	                 "n_islands"=n.islands,
	                 "island"=island) 
		
		if(all(is.null(rdev_startval_t))) rdev_startval_t <- rep(0, Nyears)
	    Parameters <- list("log_F_ft"=log(f_startval_ft),
	                    "log_q_f"=rep(log(qcoef), Data$n_fl),
	                    "beta"=log(R0),
	                    "log_sigma_R"=log(SigmaR),
	                    "log_S50_f"=log(SL50),
	                    "log_Sdelta_f"=log(SL95 - SL50), 
	                    "log_sigma_F"=log(SigmaF), 
	                    "log_sigma_C"=log(SigmaC),
	                    "log_sigma_I"=log(SigmaI),
	                    "log_CV_L"=log(CVlen_vec[cvi]),
	                    "log_theta"=log(rep(theta, Data$n_fl)),
	                    "growth_mu"=c(linf, log(vbk), t0),
	                    "log_sigma_mu"=c(rep(log(0.1),n.islands)),
	                    "rho_untrans"=rep(0,n.islands),
	                    "emat"=matrix(0,nrow=3,ncol=n.islands),
	                    "Nu_input"=rdev_startval_t)

	    Map = list()

	    if(grepl("Catch",data_avail)==FALSE){        
	        Map[["beta"]] <- NA
	        Map[["beta"]] <- factor(Map[["beta"]])
	    }
	    if(all(fix_more!=FALSE)){
	        for(i in 1:length(fix_more)){
	            Map[[fix_more[i]]] <- rep(NA, length(Parameters[[fix_more[i]]]))
	            Map[[fix_more[i]]] <- factor(Map[[fix_more[i]]])
	        }
	    }
	    if("log_CV_L" %in% est_more==FALSE){
	        Map[["log_CV_L"]] <- NA
	        Map[["log_CV_L"]] <- factor(Map[["log_CV_L"]])
	    }
	    if(all(est_more==FALSE)){
	        Map[["log_sigma_F"]] <- NA
	        Map[["log_sigma_F"]] <- factor(Map[["log_sigma_F"]])
	        Map[["log_CV_L"]] <- NA
	        Map[["log_CV_L"]] <- factor(Map[["log_CV_L"]]) 
	        Map[["log_sigma_C"]] <- NA
	        Map[["log_sigma_C"]] <- factor(Map[["log_sigma_C"]])  
	        Map[["log_sigma_I"]] <- NA
	        Map[["log_sigma_I"]] <- factor(Map[["log_sigma_I"]])
	    }
	    if(LFdist==0){
	        Map[["log_theta"]] <- rep(NA, length(Parameters[["log_theta"]]))
	        Map[["log_theta"]] <- factor(Map[["log_theta"]])
	    }

	    if(any(est_selex_f == FALSE)){
	        Map[["log_S50_f"]] <- Parameters$log_S50_f
	        if(all(est_selex_f==FALSE) & nfleets > 1) Map[["log_S50_f"]][which(est_selex_f==FALSE)] <- NA
	        if(all(est_selex_f==FALSE)) Map[["log_S50_f"]] <- rep(NA, length(Parameters[["log_S50_f"]]))
	        Map[["log_S50_f"]] <- factor(Map[["log_S50_f"]])

	        Map[["log_Sdelta_f"]] <- Parameters$log_Sdelta_f
	        if(all(est_selex_f==FALSE) & nfleets > 1) Map[["log_Sdelta_f"]][which(est_selex_f==FALSE)] <- NA
	        if(all(est_selex_f==FALSE)) Map[["log_Sdelta_f"]] <- rep(NA, length(Parameters[["log_Sdelta_f"]]))
	        Map[["log_Sdelta_f"]] <- factor(Map[["log_Sdelta_f"]])
	    }
	    if(any(vals_selex_ft >= 0)){
	        Map[["log_S50_f"]] <- Parameters$log_S50_f
	        Map[["log_S50_f"]][which(vals_selex_ft[,1] >= 0)] <- NA
	        Map[["log_S50_f"]] <- factor(Map[["log_S50_f"]])

	        Map[["log_Sdelta_f"]] <- Parameters$log_Sdelta_f
	        Map[["log_Sdelta_f"]][which(vals_selex_ft[,1] >= 0)] <- NA
	        Map[["log_Sdelta_f"]] <- factor(Map[["log_Sdelta_f"]])                
	    }
	    if(any(est_rdev_t_inp == 0)){
	        Map[["Nu_input"]] <- seq_along(Parameters$Nu_input)
	        Map[["Nu_input"]][which(est_rdev_t_inp == 0)] <- NA
	        Map[["Nu_input"]] <- factor(Map[["Nu_input"]])
	    }
	    if(all(est_rdev_t_inp == 0)) Map[["log_sigma_R"]] <- factor(NA)
	    if(any(est_F_ft == 0)){
	        fstart <- Parameters$log_F_ft #+ matrix(rnorm(length(Parameters$log_F_ft), 0, 2), nrow=nrow(Parameters$log_F_ft), ncol=ncol(Parameters$log_F_ft))
	        for(i in 1:nfleets){
	            fstart[i,which(est_F_ft[i,]==0)] <- NA
	            fstart[i,which(est_F_ft[i,]==1)] <- 1:length(which(est_F_ft[1,]==1))
	        }
	        Map[["log_F_ft"]] <- factor(fstart)

	    }
	    if(length(Map)==0) Map <- NULL

	    # only allow k to vary by island
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

	    Random <- c("Nu_input","emat")
		TmbList <- list("Parameters"=Parameters, "Data"=Data, "Random"=Random, "Map"=Map)
		return(TmbList)
	})
	output$Inputs <- tmblist
	obj <- TMB::MakeADFun(data=tmblist[["Data"]], parameters=tmblist[["Parameters"]], random=tmblist[["Random"]], map=tmblist[["Map"]], inner.control=list(maxit=1e3), DLL="LIME_growth_REisland", checkParameterOrder = FALSE)

	## Settings
	Upr = rep(Inf, length(obj$par))
	Upr[match("log_sigma_R",names(obj$par))] = log(2)
	if(is.null(S50_up)==FALSE) Upr[which(names(obj$par)=="log_S50_f")] <- log(S50_up)
	if(is.null(S50_up)) Upr[which(names(obj$par)=="log_S50_f")] <- tmblist$Parameters$growth_mu[1]
	Upr[which(names(obj$par)=="log_F_ft")] = log(F_up)
	Upr[match("log_sigma_F", names(obj$par))] <- log(2)

	Lwr <- rep(-Inf, length(obj$par))
	Lwr[match("log_CV_L",names(obj$par))] = log(0.001)
	Lwr[match("log_sigma_C",names(obj$par))] = log(0.001)
	Lwr[match("log_sigma_I",names(obj$par))] = log(0.001) 
	Lwr[which(names(obj$par)=="log_S50_f")] = log(1)

	opt <- fit_tmb(obj, n.newton=newtonsteps, do.sdrep=TRUE, do.check=TRUE)
		output$Report <- tryCatch( obj$report(), error=function(x) NA)
		output$Sdreport <- tryCatch(sdreport(obj, bias.correct=TRUE), error=function(x) NA )		
		output$obj <- obj
		output$opt <- opt	
		output$years <- yrs	

		saveRDS(output, file.path(out_dir, paste0("output_CV",CVlen_vec[cvi],"_integrated.rds")))

	if(!is.null(opt)){
		png(paste0(plots_dir, "/integrated_dirichlet_CV",CVlen_vec[cvi],"_LIME_plots.png"),res=300,units='in',width=11,height=5.5)
		plot_output(Inputs=output$Inputs, Report=output$Report, Sdreport=output$Sdreport, true_years=yrs,
			plot=c("ML","Rec","Ind","SB","Fish","Selex"), set_ylim=list("Ind" = c(0,9500), "Rec"=c(0,40000)), relative=FALSE, est=list(Rec=est_rdev_t, Fish=est_F_ft))
		dev.off()

		png(paste0(fits_dir,"/integrated_dirichlet_CV",CVlen_vec[cvi],"_LIME_LCfit.png"),width=10,height=7,units="in",res=300)
		g <- plot_LCfits(Inputs=output$Inputs, Report=output$Report, n=TRUE, year_labels=yrs)
		print(g + theme_sleek() + theme(legend.position='none'))
		dev.off()
	}
# } # end CV_L loop if used
