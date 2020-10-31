# Figure 7 
# LBSPR using island-specific growth (from fig3_growth_REisland.R, results saved as results/growth_REisland.rds)

dir <- here()
setwd(dir)

library(tidyverse)
library(LBSPR)
library(ggsidekick)

# ----------------------------------------------
# get growth pars from m2 (k varies by island)
# use global/mean pars for Pickle Bank
mods <- readRDS("results/revision_v2/growth_REisland.rds")
vbk <- summary(mods[[2]]$sdrep)[rownames(summary(mods[[2]]$sdrep))=="theta",1][c(2,5,8)]
names(vbk) <- c("Cayman Brac","Grand Cayman","Little Cayman")
meanpars <- summary(mods[[2]]$sdrep)[rownames(summary(mods[[2]]$sdrep))=="theta_global",1]
names(meanpars) <- c("linf","vbk","t0")
CVlen <- summary(mods[[2]]$sdrep)[rownames(summary(mods[[2]]$sdrep))=="CV_L",1]

# ----------------------------------------------------
# read in catch data
dat <- read.csv("data/CaymansNG_DOE_biological_data.csv", header=T)

# get length-weight pars
fit.wl <- lm(log(Weight_kg) ~ log(Length_cm), data=dat)
lwa = exp(fit.wl$coefficients[1]) # 3.726e-06
lwb = fit.wl$coefficients[2] # 3.385

dat <- filter(dat, Age > 2 | is.na(Age)) # Remove 0-1 yo fish from assessment
dat <- select(dat, Year, Island, Length_cm) # only need length, year, island

# group island as factor with these levels
islands <- c("Little Cayman 1978", "Little Cayman 1987-1995", "Little Cayman 2002", "Cayman Brac", "Grand Cayman", "Pickle Bank")
dat$Island <- as.character(dat$Island)
dat$Island[dat$Island=="PickleBank"] <- "Pickle Bank"
dat$Island[dat$Island=="Grand"] <- "Grand Cayman"
dat$Island[dat$Island=="Little" & dat$Year < 2001] <- "Little Cayman 1987-1995"
dat$Island[dat$Island=="Little" & dat$Year > 2001] <- "Little Cayman 2002"
dat$Island[dat$Island=="Brac"] <- "Cayman Brac"
dat$Island <- factor(dat$Island, levels=islands)

# Add Pat Colin's data from 1978 EELC Jan 24 
# Standard length, need to convert to TL (so kept in separate data file)
dat78 <- read.csv("data/NG_pre2001_biological_data_1978only.csv", header=T)
dat78$Length_cm = (2.24 + 1.11 * dat78$SL*10)/10 # SL-TL equation (in mm) from Sadovy & Eklund p.31 (Cuba, Claro 1990, largest sample size and closest to Caymans)
dat78$Island <- as.character(dat78$Island)
dat78$Island[dat78$Island=="Little"] <- "Little Cayman 1978"
dat78 <- select(dat78, Year, Island, Length_cm)
dat <- rbind(dat, dat78)

# save each island/group as its own csv
for(i in 1:length(islands)){
  thedat <- filter(dat, Island==islands[i]) %>% select(Length_cm)
  write.table(thedat, file=paste0("data/df_",islands[i],".csv"), sep=",",col.names=FALSE,quote=FALSE,row.names=FALSE)
}

# -----------------------------------------------------
# Fit LBSPR by island to catch data before protection
M50 = (2.24 + 1.11 * 425)/10 # 425 and 500 mm SL from Sadovy & Eklund p.17
M95 = (2.24 + 1.11 * 500)/10 # SL-TL equation from Sadovy & Eklund p.31 (Cuba, Claro 1990, largest sample size and closest to Caymans)
sigmaF <- 0.3
sigmaC <- 0.2
sigmaR <- 1
M = 0.276 # from fig7a_estimate_M_PickleBank.R
l.min = 1
l.max = 100
binwidth = 1

MyPars <- new("LB_pars")
MyPars@Linf <- meanpars[1]
MyPars@L50 <- M50 # 425 and 500 mm SL from Sadovy & Eklund p.17
MyPars@L95 <- M95 # SL-TL equation from Sadovy & Eklund p.31 (Cuba, Claro 1990, largest sample size and closest to Caymans)
MyPars@MK <- M/meanpars[2]
MyPars@M <- M
MyPars@BinWidth <- binwidth
MyPars@L_units <- "cm"
MyPars@BinMax <- l.max
MyPars@BinMin <- l.min
MyPars@CVLinf = CVlen
MyPars@Walpha = lwa
MyPars@Wbeta = lwb

# Using island-specific growth parameters for Little, Brac, Grand
MyPars_LC <- MyPars_CB <- MyPars_GC <- MyPars
MyPars_LC@MK <- M/vbk[3]
MyPars_CB@MK <- M/vbk[1]
MyPars_GC@MK <- M/vbk[2]

# Fit LBSPR
df.SPR <- data.frame(SPR=NULL, SPR.low=NULL, SPR.high=NULL, fitLog=NULL, nll=NULL, SL50=NULL, SL95=NULL, Island=NULL, ProtectedYrs=NULL, Yrs=NULL)
for(i in 1:length(islands)){
  if(i < 4) thepars <- MyPars_LC # use island-specific growth parameters for Little, Brac, Grand
  if(i == 4) thepars <- MyPars_CB
  if(i == 5) thepars <- MyPars_GC
  if(i == 6) thepars <- MyPars # use global/mean pars for Pickle Bank
  Len1 <- new("LB_lengths", dataType="raw", LB_pars=thepars, file=paste0("data/df_",islands[i],".csv"))
  myFit1 <- LBSPRfit(thepars, Len1)
  df.SPR <- rbind(df.SPR, data.frame(SPR=myFit1@SPR, 
                                    SPR.low=myFit1@SPR - 1.96 * sqrt(myFit1@Vars[1,4]),
                                    SPR.high=myFit1@SPR + 1.96 * sqrt(myFit1@Vars[1,4]),
                                    fitLog=myFit1@fitLog,
                                    nll=myFit1@NLL,
                                    SL50=myFit1@SL50,
                                    SL95=myFit1@SL95,
                                    Island=c(rep("Little Cayman",3),"Cayman Brac","Grand Cayman","Pickle Bank")[i],
                                    ProtectedYrs=c("Pre A","Pre B","Pre C",rep("Pre",3))[i],
                                    Yrs=c("1978", "1987-1995", "2002", "1990-2000","1988-1997","2000")[i]))
}
df.SPR$Island <- factor(c("Little Cayman", "Little Cayman", "Little Cayman", "Cayman Brac", "Grand Cayman", "Pickle Bank"),
                        levels=c("Little Cayman", "Cayman Brac", "Grand Cayman", "Pickle Bank"))

# -------------------------------------------------------
# Estimate SPR after protection (1-5, 6-10, 11-15 years pooled), Little Cayman + Cayman Brac
x <- read.csv("data/NG_LC_lengths_2003_2019.csv")
dat.LC <- x %>% filter(usable!="n")
dat.LC$Island <- "Little Cayman"
dat.LC$ProtectedYrs <- cut(dat.LC$year, breaks=c(2002.5,2007.5,2012.5,2016.5,2019.5))
levels(dat.LC$ProtectedYrs) <- c("0-4", "5-9", "10-13","14-16")
dat.LC$Length <- dat.LC$length/10 # convert mm to cm

# Read in Cayman Brac data
dat.CB <- read.csv("data/NG_CB_lengths_2017_2019.csv")
dat.CB$ProtectedYrs <- "14-16"

# Merge Little and Brac data
dat.lasers <- rbind(dat.LC[,c("Length", "Island","ProtectedYrs")],
                    dat.CB[,c("Length", "Island","ProtectedYrs")])
tmp <- dat.lasers %>% group_by(Island, ProtectedYrs) %>% summarize(n=n()) %>% as.data.frame
tmp$Yrs = c("2017-2019", "2003-2007","2008-2012","2013-2016","2017-2019")

for(i in 1:dim(tmp)[1]){
  thedat <- filter(dat.lasers, Island==tmp[i,"Island"] & ProtectedYrs==tmp[i,"ProtectedYrs"]) %>% select(Length)
  write.table(thedat, file=paste0("data/tmp_",i,".csv"), sep=",",col.names=FALSE,quote=FALSE,row.names=FALSE)  
  Len1 <- new("LB_lengths", dataType="raw", LB_pars=MyPars, file=paste0("data/tmp_",i,".csv"))
  if(tmp[i,"Island"] == "Cayman Brac") thepars <- MyPars_CB else thepars <- MyPars_LC
  myFit1 <- LBSPRfit(MyPars, Len1)
  df.SPR <- rbind(df.SPR, data.frame(SPR=myFit1@SPR, 
                                    SPR.low=myFit1@SPR - 1.96 * sqrt(myFit1@Vars[1,4]),
                                    SPR.high=myFit1@SPR + 1.96 * sqrt(myFit1@Vars[1,4]),
                                    fitLog=myFit1@fitLog,
                                    nll=myFit1@NLL,
                                    SL50=myFit1@SL50,
                                    SL95=myFit1@SL95,
                                    Island=tmp[i,"Island"],
                                    ProtectedYrs=tmp[i,"ProtectedYrs"],
                                    Yrs=tmp[i,"Yrs"]))
}
df.SPR$Island <- dplyr::recode(df.SPR$Island, `Grand Cayman`="Grand", `Pickle Bank`="Pickle")
df.SPR$Yrs <- as.character(df.SPR$Yrs)
df.SPR$Yrs <- factor(df.SPR$Yrs, levels=c("1978","1987-1995","1990-2000","1988-1997","2000","2002","2003-2007","2008-2012","2013-2016","2017-2019"))
df.SPR$border <- as.character(ifelse(df.SPR$ProtectedYrs %in% c("Pre","Pre A","Pre B","Pre C"), 0,1))
addline_format <- function(x,...){
    gsub('-','–\n',x)
}

png("plots/fig7_SPR.png",res=300,units='in',width=8,height=4.5)
ggplot(df.SPR, aes(x=Yrs, y=SPR, fill=Island, group=interaction(Island,Yrs))) +
  geom_hline(yintercept=0.4, linetype=2) +
  geom_linerange(aes(ymin=SPR.low,ymax=SPR.high), size=0.6) +  
  geom_point(aes(fill=Island,y=SPR,shape=border), size=5) +   
  scale_shape_manual(values=c(21,22), guide=FALSE) +
  scale_fill_brewer(palette="Set2", guide=FALSE) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6), limits=c(0.2,1.08)) +
  scale_x_discrete(breaks=unique(df.SPR$Yrs), 
    labels=addline_format(unique(df.SPR$Yrs))) +
  ylab("Spawning potential ratio (SPR)") +
  xlab("Year") +
  facet_grid (.~ Island, scales = "free_x", space = "free_x") +
  theme_bw() +
  # theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), 
  theme(panel.grid.major = ggplot2::element_blank(),   
    panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
    axis.line = ggplot2::element_line(colour = "black"), axis.title=ggplot2::element_text(size=16), 
    axis.text=ggplot2::element_text(size=12), legend.text=ggplot2::element_text(size=14),
    strip.text=ggplot2::element_text(size=14))  
dev.off() 

cbind(round(df.SPR[,1:3],3),df.SPR[,6:11])
#         SPR SPR.low SPR.high  SL50  SL95        Island ProtectedYrs       Yrs
# SPR   0.446   0.278    0.613 52.91 60.51 Little Cayman        Pre A      1978
# SPR1  0.484   0.449    0.519 54.29 61.17 Little Cayman        Pre B 1987-1995
# SPR2  1.000   1.000    1.000 50.34 53.72 Little Cayman        Pre C      2002
# SPR3  0.498   0.438    0.557 52.37 61.16   Cayman Brac          Pre 1990-2000
# SPR4  0.529   0.484    0.574 58.50 67.68         Grand          Pre 1988-1997
# SPR5  0.731   0.589    0.873 72.98 83.84        Pickle          Pre      2000
# SPR6  1.000   1.000    1.000 58.49 72.52   Cayman Brac        14-16 2017-2019
# SPR7  0.902   0.846    0.959 54.52 62.83 Little Cayman          0-4 2003-2007
# SPR8  0.596   0.559    0.633 57.60 68.65 Little Cayman          5-9 2008-2012
# SPR9  0.717   0.684    0.750 60.47 70.61 Little Cayman        10-13 2013-2016
# SPR10 0.937   0.863    1.010 49.03 55.18 Little Cayman        14-16 2017-2019
