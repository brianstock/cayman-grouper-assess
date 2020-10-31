# Figure 6
# Compare bimodal length distributions between Little Cayman and Cayman Brac, 2017-2019 

dir <- here()
setwd(dir)

library(TMB)
library(tidyverse)
library(ggridges)
library(zoo)
library(RColorBrewer)
library(ggthemes)

# --------------------------------------------------------------
# Read in data
# Little Cayman length data
x <- read.csv("data/NG_LC_lengths_2003_2019.csv")
dat.LC <- x %>% filter(usable!="n") %>% filter(year %in% 2017:2019)
dat.LC$year <- factor(dat.LC$year)
dat.LC$Island <- "Little Cayman"
dat.LC <- data.frame(Length=dat.LC$length/10, Island=dat.LC$Island, Year=dat.LC$year)

# Read in Cayman Brac data
dat.CB <- read.csv("data/NG_CB_lengths_2017_2019.csv")
dat <- rbind(dat.LC, dat.CB)

# ---------------------------------------------------------------
# find bimodal maxima -- update to use different 'adjust' param values
get_max <- function(DAT, YEAR, ADJ){
  dens.18 <- density(subset(DAT, Year==YEAR)$Length, adjust=ADJ)
  zoo.18 <- as.zoo(dens.18$y)
  rxz <- rollapply(zoo.18, 3, function(x) which.max(x)==2)
  ind.max <- index(rxz)[coredata(rxz)]
  return(dens.18$x[ind.max])
}
mode.all <- data.frame(Length=NULL, Year=NULL, Island=NULL, adj=NULL, yend=NULL, y=NULL)

# density smoother adjust = 1 (default)
adj=1
max.CB <- c(get_max(dat.CB, YEAR=2017, ADJ=adj)[1], get_max(dat.CB, YEAR=2018, ADJ=adj)[1], get_max(dat.CB, YEAR=2019, ADJ=adj)[1])
max.LC <- c(get_max(dat.LC, YEAR=2017, ADJ=adj)[1], get_max(dat.LC, YEAR=2018, ADJ=adj)[1], get_max(dat.LC, YEAR=2019, ADJ=adj)[1])
mode.tmp <- data.frame(Length=c(rev(max.CB), rev(max.LC)),
                      Year=c(2019:2017, 2019:2017),
                      Island=c(rep("Cayman Brac",3), rep("Little Cayman",3)),
                      adj=rep(adj,6),
                      yend=c(4.4,4,2.2,5,3.5,2),
                      y=c(3.8,2,1,3.5,2,1))

# get sample size
sampsizes <- dat %>% group_by(Island, Year) %>% summarize(n=n()) %>% mutate(nlab=paste0("n = ",n)) %>% as.data.frame
sampsizes$Length = 86.5 # x-axis position

# ---------------------------------------------------------------------------
# get growth curve length estimates at ages 5-8 for LC and CB
#   calculated in fig3_growth_REisland.R, use m2, only k varies by island
# first get island-specific growth pars
mods <- readRDS("results/growth_REisland.rds")
LCpars <- summary(mods[[2]]$sdrep)[rownames(summary(mods[[2]]$sdrep))=="theta",1][7:9]
LCpars <- c(LCpars, summary(mods[[2]]$sdrep)[rownames(summary(mods[[2]]$sdrep))=="CV_L",1])
names(LCpars) <- c("linf","vbk","t0","CVlen")
CBpars <- summary(mods[[2]]$sdrep)[rownames(summary(mods[[2]]$sdrep))=="theta",1][1:3]
CBpars <- c(CBpars, summary(mods[[2]]$sdrep)[rownames(summary(mods[[2]]$sdrep))=="CV_L",1])
names(CBpars) <- c("linf","vbk","t0","CVlen")

# calculate expected length at ages 5-8
vonB <- function(params,ages){
  Linf <- params[1]
  k <- params[2]
  a0 <- params[3]
  return(Linf*(1-exp(-k*(ages-a0))))
}
ages <- 5:8
len.at.age.LC <- round(vonB(params=LCpars, ages=ages),1)
len.at.age.CB <- round(vonB(params=CBpars, ages=ages),1)
len.df <- data.frame(Length=c(len.at.age.CB, len.at.age.LC), Age=rep(ages,2), 
                    Island=c(rep("Cayman Brac",length(ages)),rep("Little Cayman",length(ages))), y=1.08)

# stacked histogram with density overlay
sampsizes$y <- c(.037,.043,.043,.043,.045,.043)
sampsizes$x_yr <- 40
len.df$y=0.002
len.df$Year=factor(2017)
mode.tmp$y = 0
mode.tmp$yend = c(.03,.044,.025,.048,.043,.027)
mode.tmp$Year <- factor(mode.tmp$Year, levels=2019:2017, labels=2019:2017)
mode.tmp$Island <- factor(as.character(mode.tmp$Island), levels=c("Little Cayman","Cayman Brac"), labels=c("Little Cayman","Cayman Brac"))
dat$Year <- factor(as.numeric(as.character(dat$Year)), levels=2019:2017, labels=2019:2017)
dat$IslYr <- interaction(dat$Year, dat$Island)

png("plots/fig6_LC_CB_2017_2019.png",width=5,height=7,units="in",res=300)
dat %>% ggplot(aes(x=Length)) +
  geom_histogram(aes(fill = Island, x=Length, y = ..density..), color = "grey50", bins=22, size=0.2,alpha=0.6, inherit.aes=FALSE) +
  geom_density(color='grey20', size=.5, alpha = 1) +  
  scale_fill_manual(values=brewer.pal(4,"Set2"), guide=FALSE) +
  scale_color_manual(values=brewer.pal(4,"Set2"), guide=FALSE) +
  geom_text(data=sampsizes, aes(x=Length, y=y, label=nlab), size = 4,inherit.aes=FALSE) +
  geom_text(data=sampsizes, aes(x=x_yr ,y=y, label=Year), hjust = 0, size=4.5) +
  facet_wrap(~Island+Year, ncol=1, scales='free_y', strip.position="right") +
  geom_segment(data = mode.tmp, aes(xend=Length, y=y, yend=yend), color = "grey20", lwd=1.5, linetype=3) +
  geom_point(data=len.df, aes(y=y, shape=factor(Age)), color="grey20",size=3) +
  guides(shape = "none") +  
  coord_cartesian(xlim=c(41, 89)) +
  scale_y_continuous(expand = c(0.01, 0.001)) +
  theme_pander() +
  theme(plot.title = element_text(size = 16),
          axis.text.x= element_text(size = 12),
          axis.title.x= element_text(size = 14),
          axis.title.y= element_text(size = 14),
          axis.text.y=element_blank(),
          axis.ticks.x=element_line(colour="black"),
          axis.ticks.y = element_blank(),
          strip.background = element_blank(),
          strip.text.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
  labs(x = "Total length (cm)", y = "Frequency") 
dev.off()
