# Figure 5
# Ridge density Little Cayman (pre2001 catch + post2001 lasers)

dir <- here()
setwd(dir)
plots_dir <- file.path(dir, "plots")

library(tidyverse)
library(ggthemes)
library(RColorBrewer)
library(viridis)
library(ggridges)
dat <- read.csv("data/NG_pre2001_biological_data_2002.csv", header=T)

dat$Island <- as.character(dat$Island)
dat$Island[dat$Island=="Grand_12MileBanks"] <- "Grand"
dat$Island[dat$Island=="PickleBank"] <- "Pickle Bank"
dat$Island[dat$Island=="Grand"] <- "Grand Cayman"
dat$Island[dat$Island=="Little"] <- "Little Cayman"
dat$Island[dat$Island=="Brac"] <- "Cayman Brac"
dat$Island <- factor(dat$Island, levels=c("Little Cayman", "Cayman Brac", "Grand Cayman", "Pickle Bank"))

# convert lengths in inches to cm
dat$Length_cm <- dat$Length
dat$Length_cm[dat$Length_units == "in"] <- dat$Length[dat$Length_units == "in"] * 2.54

# Leave in the 1 0-yr-old and 7 1-yr-old fish for growth curve
# Remove for analysis of length distribution of spawning population (not caught on SPAG, summer instead of winter)
dat <- filter(dat, Age > 1 | is.na(Age))
dat_LC_Brac <- filter(dat, Island %in% c("Little Cayman","Cayman Brac"))

# Add Pat Colin's data from 1978 EELC Jan 24 
# Standard length, need to convert to TL (so kept in separate data file)
# SL-TL equation (in mm) from Sadovy & Eklund p.31 (Cuba, Claro 1990, largest sample size and closest to Caymans)
dat78 <- read.csv("data/NG_pre2001_biological_data_1978only.csv", header=T)
dat78$Length_cm = (2.24 + 1.11 * dat78$SL*10)/10 
dat78$Island <- as.character(dat78$Island)
dat78$Island[dat78$Island=="Little"] <- "Little Cayman"

tmp <- rbind(dat_LC_Brac[,c("Length_cm","Island","Year")], dat78[,c("Length_cm","Island","Year")])
dat.catch <- data.frame(Length=tmp$Length_cm, Island=tmp$Island, Year=tmp$Year, type="Catch")

# ------------------------------------------------------------------
# read in Little Cayman laser lengths, 2003-2019
x <- read.csv("data/NG_LC_lengths_2003_2019.csv")
x <- x %>% filter(usable!="n")
dat.LC <- data.frame(Length=x$length/10, Island="Little Cayman", Year=x$year, type="Laser calipers (in situ)",stringsAsFactors=FALSE)

# --------------------------------------------------------------
# Read in Cayman Brac data
dat.CB <- read.csv("data/NG_CB_lengths_2017_2019.csv")
dat.CB$type = "Laser calipers (in situ)"

# --------------------------------------------------------
# Merge data frames
dat.catch[] <- lapply(dat.catch, as.character) # convert all columns to character
dat.LC[] <- lapply(dat.LC, as.character)
dat.CB[] <- lapply(dat.CB, as.character)
dat.all <- rbind(dat.catch, dat.LC, dat.CB)

# reconvert columns to appropriate data classes
dat.all$Length <- as.numeric(dat.all$Length)
dat.all$Year <- factor(dat.all$Year, levels=1978:2019)
dat.all$type <- factor(dat.all$type)
dat.all$Island <- factor(dat.all$Island)

# ------------------------------------------------------------------
# Ridge density Little Cayman only, pre2001 catch + lasers
sampsizes <- dat.all %>% filter(Island=="Little Cayman") %>% group_by(Year) %>% summarize(n=n()) %>% 
                         mutate(nlab=paste0("n = ",as.character(n))) %>% complete(Year, fill = list(n=0,nlab="")) %>% as.data.frame
sampsizes$Length = 85.5 # x-axis position

dev.new(width = 3.35, height = 7.81)
dat.all %>% filter(Island=="Little Cayman") %>% ggplot(aes(x=Length, y=Year, fill=type, height=..density..)) +
  geom_density_ridges_gradient(scale = 2, gradient_lwd = 1.) +  
  geom_text(data=sampsizes, aes(x=Length, y=Year, label=nlab, height=0), size = 2, vjust=-1.15, hjust=1, inherit.aes=FALSE) +
  scale_fill_manual(values=c("grey60","grey80"), guide=FALSE) +  
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = expansion(add = c(0, 2.5)), drop=FALSE) +
  labs(x = "Total length (cm)", y = "") +
  coord_cartesian(xlim=c(38, 85)) +
  theme_ridges(center=TRUE, grid = FALSE) +
  theme(axis.text= element_text(size = 8),
          axis.title.x= element_text(size = 10),
          axis.ticks.x=element_line(colour="black"),          
          axis.title.y=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0,0.1,0,0.1), "cm"))
ggsave("plots/fig5_LClengths.pdf", device=cairo_pdf,width = 3.35, height = 7.81, units = "in")

# dev.off()
