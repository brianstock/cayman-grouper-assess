# Fig 4: pre-protection length by island

dir <- here()
setwd(dir)
plots_dir <- file.path(dir, "plots")

library(tidyverse)
library(ggthemes)
library(RColorBrewer)
library(viridis)
library(ggridges)

dat <- read.csv("data/NG_pre2001_biological_data_2002.csv", header=T)
islands <- c("Little Cayman 1987-1995", "Little Cayman 2002", "Cayman Brac", "Grand Cayman", "Pickle Bank")
dat$Island <- as.character(dat$Island)
dat$Island[dat$Island=="Grand_12MileBanks"] <- "Grand"
dat$Island[dat$Island=="PickleBank"] <- "Pickle Bank"
dat$Island[dat$Island=="Grand"] <- "Grand Cayman"
# dat$Island[dat$Island=="Little"] <- "Little Cayman"
dat$Island[dat$Island=="Little" & dat$Year < 2001] <- "Little Cayman 1987-1995"
dat$Island[dat$Island=="Little" & dat$Year > 2001] <- "Little Cayman 2002"
dat$Island[dat$Island=="Brac"] <- "Cayman Brac"
dat$Island <- factor(dat$Island, levels=islands)

# convert lengths in inches to cm
dat$Length_cm <- dat$Length
dat$Length_cm[dat$Length_units == "in"] <- dat$Length[dat$Length_units == "in"] * 2.54

# 0-2 yr-old fish were used in growth curve but should be removed for analysis of 
# length distributions of spawning populations (not caught on SPAGs, summer instead of winter)
table(dat$Age)
dat <- filter(dat, Age > 1 | is.na(Age))

# sample size by island
sampsizes <- dat %>% group_by(Island) %>% summarize(n=n()) %>% 
                     mutate(nlab=paste0("n = ",as.character(n))) %>% as.data.frame
sampsizes$Length = 80 # x-axis position
# sampsizes$y = 0.8
sampsizes$y = c(0.065,0.049,0.052,0.074,0.071)

# -------------------------------------------------------------------------
# Histograms by island, color to match SPR plot
# first calculate mean by island (to add with vert line)
island.means <- by(dat$Length_cm, dat$Island, mean, na.rm=T)
df <- data.frame(islandmean = unclass(island.means), Island = names(island.means), y=0.3)
df$islandmean <- round(df$islandmean, 1)
# island.labs <- data.frame(x = 39, y = 0.8, Island = names(island.means))
dat <- dat %>% group_by(Island) %>% mutate(islandmean = mean(Length_cm,na.rm=T)) 
mode.tmp <- data.frame(Length_cm=as.numeric(df$islandmean),
                      Island=islands,
                      yend=c(0.06,0.033,0.05,0.074,0.07),
                      y=rep(0,5))
island.labs <- data.frame(x = 39, y = c(0.06,0.044,0.052,0.074,0.071), Island=islands, lab = c("Little Cayman\n1987-1995", "Little Cayman\n2002", "Cayman Brac", "Grand Cayman", "Pickle Bank"))

dev.new(width = 3.34645669291339, height = 3.90419947506562)
ggplot(data = dat, aes(x=Length_cm)) +
  geom_histogram(color = "grey20", bins=22, size=0.2,alpha=0.6,mapping = aes(fill = Island, y = ..density..)) +
  geom_density(color='grey20', size=.5, alpha = 1) +
  scale_fill_manual(values=c(brewer.pal(4,"Set2")[1], brewer.pal(4,"Set2")), guide=FALSE) +
  scale_color_manual(values=c(brewer.pal(4,"Set2")[1], brewer.pal(4,"Set2")), guide=FALSE) +
  geom_segment(data = mode.tmp, aes(x=Length_cm, xend=Length_cm, y=y, yend=yend), color='grey20', linetype=2, alpha=1,lwd=0.8, inherit.aes = F) +
  geom_text(data=sampsizes, aes(x=Length, y=y, label=nlab), size = 2.8,inherit.aes=FALSE) +
  geom_text(data=island.labs, aes(x=x ,y=y, label=lab), lineheight=.8, hjust = 0, size=3) +
  facet_grid(Island ~ ., scales="free_y") +
  coord_cartesian(xlim=c(41, 83)) +
  scale_y_continuous(expand = c(0.01, 0.001)) +
  theme_pander() +
  theme(axis.text.x= element_text(size = 9),
          axis.title.x= element_text(size = 10),
          axis.title.y= element_text(size = 10),
          axis.text.y=element_blank(),
          axis.ticks.x=element_line(colour="black"),
          axis.ticks.y = element_blank(),
          strip.background = element_blank(),
          strip.text.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
  labs(x = "Total length (cm)", y = "Proportional density")
ggsave(file.path(plots_dir, "fig4_growth.pdf"), device=cairo_pdf,width = 85, height = 99.167, units = "mm")
dev.off()

