# Figure 9
# Compare recruitment from 4 LIME models with differing amount of uncertainty in length/growth
# More weight given to the length data = less uncertainty in 2011 recruitment spike:
#   Francis weighting (vs. Dirichlet-multinomial)
#   Using growth parameter values fit outside assessment (LIME-fixed-growth) vs. estimating them in the assessment (LIME-integrated)
# Must first run the 4 models:
#   fig8a_run_LIME_integrated_francis.R
#   fig9b_run_LIME_fixed_growth_francis.R
#   fig9c_run_LIME_integrated_dirichlet.R
#   fig9d_run_LIME_fixed_growth_dirichlet.R

dir <- here()
setwd(dir)

library(LIME)

plot_rec_4panel <- function(mods, ylim=NULL){
  Inputs=mods[[1]]$Inputs
  Report=mods[[1]]$Report
  Sdreport=mods[[1]]$Sdreport

  LBSPR=NULL
  true_years=mods[[1]]$years
  True=NULL
  legend=FALSE
  relative=FALSE
  est_rdev <- rep(1,length(true_years))
  est_rdev[which(is.na(Inputs$Map$Nu_input))] = 0

  xlim.offset = 0.5 # add buffer to x-axis (year)
  nf <- Inputs$Data$n_f
  ns <- Inputs$Data$n_s

  all_years <- 1:Inputs$Data$n_t
  lc_years <- lapply(1:nf, function(x){
    sub <- Inputs$Data$LF_tlf[,,x]
    if(is.vector(sub) & sum(sub)>0) out <- all_years
    if(is.matrix(sub)) out <- which(rowSums(sub) > 0)
    return(out)
  })
  if(all(is.null(true_years))) true_years <- all_years

  if(length(mods)==1) dim <- c(1,1)
  if(length(mods)==2) dim <- c(2,1)
  if(length(mods)==3) dim <- c(3,1)
  if(length(mods)==4) dim <- c(2,2)
  if(length(mods)==5 | length(mods)==6) dim <- c(2,3)
  by <- round(length(true_years)/5)
  lab <- rev(seq(from=true_years[length(true_years)], to=min(true_years), by=-by))
  ilab <- which(true_years %in% lab)

  if(all(is.null(Inputs))==FALSE){
    if(ns==1){
      xY <- seq_along(all_years)
      xLC <- lapply(1:nf, function(x) which(all_years %in% lc_years[[x]]))
    }
    if(ns>1){
      xY <- 1:Inputs$Data$n_y
      xLC <- lapply(1:nf, function(x) unique(Inputs$Data$S_yrs[which(all_years %in% lc_years[[x]])]))
    }
    ilab2 <- sapply(1:length(ilab), function(x){
      sub <- which(Inputs$Data$S_yrs %in% ilab[x])
      return(sub[length(sub)])
    })
  }

  par(mfrow=dim, mgp=c(1.8,0.5,0), oma = c(0,0,0,0) + 0.1)
  col_total <- "#228B22"
  cols <- col_total

  for(i in 1:length(mods)){
    Inputs=mods[[i]]$Inputs
    Report=mods[[i]]$Report
    Sdreport=mods[[i]]$Sdreport

    if(all(is.na(Sdreport))==FALSE){
        sd <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lR_t"),]
        sd[,2][which(is.na(sd[,2]))] <- 0
        sd <- sd[seq(1,by=ns,length.out=Inputs$Data$n_y),]
        r_est <- exp(sd[,1])/10000
        doplot <- as.logical(est_rdev)
        if(is.null(ylim)) ylim <- c(0, max(read_sdreport(sd, log=TRUE)))
        
        if(i < dim[2]){
          par(mar = c(1.5,3.5,0.5,.5) + 0.1) # top row, no x-axis label
          plot(x=1, y=1, type="n", xaxt="n", xaxs="i", yaxs="i", ylim=ylim, cex.axis=1.1, cex.lab=1.3, ylab=expression("Recruitment (x"*10^4*")"), xlab="", xlim=c(min(xY)-xlim.offset, max(xY)+xlim.offset))          
        }
        if(i == dim[2]){
          par(mar = c(1.5,1.5,0.5,1) + 0.1) # top row, no x-axis label
          plot(x=1, y=1, type="n", xaxt="n", xaxs="i", yaxs="i", ylim=ylim, cex.axis=1.1, cex.lab=1.3, ylab="", xlab="", xlim=c(min(xY)-xlim.offset, max(xY)+xlim.offset))          
        }
        if(i > dim[2] & i < length(mods)){
          par(mar = c(3,3.5,0,.5) + 0.1) # bottom row, space for x-axis label
          plot(x=1, y=1, type="n", xaxt="n", xaxs="i", yaxs="i", ylim=ylim, cex.axis=1.1, cex.lab=1.3, ylab=expression("Recruitment (x"*10^4*")"), xlab="Year", xlim=c(min(xY)-xlim.offset, max(xY)+xlim.offset))
        }
        if(i == length(mods)){
          par(mar = c(3,1.5,0,1) + 0.1) # bottom row, space for x-axis label
          plot(x=1, y=1, type="n", xaxt="n", xaxs="i", yaxs="i", ylim=ylim, cex.axis=1.1, cex.lab=1.3, ylab="", xlab="Year", xlim=c(min(xY)-xlim.offset, max(xY)+xlim.offset))
        }
        polygon( y=read_sdreport(sd, log=TRUE)/10000, x=c(which(is.na(sd[,2])==FALSE), rev(which(is.na(sd[,2])==FALSE))), col=paste0(col_total, "40"), border=NA)
        lines(x=seq_along(xY), y=r_est, lwd=2, col=col_total, ylim=ylim, xpd=NA)
        points(x=all_years[doplot], y=r_est[doplot], col=col_total, pch=19, cex=1.5, xpd=NA)
        axis(1, cex.axis=1.1, cex.lab=1.3, at=ilab2, labels=lab)
        if(all(is.null(True))==FALSE) lines(True$R_t[seq(1,by=ns,length.out=Inputs$Data$n_y)], lwd=2)
    }
    mtext(LETTERS[i], side = 3, line = -2.2, adj = 0.03, cex = 2) # add letters to plots
  } # end plot loop
}

mods <- list(readRDS("results/LIME_integrated_francis/output/output_CV0.096_sigI0.175_sigC0.2_francis_2.rds"),
            readRDS("results/LIME_fixed_growth_francis/output/output_CV0.096_francis_2.rds"),
            readRDS("results/LIME_integrated_dirichlet/output/output_CV0.096_integrated.rds"),
            readRDS("results/LIME_fixed_growth_dirichlet/output/output_CV0.096.rds"))
cairo_pdf(filename="plots/fig9_rec_4panel.pdf", width=8, height=6)
plot_rec_4panel(mods, ylim=c(0,4.9))
dev.off()

