plot_stdresiduals <- function(mod){
  if(!mod$opt$na_sdrep){
    # get index residuals
    years <- mod$years
    ny = length(years)
    ni = dim(mod$Inputs$Data$I_ft)[1]
    temp = summary(mod$Sdreport)
    ind = rownames(temp) == "log_index_resid"
    ind2 = which(!is.nan(temp[ind,1]))
    temp = temp[ind,1][ind2]/temp[ind,2][ind2]
    # ind = which(!is.nan(temp[,1]))
    xi = data.frame(Label = rep("Index",length(ind2)),
        Year = years[ind2],
        Stdres = temp)

    # get mean length residuals
    temp = summary(mod$Sdreport)
    ind = rownames(temp) == "ML_resid"
    temp = matrix(temp[ind,1]/temp[ind,2], ny, 1)
    ind = which(temp[,1] != 0)
    xl = data.frame(Label = rep("Mean length",length(ind)),
      Year = years[ind],
      Stdres = temp[ind,1])

    x <- rbind(xi, xl)
    x$Label = factor(x$Label)
    ggp = ggplot(x, aes(x=Year, y = Stdres)) +
      # geom_ribbon(aes(ymin=lo, ymax=hi, fill=Label), alpha=0.3, linetype = 0) +
      # geom_line(size=1.1) +
      # geom_ribbon(aes(ymin=-1.96, ymax=1.96), alpha=0.2, fill='grey', color=NA) +
      geom_hline(yintercept = 0, linetype=2) +
      # geom_pointrange(aes(ymin=lo, ymax=hi)) +
      geom_point(size=1.1) +
      ylab("Standardized residual") +
      ylim(-4,4) +
      xlim(2000,2019) +
      theme_bw() +
      theme(legend.position = "none") +
      facet_wrap(vars(Label))  
  } else { # plot raw residuals
    # get index residuals
    years <- mod$years
    ny = length(years)
    ni = dim(mod$Inputs$Data$I_ft)[1]
    temp = summary(mod$Sdreport)
    ind = rownames(temp) == "log_index_resid"
    temp = matrix(temp[ind,1], ny, 1)
    ind = which(!is.nan(temp[,1]))
    xi = data.frame(Label = rep("Index",length(ind)),
        Year = years[ind],
        Stdres = temp[ind,1])

    # get mean length residuals
    temp = summary(mod$Sdreport)
    ind = rownames(temp) == "ML_resid"
    temp = matrix(temp[ind,1], ny, 1)
    ind = which(temp[,1] != 0)
    xl = data.frame(Label = rep("Mean length",length(ind)),
      Year = years[ind],
      Stdres = temp[ind,1])

    x <- rbind(xi, xl)
    x$Label = factor(x$Label)
    ggp = ggplot(x, aes(x=Year, y = Stdres, color=Label)) +
      # geom_ribbon(aes(ymin=lo, ymax=hi, fill=Label), alpha=0.3, linetype = 0) +
      # geom_line(size=1.1) +
      geom_ribbon(aes(ymin=-1.96, ymax=1.96), alpha=0.2, fill='grey', color=NA) +
      geom_hline(yintercept = 0, linetype=2) +
      # geom_pointrange(aes(ymin=lo, ymax=hi)) +
      geom_point(size=1.1) +
      ylab("Raw residual (no SE estimate)") +
      ylim(-4,4) +
      theme_bw() +
      theme(legend.position = "none") +
      facet_wrap(vars(Label))     
  }

  # # get index residuals
  # years <- mod$years
  # ny = length(years)
  # ni = dim(mod$Inputs$Data$I_ft)[1]
  # temp = summary(mod$Sdreport)
  # ind = rownames(temp) == "log_index_resid"
  # templo = matrix(temp[ind,1] - qnorm(0.975)*temp[ind,2], ny, ni)
  # temphi = matrix(temp[ind,1] + qnorm(0.975)*temp[ind,2], ny, ni)
  # temp = matrix(temp[ind,1], ny, ni)
  # xi = data.frame(Label = integer(),
  #   Year = numeric(),
  #   Stdres = numeric(),
  #   lo = numeric(),
  #   hi = numeric())
  # for(i in 1:ni)
  # {
  #   ind = which(!is.nan(temp[,i]))
  #   td = data.frame(Label = rep("Index",length(ind)),
  #     Year = years[ind],
  #     Stdres = temp[ind,i],
  #     lo = templo[ind,i],
  #     hi = temphi[ind,i])
  #   xi <- rbind(xi, td)
  # }

  # # get mean length residuals
  # temp = summary(mod$Sdreport)
  # ind = rownames(temp) == "ML_resid"
  # templo = matrix(temp[ind,1] - qnorm(0.975)*temp[ind,2], ny, 1)
  # temphi = matrix(temp[ind,1] + qnorm(0.975)*temp[ind,2], ny, 1)
  # temp = matrix(temp[ind,1], ny, 1)
  # ind = which(temp[,1] != 0)
  # xl = data.frame(Label = rep("Mean length",length(ind)),
  #   Year = years[ind],
  #   Stdres = temp[ind,1],
  #   lo = templo[ind,1],
  #   hi = temphi[ind,1])  


  return(list(ggp,x))
}
