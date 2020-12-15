#' Authors: Thomas Bryce Kelly
#' Contact: tbk14 (at) fsu.edu or tom (at) tkelly.org

source('R/MCMC Functions.R')
col = pals::alphabet(8)


{
  ## Normal
  model = list(
    obs.spreadsheet = read.xlsx('Data/NF_MCMC Data.xlsx', startRow = 2)[1:96,], ## First row is units, 97 onward is deckboard data that is not used!
    param.spreadsheet = read.xlsx('Data/NF_Model Parameters.xlsx', startRow = 1), ## First row is units
    params = init.params()
    )
  
  model$obs.spreadsheet$SYN.Biomass = model$obs.spreadsheet$SYN.Biomass * 1.2 ## Biomass adjustment based on difference between 13CPP and dilution-based PP (microscopic estimate of biomass may be low)
  model$obs.spreadsheet$PRO.Biomass = model$obs.spreadsheet$PRO.Biomass * 1.2 ## Biomass adjustment based on difference between 13CPP and dilution-based PP (microscopic estimate of biomass may be low)
  model$obs.spreadsheet$PRYM.Biomass = model$obs.spreadsheet$PRYM.Biomass * 1.2 ## Biomass adjustment based on difference between 13CPP and dilution-based PP (microscopic estimate of biomass may be low)
  model$obs.spreadsheet$DINO.Biomass = model$obs.spreadsheet$DINO.Biomass * 1.2 ## Biomass adjustment based on difference between 13CPP and dilution-based PP (microscopic estimate of biomass may be low)
  model$obs.spreadsheet$DIA.Biomass = model$obs.spreadsheet$DIA.Biomass * 1.2 ## Biomass adjustment based on difference between 13CPP and dilution-based PP (microscopic estimate of biomass may be low)
  model$obs.spreadsheet$OTHER.Biomass = model$obs.spreadsheet$OTHER.Biomass * 1.2 ## Biomass adjustment based on difference between 13CPP and dilution-based PP (microscopic estimate of biomass may be low)
}

## Setup initial model paramter values from spreadsheet
for (i in 1:length(model$params$current)) {
  l = which(model$param.spreadsheet$Param == names(model$params$current)[i])
  
  if (length(l) == 1) {
    model$params$current[[i]] = model$param.spreadsheet$Initial[l]
    model$params$Param1[[i]] = model$param.spreadsheet$Param1[l]
    model$params$Param2[[i]] = model$param.spreadsheet$Param2[l]
    model$params$Param3[[i]] = model$param.spreadsheet$Param3[l]
    model$params$step[[i]]   = model$param.spreadsheet$Step[l]
    
    if (model$param.spreadsheet$Distribution[l] == 'delta') { model$params$prior[i] = 'delta' }
    if (model$param.spreadsheet$Distribution[l] == 'lognormal') { model$params$prior[i] = 'dlnorm'}
    if (model$param.spreadsheet$Distribution[l] == 'normal') { model$params$prior[i] = 'dnorm'}
    if (model$param.spreadsheet$Distribution[l] == 'chisq') { model$params$prior[i] = 'dchisq' }
    if (model$param.spreadsheet$Distribution[l] == 'uniform') { model$params$prior[i] = 'dunif' }
    
  } else {
    warning('Parameter ', names(model$params$current)[i], ' does not match available initialization file.')
  }
}


## Initial Model Run -- start from scratch!
model$params$hist = data.frame(model$params$current)
state = init.state(); state = get.state(model)

{
  ## Random walk Routine
  # This loop will run for N iterations and perform the MCMC algorithm.
  # After the random walk, the solutions will be cut down to the last 50%
  # and then subsampled every K-th solution in order to maintain solution
  # independence. 
  a = Sys.time()
  N = 2e6
  
  for (i in 1:N) {
    model.temp = get.new.model(model, scale = 0.2) ## New model parameter set
    state.temp = get.state(model.temp) ## New sets of growths and scores
    
    if (accept.sol(state$score$Total, state.temp$score$Total)) {
      model = model.temp
      state = state.temp
    }
    if (i %% floor(N/1000) == 0) {
      b = difftime(Sys.time(), a, units = 'sec')
      message(Sys.time(), ': Completed:\t', floor(i/N*1000)/10, '%\t\tETA: ', Sys.time() + (b*1000*(1 - i/N)), '\t\tAccepted Sol: ', length(which(!is.na(state$score$HIST))), '\tScore: ', state$score$Total)
      a = Sys.time()
    }
    if (i %% floor(N/20) == 0) { save.image(file = '_rdata/model.crumb.rdata'); message(Sys.time(), ': Crumb saved.') }
  }
  message()
  message(Sys.time(), ': The ratio of accepted solutions was: ', round(1000 * length(which(!is.na(state$score$HIST))) / i) / 10, '%')
  message(Sys.time(), ': Number of solutions: ', length(which(!is.na(state$score$HIST))))
  
  
  # Remove first half of solutions and sample every K solution.
  l = length(which(!is.na(state$score$HIST)))
  BI = floor(l/2)
  K = 50
  model$params$hist = model$params$hist[c(1,seq(BI, l, by = K)),] ## Save first one since it == NEMURO
  state$score$HIST = state$score$HIST[c(1,seq(BI, l, by = K))] ## TODO See if initial parameterizations for DINO and PRYM = NEMURO-GOM
  state = get.state(model, entry = 0) ## Calculates the average uptake/growth across all solutions
  message(Sys.time(), ': Number of solutions after trimming: ', length(state$score$HIST))
}

plot(state$score$HIST, type = 'l')
plot(state$score$HIST[-1], type = 'l')


#### Plotting
## Analyze results from model and state objects

{ ## KNO3
  par(mfrow=c(1,1), plt = c(0.2, 0.8, 0.2, 0.8))
  
  plot(model$params$hist$KNO3_PRO, type = 'l', yaxt = 'n', ylab = 'KNO3', lwd = 2, ylim = c(-3.3, 2), col = col[1])
  add.log.axis(2, grid.major = T,)
  lines(model$params$hist$KNO3_SYN, col = col[2], lwd = 2)
  lines(model$params$hist$KNO3_OTHER, col = col[3], lwd = 2)
  lines(model$params$hist$KNO3_DIA, col = col[4], lwd = 2)
  lines(model$params$hist$KNO3_DINO, col = col[5], lwd = 2)
  lines(model$params$hist$KNO3_PRYM, col = col[6], lwd = 2)
  
  legend('bottomleft', legend = c('PRO', 'SYN', 'OTHER', 'DIA', 'DINO', 'PRYM'), col = col[c(1:6)], ncol = 2, lwd = 4)
  
  par(new = T, plt = c(0.805, 0.95, 0.2, 0.8))
  d = density(rnorm(1e5, model$params$Param1[names(model$params$current) == 'KNO3_PRO'], model$params$Param2[names(model$params$current) == 'KNO3_PRO']))
  d2 = density(rnorm(1e5, model$params$Param1[names(model$params$current) == 'KNO3_DIA'], model$params$Param2[names(model$params$current) == 'KNO3_DIA']))
  plot(d$y, d$x, ylim = c(-3.3,2), type = 'l', yaxt = 'n', ylab = '', xlab = '', xaxt = 'n', col = col[1], main = 'Prior')
  lines(d2$y, d2$x, col = col[4])
}

{ ## NH4
  par(mfrow=c(1,1), plt = c(0.2, 0.8, 0.2, 0.8))
  
  plot(model$params$hist$KNH4_PRO, type = 'l', yaxt = 'n', ylab = 'KNH4', lwd = 2, ylim = c(-4, 0), col = col[1]); add.log.axis(2, grid = T)
  lines(model$params$hist$KNH4_SYN, col = col[2], lwd = 2)
  lines(model$params$hist$KNH4_OTHER, col = col[3], lwd = 2)
  lines(model$params$hist$KNH4_DIA, col = col[4], lwd = 2)
  lines(model$params$hist$KNH4_DINO, col = col[5], lwd = 2)
  lines(model$params$hist$KNH4_PRYM, col = col[6], lwd = 2)
  
  legend('bottomleft', legend = c('PRO', 'SYN', 'OTHER', 'DIA', 'DINO', 'PRYM'), col = col[c(1:6)], ncol = 2, lwd = 4)
  
  par(new = T, plt = c(0.805, 0.95, 0.2, 0.8))
  d = density(rnorm(1e5, model$params$Param1[names(model$params$current) == 'KNH4_SYN'], model$params$Param2[names(model$params$current) == 'KNH4_PRO']))
  d2 = density(rnorm(1e5, model$params$Param1[names(model$params$current) == 'KNH4_DIA'], model$params$Param2[names(model$params$current) == 'KNH4_DIA']))
  plot(d$y, d$x, ylim = c(-4,0), type = 'l', yaxt = 'n', ylab = '', xlab = '', xaxt = 'n', col = col[1])
  lines(d2$y, d2$x, col = col[4])
}

{ ## V
  par(mfrow=c(1,1), plt = c(0.2, 0.8, 0.2, 0.8))
  
  plot(model$params$hist$V_PRO, type = 'l', yaxt = 'n', ylab = 'V_max', lwd = 3, ylim = c(-1, 1), col = col[1]); add.log.axis(2, grid = T)
  lines(model$params$hist$V_SYN, col = col[2], lwd = 3)
  lines(model$params$hist$V_OTHER, col = col[3], lwd = 3)
  lines(model$params$hist$V_DIA, col = col[4], lwd = 3)
  lines(model$params$hist$V_DINO, col = col[5], lwd = 3)
  lines(model$params$hist$V_PRYM, col = col[6], lwd = 3)
  
  legend('topleft', legend = c('PRO', 'SYN', 'OTHER', 'DIA', 'DINO', 'PRYM'), col = col[c(1:6)], ncol = 2, lwd = 4)
  
  par(new = T, plt = c(0.805, 0.95, 0.2, 0.8))
  d = density(rnorm(1e5, model$params$Param1[names(model$params$current) == 'V_SYN'], model$params$Param2[names(model$params$current) == 'V_PRO']))
  d2 = density(rnorm(1e5, model$params$Param1[names(model$params$current) == 'V_DIA'], model$params$Param2[names(model$params$current) == 'V_DIA']))
  plot(d$y, d$x, ylim = c(-1,1), type = 'l', yaxt = 'n', ylab = '', xlab = '', xaxt = 'n', col = col[1])
  lines(d2$y, d2$x, col = col[4])
}

{ ## R
  par(mfrow=c(1,1), plt = c(0.2, 0.8, 0.2, 0.8))
  
  plot(model$params$hist$R_PRO, type = 'l', ylab = 'R', lwd = 3, ylim = c(-3, 0), col = col[1], yaxt = 'n')
  add.log.axis(2)
  lines(model$params$hist$R_SYN, col = col[2], lwd = 3)
  lines(model$params$hist$R_OTHER, col = col[3], lwd = 3)
  lines(model$params$hist$R_DIA, col = col[4], lwd = 3)
  lines(model$params$hist$R_DINO, col = col[5], lwd = 3)
  lines(model$params$hist$R_PRYM, col = col[6], lwd = 3)
  legend('topleft', legend = c('PRO', 'SYN', 'OTHER', 'DIA', 'DINO', 'PRYM'), col = col[c(1:6)], ncol = 2, lwd = 4)
  
  par(new = T, plt = c(0.805, 0.95, 0.2, 0.8))
  d = density(rnorm(1e5, model$params$Param1[names(model$params$current) == 'R_SYN'], model$params$Param2[names(model$params$current) == 'R_PRO']))
  d2 = density(rnorm(1e5, model$params$Param1[names(model$params$current) == 'R_DIA'], model$params$Param2[names(model$params$current) == 'R_DIA']))
  plot(d$y, d$x, ylim = c(-3,0), type = 'l', yaxt = 'n', ylab = '', xlab = '', xaxt = 'n', col = col[1])
  lines(d2$y, d2$x, col = col[4])
}

#plot.limitations(state)

{
  #png('_figures/Nutrient Limitation and f-ratio Plot V2.png', width = 1200, height = 800)
  pdf('_figures/Nutrient Limitation and f-ratio Plot V2.pdf', width = 4.5, height = 7)
  par(mfrow = c(1,1), plt = c(0.2, 0.9, 0.65, 0.9))
  
  ## Nutrient Availability
  plot(NULL, NULL, xlim = c(0.9, 6.6), ylim = c(0,1.1), ylab = 'Nutrient Availability', yaxs = 'i', xaxt = 'n', xlab = '')
  abline(v = c(1:7) - 0.275, col = 'lightgrey', lty = 3)
  l = which(model$obs.spreadsheet$Depth < 50)
  w = 0.4
  col1 = '#D4E6F1'
  col2 = '#2471A3'
  add.boxplot.box(1, state$growth$SYN$NL[l]+state$growth$SYN$AL[l], width = w, col = col1, outliers = T)
  add.boxplot.box(2, state$growth$PRO$NL[l]+state$growth$PRO$AL[l], width = w, col = col1)
  add.boxplot.box(3, state$growth$OTHER$NL[l]+state$growth$OTHER$AL[l], width = w, col = col1)
  add.boxplot.box(4, state$growth$DIA$NL[l]+state$growth$DIA$AL[l], width = w, col = col1)
  add.boxplot.box(5, state$growth$DINO$NL[l]+state$growth$DINO$AL[l], width = w, col = col1)
  add.boxplot.box(6, state$growth$PRYM$NL[l]+state$growth$PRYM$AL[l], width = w, col = col1)
  
  add.boxplot.box(1.45, state$growth$SYN$NL[-l]+state$growth$SYN$AL[-l], width = w, col = col2)
  add.boxplot.box(2.45, state$growth$PRO$NL[-l]+state$growth$PRO$AL[-l], width = w, col = col2)
  add.boxplot.box(3.45, state$growth$OTHER$NL[-l]+state$growth$OTHER$AL[-l], width = w, col = col2)
  add.boxplot.box(4.45, state$growth$DIA$NL[-l]+state$growth$DIA$AL[-l], width = w, col = col2)
  add.boxplot.box(5.45, state$growth$DINO$NL[-l]+state$growth$DINO$AL[-l], width = w, col = col2)
  add.boxplot.box(6.45, state$growth$PRYM$NL[-l]+state$growth$PRYM$AL[-l], width = w, col = col2)
  #legend('topright', legend = c('Shallow (<50 m)', 'Deep (>50 m)'), col = c(col1, col2), pch = 15, cex = 0.9)
 
  ## f-ratio
  par(new = T, plt = c(0.2, 0.9, 0.4, 0.64))
  plot(NULL, NULL, xlim = c(0.9,6.6), ylim = c(0,0.5), ylab = 'f-ratio', yaxs = 'i', xlab = '', xaxt = 'n')
  abline(v = c(1:7) - 0.275, col = 'lightgrey', lty = 3)
  add.boxplot.box(1, state$growth$SYN$NL[l] / (state$growth$SYN$AL[l] + state$growth$SYN$NL[l]), width = w, col = col1)
  add.boxplot.box(2, state$growth$PRO$NL[l] / (state$growth$PRO$AL[l] + state$growth$PRO$NL[l]), width = w, col = col1)
  add.boxplot.box(3, state$growth$OTHER$NL[l] / (state$growth$OTHER$AL[l] + state$growth$OTHER$NL[l]), width = w, col = col1)
  add.boxplot.box(4, state$growth$DIA$NL[l] / (state$growth$DIA$AL[l] + state$growth$DIA$NL[l]), width = w, col = col1)
  add.boxplot.box(5, state$growth$DINO$NL[l] / (state$growth$DINO$AL[l] + state$growth$DINO$NL[l]), width = w, col = col1)
  add.boxplot.box(6, state$growth$PRYM$NL[l] / (state$growth$PRYM$AL[l] + state$growth$PRYM$NL[l]), width = w, col = col1)
  
  add.boxplot.box(1.45, state$growth$SYN$NL[-l] / (state$growth$SYN$AL[-l] + state$growth$SYN$NL[-l]), width = w, col = col2)
  add.boxplot.box(2.45, state$growth$PRO$NL[-l] / (state$growth$PRO$AL[-l] + state$growth$PRO$NL[-l]), width = w, col = col2)
  add.boxplot.box(3.45, state$growth$OTHER$NL[-l] / (state$growth$OTHER$AL[-l] + state$growth$OTHER$NL[-l]), width = w, col = col2)
  add.boxplot.box(4.45, state$growth$DIA$NL[-l] / (state$growth$DIA$AL[-l] + state$growth$DIA$NL[-l]), width = w, col = col2)
  add.boxplot.box(5.45, state$growth$DINO$NL[-l] / (state$growth$DINO$AL[-l] + state$growth$DINO$NL[-l]), width = w, col = col2)
  add.boxplot.box(6.45, state$growth$PRYM$NL[-l] / (state$growth$PRYM$AL[-l] + state$growth$PRYM$NL[-l]), width = w, col = col2)
  
  ## Light
  par(new = T, plt = c(0.2, 0.9, 0.15, 0.39))
  plot(NULL, NULL, xlim = c(0.9, 6.6), ylim = c(0,1), ylab = 'Light Availability', yaxs = 'i', xaxt = 'n', xlab = '')
  axis(1, at = c(1:7) - 0.275, labels = NA)
  axis(1, at = c(1:6) + 0.225, labels = c('SYN', 'PRO', 'OTHER', 'DIA', 'DINO', 'PRYM'), tick = F)
  abline(v = c(1:7) - 0.275, col = 'lightgrey', lty = 3)
  l = which(model$obs.spreadsheet$Depth < 50)
  add.boxplot.box(1, state$growth$SYN$LL[l], width = w, col = col1)
  add.boxplot.box(2, state$growth$PRO$LL[l], width = w, col = col1)
  add.boxplot.box(3, state$growth$OTHER$LL[l], width = w, col = col1)
  add.boxplot.box(4, state$growth$DIA$LL[l], width = w, col = col1)
  add.boxplot.box(5, state$growth$DINO$LL[l], width = w, col = col1)
  add.boxplot.box(6, state$growth$PRYM$LL[l], width = w, col = col1)
  
  add.boxplot.box(1.45, state$growth$SYN$LL[-l], width = w, col = col2)
  add.boxplot.box(2.45, state$growth$PRO$LL[-l], width = w, col = col2)
  add.boxplot.box(3.45, state$growth$OTHER$LL[-l], width = w, col = col2)
  add.boxplot.box(4.45, state$growth$DIA$LL[-l], width = w, col = col2)
  add.boxplot.box(5.45, state$growth$DINO$LL[-l], width = w, col = col2)
  add.boxplot.box(6.45, state$growth$PRYM$LL[-l], width = w, col = col2)
  dev.off()
}


{ ## Compare NEMURO answer to ours vs observations
  nemuro.state = get.state(model, 1)
  
  no3.uptake = matrix(NA, nrow = nrow(model$obs.spreadsheet), ncol = 3)
  no3.uptake[,1] = model$obs.spreadsheet$NO3.Uptake
  no3.uptake[,2] = nemuro.state$growth$TOTAL$GPP.NO3
  no3.uptake[,3] = state$growth$TOTAL$GPP.NO3
  
  nh4.uptake = matrix(NA, nrow = nrow(model$obs.spreadsheet), ncol = 3)
  nh4.uptake[,1] = model$obs.spreadsheet$NH4.Uptake
  nh4.uptake[,2] = nemuro.state$growth$TOTAL$GPP.NH4
  nh4.uptake[,3] = state$growth$TOTAL$GPP.NH4
  
  NPP = matrix(NA, nrow = nrow(model$obs.spreadsheet), ncol = 3)
  NPP[,1] = model$obs.spreadsheet$NPP.N
  NPP[,2] = nemuro.state$growth$TOTAL$NPP
  NPP[,3] = state$growth$TOTAL$NPP
}


{
  par(plt = c(0.2, 0.95, 0.1, 0.8), mfrow = c(2,1))
  boxplot(cbind(nemuro.state$score$NO3, state$score$NO3,
                nemuro.state$score$NH4, state$score$NH4,
                nemuro.state$score$NPP, state$score$NPP), names = c('NEM', 'NO3', 'NEM', 'NH4', 'NEM', 'NPP'),
          main = 'Costs', ylim = c(0,300), col = c('#ff000090', 'light grey'), yaxs = 'i', outline = T)
  grid(); box()
  
  par(plt = c(0.2, 0.95, 0.4, 0.9))
  barplot(apply(cbind(nemuro.state$score$NO3, state$score$NO3,
                nemuro.state$score$NH4, state$score$NH4,
                nemuro.state$score$NPP, state$score$NPP), 2, function(x) {sum(x, na.rm = T)}),
          col = c('#ff000090', 'light grey'), ylim = c(0,20000), names = c('NEM', 'NO3', 'NEM', 'NH4', 'NEM', 'NPP'),)
  grid(); box()
}

pdf('_figures/NPP and uptake.pdf', width = 7, height = 4)
{
  
  plot(model$obs.spreadsheet$NO3.Uptake, state$growth$TOTAL$GPP.NO3, pch = 16, xlim = c(0,0.08), ylim = c(0,0.08),
       ylab = 'Model NO3 Uptake', xlab = 'Observed NO3 Uptake', xaxs = 'i', yaxs = 'i', col = '#00000090',
       cex = make.cex(model$obs.spreadsheet$NO3.Uptake / model$obs.spreadsheet$NO3.Uptake.Sigma, max = 2, min = 0.3),
       main = paste('NO3 Uptake ~~ n =', length(which(!is.na(model$obs.spreadsheet$NO3.Uptake)))))
  abline(a=0, b=1, lty = 2)
  grid(); box()
  points(model$obs.spreadsheet$NO3.Uptake, nemuro.state$growth$TOTAL$GPP.NO3, pch = 16, col = '#ff000080',
         cex = make.cex(model$obs.spreadsheet$NO3.Uptake / model$obs.spreadsheet$NO3.Uptake.Sigma, max = 2, min = 0.3))
  
  #mtext(side = 3, line = -1, adj = 0.05, paste('Cost:', round(-0.5*sum(state$score$NO3, na.rm = T), digits = 2)))
  #mtext(side = 3, line = -1, adj = 0.3, paste('Cost:', round(-0.5*sum(nemuro.state$score$NO3, na.rm = T),digits = 2)), col = 'red')
}


{ ## Compare model-obs NPP (tuned and untuned)
  par(mfrow = c(1,2))
  plot(model$obs.spreadsheet$NPP.N, state$growth$TOTAL$NPP, ylab = 'Model NPP', xlab = 'Observed NPP',
       pch = 16, col = '#00000090', xlim = c(-0.05, 0.15), ylim = c(-0.05,0.15),
       cex = make.cex(model$obs.spreadsheet$NPP.N / model$obs.spreadsheet$NPP.N.Sigma, max = 2, min = 0.4),
       main = paste('NPP ~~ n =', length(which(!is.na(model$obs.spreadsheet$NPP.N)))))
  grid(); box()
  points(model$obs.spreadsheet$NPP.N, nemuro.state$growth$TOTAL$NPP, pch = 16, col = '#ff000080',
         cex = make.cex(model$obs.spreadsheet$NPP.N / model$obs.spreadsheet$NPP.N.Sigma, max = 2, min = 0.4))
  abline(a=0, b=1, lty = 2)
  
  #mtext(side = 3, line = -1, adj = 0.05, paste('Cost:', round(-0.5*sum(state$score$NPP, na.rm = T))))
  #mtext(side = 3, line = -1, adj = 0.25, paste('Cost:', round(-0.5*sum(nemuro.state$score$NPP, na.rm = T))), col = 'red')
}
dev.off()

pdf('_figures/Growth Rates.pdf', width = 8, height = 5.75)
{ ## Compare model-obs NPP (tuned and untuned)
  par(mfrow = c(2,3))
  
  ## PRO
  plot(model$obs.spreadsheet$PRO.Growth.Rate, state$growth$PRO$Rate, ylab = 'Model mu (d-1)', xlab = 'Observed mu (d-1)',
       main = paste('PRO ~~ n =', length(which(!is.na(model$obs.spreadsheet$PRO.Growth)))),
       pch = 16, col = '#00000070', xlim = c(-1.0, 2), ylim = c(-1.0,2))
  grid()
  points(model$obs.spreadsheet$PRO.Growth.Rate, nemuro.state$growth$PRO$Rate, pch = 16, col = '#ff000070')
  abline(a=0, b=1, lty = 2)
  #mtext(side = 3, line = -1.5, adj = 0.05, paste('Cost:', round(-0.5*sum(state$score$Mu.PRO, na.rm = T), digits = 1)), cex = 0.7)
  #mtext(side = 3, line = -1.5, adj = 0.5, paste('Cost:', round(-0.5*sum(nemuro.state$score$Mu.PRO, na.rm = T), digits = 1)), col = 'red', cex = 0.7)
  
  ## SYN
  plot(model$obs.spreadsheet$SYN.Growth.Rate, state$growth$SYN$Rate, ylab = 'Model mu (d-1)', xlab = 'Observed mu (d-1)',
       main = paste('SYN ~~ n =', length(which(!is.na(model$obs.spreadsheet$SYN.Growth)))),
       pch = 16, col = '#00000070', xlim = c(-1.0, 2), ylim = c(-1.0,2))
  grid()
  points(model$obs.spreadsheet$SYN.Growth.Rate, nemuro.state$growth$SYN$Rate, pch = 16, col = '#ff000070')
  abline(a=0, b=1, lty = 2)
  #mtext(side = 3, line = -1.5, adj = 0.05, paste('Cost:', round(-0.5*sum(state$score$Mu.SYN, na.rm = T), digits = 1)), cex = 0.7)
  #mtext(side = 3, line = -1.5, adj = 0.5, paste('Cost:', round(-0.5*sum(nemuro.state$score$Mu.SYN, na.rm = T), digits = 1)), col = 'red', cex = 0.7)
  
  ## OTHER
  plot(model$obs.spreadsheet$OTHER.Growth.Rate, state$growth$OTHER$Rate, ylab = 'Model mu (d-1)', xlab = 'Observed mu (d-1)',
       main = paste('OTHER ~~ n =', length(which(!is.na(model$obs.spreadsheet$OTHER.Growth)))),
       pch = 16, col = '#00000070', xlim = c(-1.0, 2), ylim = c(-1.0,2))
  grid()
  points(model$obs.spreadsheet$OTHER.Growth.Rate, nemuro.state$growth$OTHER$Rate, pch = 16, col = '#ff000070')
  abline(a=0, b=1, lty = 2)
  #mtext(side = 3, line = -1.5, adj = 0.05, paste('Cost:', round(-0.5*sum(state$score$Mu.OTHER, na.rm = T), digits = 1)), cex = 0.7)
  #mtext(side = 3, line = -1.5, adj = 0.5, paste('Cost:', round(-0.5*sum(nemuro.state$score$Mu.OTHER, na.rm = T), digits = 1)), col = 'red', cex = 0.7)
  
  ## DIA
  plot(model$obs.spreadsheet$DIA.Growth.Rate, state$growth$DIA$Rate, ylab = 'Model mu (d-1)', xlab = 'Observed mu (d-1)',
       main = paste('DIA ~~ n =', length(which(!is.na(model$obs.spreadsheet$DIA.Growth)))),
       pch = 16, col = '#00000070', xlim = c(-1.0, 2), ylim = c(-1.0,2))
  grid()
  points(model$obs.spreadsheet$DIA.Growth.Rate, nemuro.state$growth$DIA$Rate, pch = 16, col = '#ff000070')
  abline(a=0, b=1, lty = 2)
  #mtext(side = 3, line = -1.5, adj = 0.05, paste('Cost:', round(-0.5*sum(state$score$Mu.DIA, na.rm = T), digits = 1)), cex = 0.7)
  #mtext(side = 3, line = -1.5, adj = 0.5, paste('Cost:', round(-0.5*sum(nemuro.state$score$Mu.DIA, na.rm = T), digits = 1)), col = 'red', cex = 0.7)
  
  ## DINO
  plot(model$obs.spreadsheet$DINO.Growth.Rate, state$growth$DINO$Rate, ylab = 'Model mu (d-1)', xlab = 'Observed mu (d-1)',
       main = paste('DINO ~~ n =', length(which(!is.na(model$obs.spreadsheet$DINO.Growth)))),
       pch = 16, col = '#00000070', xlim = c(-1.0, 2), ylim = c(-1.0,2))
  grid()
  points(model$obs.spreadsheet$DINO.Growth.Rate, nemuro.state$growth$DINO$Rate, pch = 16, col = '#ff000070')
  abline(a=0, b=1, lty = 2)
  #mtext(side = 3, line = -1.5, adj = 0.05, paste('Cost:', round(-0.5*sum(state$score$Mu.DINO, na.rm = T), digits = 1)), cex = 0.7)
  #mtext(side = 3, line = -1.5, adj = 0.5, paste('Cost:', round(-0.5*sum(nemuro.state$score$Mu.DINO, na.rm = T), digits = 1)), col = 'red', cex = 0.7)
  
  ## PRYM
  plot(model$obs.spreadsheet$PRYM.Growth.Rate, state$growth$PRYM$Rate, ylab = 'Model mu (d-1)', xlab = 'Observed mu (d-1)',
       main = paste('PRYM ~~ n =', length(which(!is.na(model$obs.spreadsheet$PRYM.Growth)))),
       pch = 16, col = '#00000070', xlim = c(-1.0, 2), ylim = c(-1.0,2))
  grid()
  points(model$obs.spreadsheet$PRYM.Growth.Rate, nemuro.state$growth$PRYM$Rate, pch = 16, col = '#ff000070')
  abline(a=0, b=1, lty = 2)
  #mtext(side = 3, line = -1.5, adj = 0.05, paste('Cost:', round(-0.5*sum(state$score$Mu.PRYM, na.rm = T), digits = 1)), cex = 0.7)
  #mtext(side = 3, line = -1.5, adj = 0.5, paste('Cost:', round(-0.5*sum(nemuro.state$score$Mu.PRYM, na.rm = T), digits = 1)), col = 'red', cex = 0.7)
}
dev.off()

{  ## Big correlation Figure
  pdf(file = '_figures/Model Parameter Correlations.pdf')
  plot.matrix(c(2:6)) ## SYN
  plot.matrix(c(9:13)) ## PRO
  plot.matrix(c(16:20)) ## OTHER
  plot.matrix(c(23:27)) ## DIA
  plot.matrix(c(30:34)) ## DINO
  plot.matrix(c(37:41)) ## PRYM
  dev.off()
}

pdf('_figures/Violin Plot of Priors.pdf')
{
  par(plt = c(0.2, 0.5, 0.2, 0.9))
  plot(NULL, NULL, ylim = c(6.5,0.5), xlim = c(-2.5,1), xaxt = 'n', xlab = 'KNH4 (uM)', ylab = '', yaxt = 'n')
  axis(2, at = c(1:6), labels = c('SYN', 'PRO', 'OTHER', 'DIA', 'DINO', 'PRYM'), las = 1)
  add.log.axis(side = 1, grid.major = T)
  
  add.violin(1, rnorm(10*nrow(model$params$hist), model$params$Param1[3], model$params$Param2[3]), col[1], 'grey', border = T, scale = 0.3)
  add.violin(2, rnorm(10*nrow(model$params$hist), model$params$Param1[10], model$params$Param2[10]), col[2], 'grey', border = T, scale = 0.3)
  add.violin(3, rnorm(10*nrow(model$params$hist), model$params$Param1[17], model$params$Param2[17]), col[3], 'grey', border = T, scale = 0.3)
  add.violin(4, rnorm(10*nrow(model$params$hist), model$params$Param1[24], model$params$Param2[24]), col[4], 'grey', border = T, scale = 0.3)
  add.violin(5, rnorm(10*nrow(model$params$hist), model$params$Param1[31], model$params$Param2[31]), col[5], 'grey', border = T, scale = 0.3)
  add.violin(6, rnorm(10*nrow(model$params$hist), model$params$Param1[38], model$params$Param2[38]), col[6], 'grey', border = T, scale = 0.3)
  
  add.violin(1, model$params$hist$KNH4_SYN[-1], col[1], 'grey', scale = 0.3)
  add.violin(2, model$params$hist$KNH4_PRO[-1], col[2], 'grey', scale = 0.3)
  add.violin(3, model$params$hist$KNH4_OTHER[-1], col[3], 'grey', scale = 0.3)
  add.violin(4, model$params$hist$KNH4_DIA[-1], col[4], 'grey', scale = 0.3)
  add.violin(5, model$params$hist$KNH4_DINO[-1], col[5], 'grey', scale = 0.3)
  add.violin(6, model$params$hist$KNH4_PRYM[-1], col[6], 'grey', scale = 0.3)
  
  par(plt = c(0.51, 0.8, 0.2, 0.9), new = T)
  plot(NULL, NULL, ylim = c(6.5,0.5), xlim = c(-2.5,1), xaxt = 'n', yaxt = 'n', ylab = '', xlab = 'KNO3 (uM)')
  add.log.axis(side = 1, grid.major = T)
  
  add.violin(1, rnorm(10*nrow(model$params$hist), model$params$Param1[2], model$params$Param2[2]), col[1], 'grey', border = T, scale = 0.3)
  add.violin(2, rnorm(10*nrow(model$params$hist), model$params$Param1[9], model$params$Param2[9]), col[2], 'grey', border = T, scale = 0.3)
  add.violin(3, rnorm(10*nrow(model$params$hist), model$params$Param1[16], model$params$Param2[16]), col[3], 'grey', border = T, scale = 0.3)
  add.violin(4, rnorm(10*nrow(model$params$hist), model$params$Param1[23], model$params$Param2[23]), col[4], 'grey', border = T, scale = 0.3)
  add.violin(5, rnorm(10*nrow(model$params$hist), model$params$Param1[30], model$params$Param2[30]), col[5], 'grey', border = T, scale = 0.3)
  add.violin(6, rnorm(10*nrow(model$params$hist), model$params$Param1[37], model$params$Param2[37]), col[6], 'grey', border = T, scale = 0.3)
  
  add.violin(1, model$params$hist$KNO3_SYN[-1], col[1], 'grey', scale = 0.3)
  add.violin(2, model$params$hist$KNO3_PRO[-1], col[2], 'grey', scale = 0.3)
  add.violin(3, model$params$hist$KNO3_OTHER[-1], col[3], 'grey', scale = 0.3)
  add.violin(4, model$params$hist$KNO3_DIA[-1], col[4], 'grey', scale = 0.3)
  add.violin(5, model$params$hist$KNO3_DINO[-1], col[5], 'grey', scale = 0.3)
  add.violin(6, model$params$hist$KNO3_PRYM[-1], col[6], 'grey', scale = 0.3)
  
  add.violin(1, rnorm(nrow(model$params$hist), model$params$Param1[2], model$params$Param2[2]), col[1], 'grey', border = T, scale = 0.3)
}
dev.off()




#### Scores
lapply(state$score, function(x) {sum(x, na.rm = T)})
