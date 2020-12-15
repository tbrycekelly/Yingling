library(TheSource) ## devtools::install_github('tbrycekelly/thesource')
library(openxlsx)
library(pals)
library(data.table)

max.iter = 1e6

delta = function(x, x0) {
  as.numeric(x == x0)
}

#' @title Initialize Parameters
#' @author Thomas Bryce Kelly
#' @export
init.params = function() {
  params = data.table(
    Q         = NA, # Temperature Limitation Factor
    
    KNO3_SYN  = NA, # Half Saturation (Nitrate) [umol N m-3]
    KNH4_SYN  = NA, # Half Saturation (Ammonium) [umol N m-3]
    Alpha_SYN = NA, # Alpha   [m2 W-1 d-1]
    Beta_SYN  = NA, # Beta    [m2 W-1 d-1]
    V_SYN     = NA, # Maximum Growth rate  [d-1]
    R_SYN     = NA, # Respiration Rate [d-1]
    E_SYN     = NA, # Excretion Factor [unitless]
    
    KNO3_PRO  = NA, # Half Saturation (Nitrate) [umol N m-3]
    KNH4_PRO  = NA, # Half Saturation (Ammonium) [umol N m-3]
    Alpha_PRO = NA, # Alpha   [m2 W-1 d-1]
    Beta_PRO  = NA, # Beta    [m2 W-1 d-1]
    V_PRO     = NA, # Maximum Growth rate  [d-1]
    R_PRO     = NA, # Respiration Rate  [d-1]
    E_PRO     = NA, # Excretion Factor [unitless]
    
    KNO3_OTHER  = NA, # Half Saturation (Nitrate) [umol N m-3]
    KNH4_OTHER  = NA, # Half Saturation (Ammonium) [umol N m-3]
    Alpha_OTHER = NA, # Alpha   [m2 W-1 d-1]
    Beta_OTHER  = NA, # Beta    [m2 W-1 d-1]
    V_OTHER     = NA, # Maximum Growth rate  [d-1]
    R_OTHER     = NA, # Respiration Rate  [d-1]
    E_OTHER     = NA, # Excretion Factor [unitless]
    
    KNO3_DIA  = NA, # Half Saturation (Nitrate) [umol N m-3]
    KNH4_DIA  = NA, # Half Saturation (Ammonium) [umol N m-3]
    Alpha_DIA = NA, # Alpha   [m2 W-1 d-1]
    Beta_DIA  = NA, # Beta    [m2 W-1 d-1]
    V_DIA     = NA, # Maximum Growth rate  [d-1]
    R_DIA     = NA, # Respiration Rate  [d-1]
    E_DIA     = NA, # Excretion Factor [unitless]
    
    KNO3_DINO  = NA, # Half Saturation (Nitrate) [umol N m-3]
    KNH4_DINO  = NA, # Half Saturation (Ammonium) [umol N m-3]
    Alpha_DINO = NA, # Alpha   [m2 W-1 d-1]
    Beta_DINO  = NA, # Beta    [m2 W-1 d-1]
    V_DINO     = NA, # Maximum Growth rate  [d-1]
    R_DINO     = NA, # Respiration Rate  [d-1]
    E_DINO     = NA, # Excretion Factor [unitless]
    
    KNO3_PRYM  = NA, # Half Saturation (Nitrate) [umol N m-3]
    KNH4_PRYM  = NA, # Half Saturation (Ammonium) [umol N m-3]
    Alpha_PRYM = NA, # Alpha   [m2 W-1 d-1]
    Beta_PRYM  = NA, # Beta    [m2 W-1 d-1]
    V_PRYM     = NA, # Maximum Growth rate  [d-1]
    R_PRYM     = NA, # Respiration Rate  [d-1]
    E_PRYM     = NA # Excretion Factor [unitless]
  )
  
  hist = matrix(params, ncol = length(params), nrow = max.iter)
  colnames(hist) = names(params)
  
  list(current = params,
       prior = rep(NA, length(params)),
       Param1 = rep(NA, length(params)),
       Param2 = rep(NA, length(params)),
       Param3 = rep(NA, length(params)),
       step = rep(NA, length(params)),
       hist = hist)
}


#' @title Initialize State
#' @author Thomas Bryce Kelly
#' 
init.state = function() {
  list(
    growth = list(
      SYN = data.table(TL = 0, NL = 0, AL = 0, LL = 0, GPP.NO3 = 0, GPP.NH4 = 0, NPP = 0, Rate = 0),
      PRO = data.table(TL = 0, NL = 0, AL = 0, LL = 0, GPP.NO3 = 0, GPP.NH4 = 0, NPP = 0, Rate = 0),
      OTHER = data.table(TL = 0, NL = 0, AL = 0, LL = 0, GPP.NO3 = 0, GPP.NH4 = 0, NPP = 0, Rate = 0),
      DIA = data.table(TL = 0, NL = 0, AL = 0, LL = 0, GPP.NO3 = 0, GPP.NH4 = 0, NPP = 0, Rate = 0),
      DINO = data.table(TL = 0, NL = 0, AL = 0, LL = 0, GPP.NO3 = 0, GPP.NH4 = 0, NPP = 0, Rate = 0),
      PRYM = data.table(TL = 0, NL = 0, AL = 0, LL = 0, GPP.NO3 = 0, GPP.NH4 = 0, NPP = 0, Rate = 0),
      TOTAL = data.table(GPP.NO3 = 0, GPP.NH4 = 0, NPP = 0)
    ),
    score = list(
      NO3 = 0,
      NH4 = 0,
      NPP = 0,
      prior = 0,
      Total = 0,
      HIST = rep(NA, max.iter)
    )
  )
}


#' @title Get Updated State
#' @author Thomas Bryce Kelly
#' @description Calculate group specific growth based on ambient conditions and model paramters
#' @param entry The index of the parameter set in the model history
get.state = function(model, entry = NULL) {
  
  ## If NULL then we just want to get current state, otherwise we either want to avg states or get a specific one:
  if (is.null(entry)) {
    params = model$params$current
  }else if (entry == 0) {
    
    ## Average all valid states
    state.orig = state
    state = init.state()
    state$score$HIST = state.orig$score$HIST
    n = max(which(!is.na(model$params$hist$Q))) - 1 ## First answer is always NEMURO-GOM
    
    ## Get first entry
    params = model$params$hist[2,] ## Skip 1 (NEMURO-GOM)
    names(params) = colnames(model$params$hist)
    
    state$growth$SYN =  GPP.SYN(params, model$obs.spreadsheet) / n
    state$growth$PRO =  GPP.PRO(params, model$obs.spreadsheet) / n
    state$growth$OTHER = GPP.OTHER(params, model$obs.spreadsheet) / n
    state$growth$DIA =  GPP.DIA(params, model$obs.spreadsheet) / n
    state$growth$DINO = GPP.DINO(params, model$obs.spreadsheet) / n
    state$growth$PRYM = GPP.PRYM(params, model$obs.spreadsheet) / n
    
    ## Run for all the other answers
    for (i in 2:n) {
      params = model$params$hist[i+1,]
      names(params) = colnames(model$params$hist)
      
      state$growth$SYN =  state$growth$SYN + GPP.SYN(params, model$obs.spreadsheet) / n
      state$growth$PRO =  state$growth$PRO + GPP.PRO(params, model$obs.spreadsheet) / n
      state$growth$OTHER = state$growth$OTHER + GPP.OTHER(params, model$obs.spreadsheet) / n
      state$growth$DIA =  state$growth$DIA + GPP.DIA(params, model$obs.spreadsheet) / n
      state$growth$DINO = state$growth$DINO + GPP.DINO(params, model$obs.spreadsheet) / n
      state$growth$PRYM = state$growth$PRYM + GPP.PRYM(params, model$obs.spreadsheet) / n
    }
    state$growth$NPP = data.table(SYN = state$growth$SYN$NPP, PRO = state$growth$PRO$NPP, OTHER = state$growth$OTHER$NPP,
                                 DIA = state$growth$DIA$NPP, DINO = state$growth$DINO$NPP, PRYM = state$growth$PRYM$NPP)
    
    state$growth$Rate = data.table(SYN = state$growth$SYN$Rate, PRO = state$growth$PRO$Rate, OTHER = state$growth$OTHER$Rate,
                                  DIA = state$growth$DIA$Rate, DINO = state$growth$DINO$Rate, PRYM = state$growth$PRYM$Rate)
    
    state$growth$TOTAL = data.table(GPP.NO3 = state$growth$SYN$GPP.NO3 + state$growth$PRO$GPP.NO3 + state$growth$OTHER$GPP.NO3 +
                                     state$growth$DIA$GPP.NO3 + state$growth$DINO$GPP.NO3 + state$growth$PRYM$GPP.NO3,
                                   GPP.NH4 = state$growth$SYN$GPP.NH4 + state$growth$PRO$GPP.NH4 + state$growth$OTHER$GPP.NH4 +
                                     state$growth$DIA$GPP.NH4 + state$growth$DINO$GPP.NH4 + state$growth$PRYM$GPP.NH4,
                                   NPP = state$growth$SYN$NPP + state$growth$PRO$NPP + state$growth$OTHER$NPP +
                                     state$growth$DIA$NPP + state$growth$DINO$NPP + state$growth$PRYM$NPP)
    state = calc.scores(model, state)
    return(state)
  } else {
    params = model$params$hist[entry,]
  }
  
  state$growth$SYN =  GPP.SYN(params, model$obs.spreadsheet)
  state$growth$PRO =  GPP.PRO(params, model$obs.spreadsheet)
  state$growth$OTHER = GPP.OTHER(params, model$obs.spreadsheet)
  state$growth$DIA =  GPP.DIA(params, model$obs.spreadsheet)
  state$growth$DINO = GPP.DINO(params, model$obs.spreadsheet)
  state$growth$PRYM = GPP.PRYM(params, model$obs.spreadsheet)
  state$growth$NPP = data.table(SYN = state$growth$SYN$NPP, PRO = state$growth$PRO$NPP, OTHER = state$growth$OTHER$NPP,
                                 DIA = state$growth$DIA$NPP, DINO = state$growth$DINO$NPP, PRYM = state$growth$PRYM$NPP)
  
  state$growth$Rate = data.table(SYN = state$growth$SYN$Rate, PRO = state$growth$PRO$Rate, OTHER = state$growth$OTHER$Rate,
                                DIA = state$growth$DIA$Rate, DINO = state$growth$DINO$Rate, PRYM = state$growth$PRYM$Rate)
  
  state$growth$TOTAL = data.table(GPP.NO3 = state$growth$SYN$GPP.NO3 + state$growth$PRO$GPP.NO3 + state$growth$OTHER$GPP.NO3 +
                                    state$growth$DIA$GPP.NO3 + state$growth$DINO$GPP.NO3 + state$growth$PRYM$GPP.NO3,
                                  GPP.NH4 = state$growth$SYN$GPP.NH4 + state$growth$PRO$GPP.NH4 + state$growth$OTHER$GPP.NH4 +
                                    state$growth$DIA$GPP.NH4 + state$growth$DINO$GPP.NH4 + state$growth$PRYM$GPP.NH4,
                                  NPP = state$growth$SYN$NPP + state$growth$PRO$NPP + state$growth$OTHER$NPP +
                                    state$growth$DIA$NPP + state$growth$DINO$NPP + state$growth$PRYM$NPP)
  
  state = calc.scores(model, state)
  state
}



#' @title Calculate Model Scores
#' @author Thomas Bryce Kelly
#' @description Calculate Scores based  on difference between the model result and observations
calc.scores = function(model, state) {
  total = state$growth$TOTAL ## Use total growth to compare to bulk measurements
  rate = state$growth$Rate
  epsilon = 1e-5
  
  state$score$NO3 = -0.5 * (abs(total$GPP.NO3 - model$obs.spreadsheet$NO3.Uptake) / (model$obs.spreadsheet$NO3.Uptake.Sigma + epsilon))^2 * -log(model$obs.spreadsheet$NO3.Uptake.Sigma * sqrt(2 * 3.14159))
  state$score$NH4 = 0 #(abs(total$GPP.NH4 - model$obs.spreadsheet$NH4.Uptake) / (model$obs.spreadsheet$NH4.Uptake.Sigma + epsilon))^2 * sum(-log(model$obs.spreadsheet$NH4.Uptake.Sigma * sqrt(2 * 3.14159)))
  state$score$NPP = -0.5 * (abs(total$NPP - model$obs.spreadsheet$NPP.N) / (model$obs.spreadsheet$NPP.N.Sigma + epsilon))^2 * -log(model$obs.spreadsheet$NPP.N.Sigma * sqrt(2 * 3.14159))
  
  ## Score for growth rates
  state$score$Mu.PRO = -0.5 * (abs(rate$PRO - model$obs.spreadsheet$PRO.Growth.Rate) / (0.16))^2 * -log(0.16 * sqrt(2 * 3.14159))
  state$score$Mu.SYN = -0.5 * (abs(rate$SYN - model$obs.spreadsheet$SYN.Growth.Rate) / (0.16))^2 * -log(0.16 * sqrt(2 * 3.14159))
  state$score$Mu.OTHER = -0.5 * (abs(rate$OTHER - model$obs.spreadsheet$OTHER.Growth.Rate) / (0.16))^2 * -log(0.16 * sqrt(2 * 3.14159))
  state$score$Mu.DIA = -0.5 * (abs(rate$DIA - model$obs.spreadsheet$DIA.Growth.Rate) / (0.16))^2 * -log(0.16 * sqrt(2 * 3.14159))
  state$score$Mu.DINO = -0.5 * (abs(rate$DINO - model$obs.spreadsheet$DINO.Growth.Rate) / (0.16))^2 * -log(0.16 * sqrt(2 * 3.14159))
  state$score$Mu.PRYM = -0.5 * (abs(rate$PRYM - model$obs.spreadsheet$PRYM.Growth.Rate) / (0.16))^2 * -log(0.16 * sqrt(2 * 3.14159))
  
  set = c(state$score$NO3, state$score$NH4, state$score$NPP, state$score$Mu.PRO, state$score$Mu.SYN, state$score$Mu.OTHER,
          state$score$Mu.DIA, state$score$Mu.DINO, state$score$Mu.PRYM)
  
  prior = rep(NA, length(model$params$prior))
  current = as.numeric(model$params$current)
  for (i in 1:length(model$params$current)) {
    if (model$params$prior[i] == 'delta') {
      prior[i] = delta(current[i], model$params$Param1[i]) - 1
    }
    if (model$params$prior[i] == 'dlnorm') {
      prior[i] = dlnorm(current[i], model$params$Param1[i], model$params$Param2[i], log = T)
    }
    if (model$params$prior[i] == 'dnorm') {
      prior[i] = dnorm(current[i], model$params$Param1[i], model$params$Param2[i], log = T)
    }
    if (model$params$prior[i] == 'dchisq') {
      prior[i] = dchisq(current[i], model$params$Param1[i], log = T)
    }
    if (model$params$prior[i] == 'dunif') {
      prior[i] = dunif(current[i], model$params$Param1[i], model$params$Param2[i], log = T)
    }
    
  }
  index = max(which(!is.na(state$score$HIST)), 0) + 1
  state$score$prior = prior
  state$score$Total = sum(set, na.rm = T) + sum(state$score$prior) ## Calculate log-likihood: P(x|theta) * P(theta)
  state$score$HIST[index] = state$score$Total
  
  state ## Return updated state
}



#' @title Accept New Solution?
#' @description A function used to determine if a new paramter set should be accepted by comparing an old to new log normal probabilities
#' @author Thomas Bryce Kelly
#' @param score.old The previous MCMC score
#' @param score.new The current, tentative MCMC score
accept.sol = function(score.old, score.new) {
  if (exp(score.new - score.old) > runif(1)) {
    return(TRUE)
  }
  FALSE
}


#' @title Get New Model Parameters
#' @author Thomas Bryce Kelly
#' @description Perform a random step from a given parameter set using the min and max bounds and step length saved in the model object
get.new.model = function(model, scale = 1) {
  
  ## Save old paramters
  index = max(which(!is.na(state$score$HIST)))
  model$params$hist[index,] = model$params$current
  for (i in 1:length(model$params$current)) {
    model$params$current[[i]] = rnorm(1, as.numeric(model$params$current)[i], model$params$step[i] * scale)
  }
  ## Return
  model
}


calc.TL = function(Temp, Q) {  # [unitless]
  exp(Temp * Q)
}

calc.NL = function(Nitrate, KNO3, Ammonium, KNH4) { #  [unitless]
  Nitrate / ((Nitrate + 10^KNO3) * (1 + Ammonium / 10^KNH4))
}

calc.AL = function(Ammonium, KNH4) { # [unitless]
  Ammonium / (Ammonium + 10^KNH4)
}

calc.LL = function(PAR, V, Alpha, Beta) { # [unitless]
  PAR = conv.einstein.to.watt(PAR * 1e-6) # conv [uE m-2 d-1] to [W m-2 d-1]
  (1 - exp(-exp(Alpha) * PAR / 10^V)) * exp(-Beta * PAR / 10^V)
}

calc.RESP = function(Temp, Q, R) { # [unitless]
  exp(R) * exp(Q * Temp)
}


#' @title Calculate Syn Growth
#' @author Thomas Bryce Kelly
#' @param params The model parameter set to use
#' @param obs The environmental observations to be used in the model (e.g. Nitrate, PAR)
GPP.SYN = function(params, obs) {
  
  # Unitless limitation factors
  TL = calc.TL(obs$Temp, params$Q)
  NL = calc.NL(Nitrate = obs$Nitrate, KNO3 = params$KNO3_SYN, Ammonium = obs$Ammonium, KNH4 = params$KNH4_SYN)
  AL = calc.AL(obs$Ammonium, params$KNH4_SYN)
  LL = calc.LL(obs$PAR, params$V_SYN, params$Alpha_SYN, params$Beta_SYN)
  
  GPP.NO3 = 10^(params$V_SYN) * TL * NL * LL * obs$SYN.Biomass# mmol N m-3 d-1 (via NO3)
  GPP.NH4 = 10^(params$V_SYN) * TL * AL * LL * obs$SYN.Biomass # mmol N m-3 d-1 (via NH4)
  GPP.NO3[is.na(GPP.NO3)] = 0 # Set NA to 0
  GPP.NH4[is.na(GPP.NH4)] = 0 # Set NA to 0
  
  RESP = calc.RESP(obs$Temp, params$Q, params$R_SYN) * obs$SYN.Biomass # mmol N m-3 d-1
  EXCR = params$E_SYN * (GPP.NO3 + GPP.NH4)  # mmol N m-3 d-1
  
  NPP = GPP.NO3 + GPP.NH4 - RESP - EXCR  # mmol N m-3 d-1
  Rate = NPP / obs$SYN.Biomass
  NPP[is.na(NPP)] = 0
  
  return(data.table(TL = TL, NL = NL, AL = AL, LL = LL, GPP.NO3 = GPP.NO3, GPP.NH4 = GPP.NH4, RESP = RESP, EXCR = EXCR, NPP = NPP, Rate = Rate))
}



#' @title Calculate Pro Growth
#' @author Thomas Bryce Kelly
#' @param params The model parameter set to use
#' @param obs The environmental observations to be used in the model (e.g. Nitrate, PAR)
GPP.PRO = function(params, obs) {
  
  # Unitless limitation factors
  TL = calc.TL(obs$Temp, params$Q)
  NL = calc.NL(obs$Nitrate, params$KNO3_PRO, obs$Ammonium, params$KNH4_PRO)
  AL = calc.AL(obs$Ammonium, params$KNH4_PRO)
  LL = calc.LL(obs$PAR, params$V_PRO, params$Alpha_PRO, params$Beta_PRO)
  
  GPP.NO3 = 10^(params$V_PRO) * TL * NL * LL * obs$PRO.Biomass
  GPP.NH4 = 10^(params$V_PRO) * TL * AL * LL * obs$PRO.Biomass
  GPP.NO3[is.na(GPP.NO3)] = 0 # Set NA to 0
  GPP.NH4[is.na(GPP.NH4)] = 0 # Set NA to 0
  
  RESP = calc.RESP(obs$Temp, params$Q, params$R_PRO) * obs$PRO.Biomass
  EXCR = params$E_PRO * (GPP.NO3 + GPP.NH4)
  
  NPP = GPP.NO3 + GPP.NH4 - RESP - EXCR
  Rate = NPP / obs$PRO.Biomass
  NPP[is.na(NPP)] = 0
  return(data.table(TL = TL, NL = NL, AL = AL, LL = LL, GPP.NO3 = GPP.NO3, GPP.NH4 = GPP.NH4, RESP = RESP, EXCR = EXCR, NPP = NPP, Rate = Rate))
}


#' @title Calculate OTHER-Euk Growth
#' @author Thomas Bryce Kelly
#' @param params The model parameter set to use
#' @param obs The environmental observations to be used in the model (e.g. Nitrate, PAR)
GPP.OTHER = function(params, obs) {
  
  # Unitless limitation factors
  TL = calc.TL(obs$Temp, params$Q)
  NL = calc.NL(obs$Nitrate, params$KNO3_OTHER, obs$Ammonium, params$KNH4_OTHER)
  AL = calc.AL(obs$Ammonium, params$KNH4_OTHER)
  LL = calc.LL(obs$PAR, params$V_OTHER, params$Alpha_OTHER, params$Beta_OTHER)
  
  GPP.NO3 = 10^(params$V_OTHER) * TL * NL * LL * obs$OTHER.Biomass
  GPP.NH4 = 10^(params$V_OTHER) * TL * AL * LL * obs$OTHER.Biomass
  GPP.NO3[is.na(GPP.NO3)] = 0 # Set NA to 0
  GPP.NH4[is.na(GPP.NH4)] = 0 # Set NA to 0
  
  RESP = calc.RESP(obs$Temp, params$Q, params$R_OTHER) * obs$OTHER.Biomass
  EXCR = params$E_OTHER * (GPP.NO3 + GPP.NH4)
  
  NPP = GPP.NO3 + GPP.NH4 - RESP - EXCR
  Rate = NPP / obs$OTHER.Biomass
  NPP[is.na(NPP)] = 0
  return(data.table(TL = TL, NL = NL, AL = AL, LL = LL, GPP.NO3 = GPP.NO3, GPP.NH4 = GPP.NH4, RESP = RESP, EXCR = EXCR, NPP = NPP, Rate = Rate))
}


#' @title Calculate Diatom Growth
#' @author Thomas Bryce Kelly
#' @param params The model parameter set to use
#' @param obs The environmental observations to be used in the model (e.g. Nitrate, PAR)
GPP.DIA = function(params, obs) {
  
  # Unitless limitation factors
  TL = calc.TL(obs$Temp, params$Q)
  NL = calc.NL(obs$Nitrate, params$KNO3_DIA, obs$Ammonium, params$KNH4_DIA)
  AL = calc.AL(obs$Ammonium, params$KNH4_DIA)
  LL = calc.LL(obs$PAR, params$V_DIA, params$Alpha_DIA, params$Beta_DIA)
  
  GPP.NO3 = 10^(params$V_DIA) * TL * NL * LL * obs$DIA.Biomass
  GPP.NH4 = 10^(params$V_DIA) * TL * AL * LL * obs$DIA.Biomass
  GPP.NO3[is.na(GPP.NO3)] = 0 # Set NA to 0
  GPP.NH4[is.na(GPP.NH4)] = 0 # Set NA to 0
  
  RESP = calc.RESP(obs$Temp, params$Q, params$R_DIA) * obs$DIA.Biomass
  EXCR = params$E_DIA * (GPP.NO3 + GPP.NH4)
  
  NPP = GPP.NO3 + GPP.NH4 - RESP - EXCR
  Rate = NPP / obs$DIA.Biomass
  NPP[is.na(NPP)] = 0
  
  return(data.table(TL = TL, NL = NL, AL = AL, LL = LL, GPP.NO3 = GPP.NO3, GPP.NH4 = GPP.NH4, RESP = RESP, EXCR = EXCR, NPP = NPP, Rate = Rate))
}


#' @title Calculate Other Growth
#' @author Thomas Bryce Kelly
#' @param params The model parameter set to use
#' @param obs The environmental observations to be used in the model (e.g. Nitrate, PAR)
GPP.DINO = function(params, obs) {
  
  # Unitless limitation factors
  TL = calc.TL(obs$Temp, params$Q)
  NL = calc.NL(obs$Nitrate, params$KNO3_DINO, obs$Ammonium, params$KNH4_DINO)
  AL = calc.AL(obs$Ammonium, params$KNH4_DINO)
  LL = calc.LL(obs$PAR, params$V_DINO, params$Alpha_DINO, params$Beta_DINO)
  
  GPP.NO3 = 10^(params$V_DINO) * TL * NL * LL * obs$DINO.Biomass
  GPP.NH4 = 10^(params$V_DINO) * TL * AL * LL * obs$DINO.Biomass
  GPP.NO3[is.na(GPP.NO3)] = 0 # Set NA to 0
  GPP.NH4[is.na(GPP.NH4)] = 0 # Set NA to 0
  
  RESP = calc.RESP(obs$Temp, params$Q, params$R_DINO) * obs$DINO.Biomass
  EXCR = params$E_DINO * (GPP.NO3 + GPP.NH4)
  
  NPP = GPP.NO3 + GPP.NH4 - RESP - EXCR
  Rate = NPP / obs$DINO.Biomass
  NPP[is.na(NPP)] = 0
  
  return(data.table(TL = TL, NL = NL, AL = AL, LL = LL, GPP.NO3 = GPP.NO3, GPP.NH4 = GPP.NH4, RESP = RESP, EXCR = EXCR, NPP = NPP, Rate = Rate))
}


#' @title Calculate PRYM Growth
#' @author Thomas Bryce Kelly
#' @param params The model parameter set to use
#' @param obs The environmental observations to be used in the model (e.g. Nitrate, PAR)
GPP.PRYM = function(params, obs) {
  
  # Unitless limitation factors
  TL = calc.TL(obs$Temp, params$Q)
  NL = calc.NL(obs$Nitrate, params$KNO3_PRYM, obs$Ammonium, params$KNH4_PRYM)
  AL = calc.AL(obs$Ammonium, params$KNH4_PRYM)
  LL = calc.LL(obs$PAR, params$V_PRYM, params$Alpha_PRYM, params$Beta_PRYM)
  
  GPP.NO3 = 10^(params$V_PRYM) * TL * NL * LL * obs$PRYM.Biomass
  GPP.NH4 = 10^(params$V_PRYM) * TL * AL * LL * obs$PRYM.Biomass
  GPP.NO3[is.na(GPP.NO3)] = 0 # Set NA to 0
  GPP.NH4[is.na(GPP.NH4)] = 0 # Set NA to 0
  
  RESP = calc.RESP(obs$Temp, params$Q, params$R_PRYM) * obs$PRYM.Biomass
  EXCR = params$E_PRYM * (GPP.NO3 + GPP.NH4)
  
  NPP = GPP.NO3 + GPP.NH4 - RESP - EXCR
  Rate = NPP / obs$PRYM.Biomass
  NPP[is.na(NPP)] = 0
  
  return(data.table(TL = TL, NL = NL, AL = AL, LL = LL, GPP.NO3 = GPP.NO3, GPP.NH4 = GPP.NH4, RESP = RESP, EXCR = EXCR, NPP = NPP, Rate = Rate))
}



#########################
## Plotting Functions ###
#########################



plot.limitations = function(state, col = get.pal(n = 8, pal = inferno)) {
  par(mfrow=c(2,2))
  
  ##TL
  plot(state$growth$SYN$TL, pch = 20, ylab = 'TL', ylim = c(0, 8), col = col[1])
  points(state$growth$PRO$TL, pch = 20, col = col[2])
  points(state$growth$OTHER$TL, pch = 20, col = col[3])
  points(state$growth$DIA$TL, pch = 20, col = col[4])
  points(state$growth$DINO$TL, pch = 20, col = col[5])
  points(state$growth$PRYM$TL, pch = 20, col = col[6])
  
  ##LL
  plot(state$growth$SYN$LL, pch = 20, ylab = 'LL', ylim = c(0, 1), yaxs = 'i', col = col[1])
  points(state$growth$PRO$LL, pch = 20, col = col[2])
  points(state$growth$OTHER$LL, pch = 20, col = col[3])
  points(state$growth$DIA$LL, pch = 20, col = col[4])
  points(state$growth$DINO$LL, pch = 20, col = col[5])
  points(state$growth$PRYM$LL, pch = 20, col = col[6])
  
  ##NL
  plot(state$growth$SYN$NL, pch = 20, ylab = 'NL', ylim = c(0, 1), yaxs = 'i', col = col[1])
  points(state$growth$PRO$NL, pch = 20, col = col[2])
  points(state$growth$OTHER$NL, pch = 20, col = col[3])
  points(state$growth$DIA$NL, pch = 20, col = col[4])
  points(state$growth$DINO$NL, pch = 20, col = col[5])
  points(state$growth$PRYM$NL, pch = 20, col = col[6])
  
  ##AL
  plot(state$growth$SYN$AL, pch = 20, ylab = 'AL', ylim = c(0, 1), yaxs = 'i', col = col[1])
  points(state$growth$PRO$AL, pch = 20, col = col[2])
  points(state$growth$OTHER$AL, pch = 20, col = col[3])
  points(state$growth$DIA$AL, pch = 20, col = col[4])
  points(state$growth$DINO$AL, pch = 20, col = col[5])
  points(state$growth$PRYM$AL, pch = 20, col = col[6])
  
}


plot.matrix = function(k){
  n = length(k)
  names = names(model$params$hist)
  par(mfcol = c(n,n), plt = c(0.25, 0.9, 0.2, 0.9))
  #par(mfcol = c(n+1,n+1), plt = c(0, 1, 0, 1))
  for (i in k) {
    for (j in k) {
      namei = strsplit(names[i], split = '_')[[1]][1]
      namej = strsplit(names[j], split = '_')[[1]][1]
      if (i == j) {
        #hist(model$params$hist[,i], xlab = names[i], main = '', yaxt = 'n', ylab = '')
        d = density(model$params$hist[-1,i], adjust = 0.2)
        plot(d$x, d$y, xlab = names[i], main = '', yaxt = 'n', ylab = '', frame.plot = F, type = 'l', col = 'darkgrey', xaxt = 'n', xlim = range(pretty(model$params$hist[,i])))
        polygon(x = c(d$x, rev(d$x)), y = c(d$y, rep(0, length(d$y))), col = '#00005530')
        if (which(i == k) == 4) {
          axis(1)
        } else {
          add.log.axis(1)
        }
        axis(1, at = c(-10,10))
        mtext(namei, line = -5, cex = 1.5, col = 'darkgrey')
        
        if (model$params$prior[i] == 'dnorm') { mc = rnorm(nrow(model$params$hist), model$params$Param1[i], model$params$Param2[i]) }
        if (model$params$prior[i] == 'dlnorm') { mc = rlnorm(nrow(model$params$hist), model$params$Param1[i], model$params$Param2[i]) }
        if (model$params$prior[i] == 'dunif') { mc = runif(nrow(model$params$hist), model$params$Param1[i], model$params$Param2[i]) }
        if (model$params$prior[i] == 'delta') { mc = rep(model$params$Param1[i], nrow(model$params$hist)) }
        d = density(mc)
        lines(d$x, d$y, col = 'blue', lwd = 2)
      }
      if (i < j) {
        smoothScatter(model$params$hist[-1,i], model$params$hist[-1,j], xlab = '', ylab = '', main = '', nrpoints = 0, xaxt = 'n', yaxt = 'n',
                      xlim = range(pretty(model$params$hist[,i])), ylim =  range(pretty(model$params$hist[,j])))
        mtext(namei, side = 1, cex = 0.8, line = 0.3)
        mtext(namej, side = 2, cex = 0.8, line = 0.3)
        correlation = cor.test(model$params$hist[-1,i], model$params$hist[-1,j])
        font = 2
        if (correlation$p.value > 0.05) { 
          font = 3
        }
        mtext(paste0(format(round(correlation$estimate*100)/100, digits = 2)), side = 1, adj = 0.9, cex = 0.8, col = 'darkred', line = -1.5, font = font)
      }
      if (i > j) {
        plot(NULL, NULL, xaxt = 'n', yaxt = 'n', ylab = '', xlab = '', frame.plot = F, xlim = c(0,1), ylim = c(0,1))
      }
    }
    #plot(density(model$params$hist[,i]), xlab = names[i], main = '', yaxt = 'n', ylab = '', frame.plot = F)
  }
  par(new = T, mfrow = c(1,1))
  plot(NULL, NULL, xaxt = 'n', yaxt = 'n', ylab = '', xlab = '', frame.plot = F, xlim = c(0,1), ylim = c(0,1))
  text(0.75, 0.75, strsplit(names[i], split = '_')[[1]][2], cex = 6)
}



add.violin = function(y, data, col1, col2, scale = 1, n = 1e4, border = F) {
  density = density(data)
  temp = approx(density$x, density$y, xout = seq(min(density$x), max(density$x), length.out = n))
  density$x = temp$x
  density$y = temp$y
  xx = c(density$x, rev(density$x))
  yy  = y + scale * c(density$y, -rev(density$y)) / mean(density$y) / 2
  if (border) {
    polygon(x = xx, y = yy, border = col1, col = NA)
  } else {
    polygon(x = xx, y = yy, col = col1, border = NA)
    add.boxplot.box(y, data, col = col2, width = scale*0.8, outliers = F)
  }
}


add.boxplot.box = function(x, y, col = 'grey', border = 'black', width = 0.7, lty = 1, lcol = 'black',
                            lwd = 1, xlwd = NULL, outliers = T, pcol = 'black', pch = 1, cex = 1) {
  if (length(x) == 1) { x = rep(x, length(y))}
  for (xx in unique(x)) {
    l = which(x == xx)
    if(is.null(xlwd)) { xlwd = width / 3 }
    
    ## statistics
    q1 = quantile(y[l], probs = 0.25, na.rm = T)
    q3 = quantile(y[l], probs = 0.75, na.rm = T)
    iqr = IQR(y[l], na.rm = TRUE)
    m = median(y[l], na.rm = TRUE)
    
    ## Box
    rect(ybottom = xx - width/2, xleft = q1, ytop = xx + width/2, xright = q3, col = col, border = border)
    lines(y = c(xx - width/2, xx + width/2), x = rep(m,2)) # Horizontal
    
    ## Add outliers
    k = which(y[l] < q1 - 1.5 * iqr | y[l] > q3 + 1.5 * iqr)
    if (length(k) > 0 & outliers) { points(y = rep(xx, length(k)), x = y[l[k]], pch = pch, col = pcol, cex = cex) }
    
    ## Add whiskers
    if (length(k) > 0) {
      lines(y = rep(xx, 2), x = c(q1, min(y[l[-k]])), col = lcol, lwd = lwd)
      lines(y = rep(xx, 2), x = c(q3, max(y[l[-k]])), col = lcol, lwd = lwd)
      lines(y = c(xx - xlwd/2, xx + xlwd/2), x = rep(min(y[l[-k]]), 2), col = lcol, lwd = lwd)
      lines(y = c(xx - xlwd/2, xx + xlwd/2), x = rep(max(y[l[-k]]), 2), col = lcol, lwd = lwd)
    } else {
      lines(y = rep(xx, 2), x = c(q1, min(y[l])), col = lcol, lwd = lwd)
      lines(y = rep(xx, 2), x = c(q3, max(y[l])), col = lcol, lwd = lwd)
      lines(y = c(xx - xlwd/2, xx + xlwd/2), x = rep(min(y[l]), 2), col = lcol, lwd = lwd)
      lines(y = c(xx - xlwd/2, xx + xlwd/2), x = rep(max(y[l]), 2), col = lcol, lwd = lwd)
    }
  }
}