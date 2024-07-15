
##1. IPW estimator for binary treatment
ipwEstimator_binaryA <- function(formula, data, numerCond = ~ 1, ...){ #If u want to fit for f(A|V)/f(A|L,V), modify numerCond.
  # IP Weighting estimator. formula should have the format ~ C1 + C2 + ... 
  denomFormula <- update(formula, A ~ .)
  numerFormula <- update(numerCond, A ~ .)
  probDenom <- predict(glm(formula = denomFormula, family = binomial(), data = data),
                       type = "response")
  probNumer <- predict(glm(formula = numerFormula, family = binomial(), data = data),
                       type = "response")
  IPW <- ifelse(data$A, probNumer/probDenom, (1 - probNumer)/(1 - probDenom))
  model <- glm(Y ~ A, family = gaussian(), weights = IPW, data = data)
  return(list(est = coef(model)[2], IPW = IPW)) #!! return a list. IPWs are returned for checking balance using SMD.
}

##2. IPW estimator for continuous treatment
ipwEstimator_continuousA <- function(formula, data, numerCond = ~ 1, ...){ #If u want to fit for f(A|V)/f(A|L,V), modify cond.
  denomFormula <- update(formula, A ~ .)
  numerFormula <- update(numerCond, A ~ .)
  probdensDenom <- getProbDensA(distr_A, denomFormula, data)
  probdensNumer <- getProbDensA(distr_A, numerFormula, data, numer = TRUE)
  IPW <- probdensNumer/probdensDenom
  model <- glm(Y ~ A, family = gaussian(), weights = IPW, data = data)
  return(list(est = coef(model)[2], IPW = IPW)) #!! return a list. IPWs are returned for checking balance using SMD.
}

##3. Outcome regression estimator.
outcomeRegEstimator <- function(formula, data){
  # Outcome regression approach. formula should have the format ~ C1 + C2 + ...
  model <- glm(formula = update(formula, Y ~ A + .),family = gaussian(),data = data)
  return(coef(model)[2])
}

##6. Doubly robust estimator (tmle::tmle), binary treatment
drEstimator_binaryA <- function(formula, data, return = c("ATE","ATT","ATC"),
                                Q.SL.library = "SL.glm", g.SL.library = "SL.glm", obsWeights = NULL,
                                ...){
  Qform <- paste0("Y ~ A + ", as.character(formula)[2])
  gform <- paste0("A ~ ", as.character(formula)[2])
  TmleEst <- tmle::tmle(Y = data$Y, A = data$A, W = subset(data, select = -c(Y,A)),
                        Qform = Qform, gform = gform, obsWeights = obsWeights,
                        Q.SL.library = Q.SL.library, g.SL.library = g.SL.library)
  return(TmleEst$estimates[[return]]$psi)
  
}

drEstimator_continuousA <- function(...){
  return(NA) # Not implemented yet
}