# Dependencies: truncnorm

getProbDensA <- function(distr, formula, data, numer = FALSE,...){ 
  # get expected probability density for a continuous variable.
  # Specially for continuous treatment A.
  
  ProbDens <- switch(distr,
                     "rnorm" = {
                       predictValuesNorm <- predict.glm(glm(formula = formula, data = data, family = gaussian()),type = "response")
                       sd <- sqrt(var_A) # use "true" sd.
                       dnorm(x = data$A,
                                     mean = predictValuesNorm,
                                     sd = sd)
                       }, # assume homoskedasticity.
                     "rtruncnorm" = {
                       predictValuesNorm <- predict.glm(glm(formula = formula, data = data, family = gaussian()),type = "response")
                       sd <- sqrt(var_A) # use "true" sd.
                       dtruncnorm(x = data$A,
                                  a = min(predictValuesNorm),
                                  b = max(predictValuesNorm),
                                  sd = sd,
                                  mean = predictValuesNorm)
                       }, # approx. 
                     "runif" = NULL, # Not implemented. Misspecified models in glm.
                     "rbeta" = NULL, # will use betareg::betareg to implement beta glm. --!need to scale the variable back to (0,1). Not implemented yet.
                     "rlnorm" = { # tbd.
                       probDensLnorm <- function(){
                         predictValues <- predict.glm(glm(formula = formula, data = data, family = gaussian(link = "log")),type = "response")
                         sdnormValues <- sqrt(var_A)
                         dlnorm() # tbd. 
                       }
                       probDensLnorm()
                     }, # gaussian(link = "log")
                     "rgamma" = {
                       probDensGammaGen <- function(){
                         predictValues <- predict.glm(glm(formula = formula,data = data, family = Gamma(link = "identity")),type = "response")
                         varValues <- sqrt(var_A) # assume homoskedasticity, because we know our gamma-shaped treatment is simulated in a homoskedastic way...?
                         # But Gamma glm assumes heteroskedasticity by definition! --How to fix it?
                         return(dgamma(x = data$A, shape = (predictValues^2)/varValues, scale = abs(varValues/predictValues))) # shape and scale taken to be as calculated from predicted values
                       }
                       probDensGammaGen()
                     }, # Gamma(link = "identity")
                     "rchisq" = NULL, #tbd
                     stop("Specified distribution for A not supported yet when estimating probability density")
  )
  return(ProbDens)
}