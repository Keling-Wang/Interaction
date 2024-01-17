# Library packages
#library(rlang)
#library(pryr)
#library(truncnorm)
#library(Hmisc)

# The DAG representing causal relationships here can be found in corresponding
# explanations.
# A (treat): binary. I: binary, 0: Level M, 1: Level N. C1, C2, V, Y: continuous.

# |                       C1
# |                      /  \
# |           I ------> /    \
# |           |        /      \
# |           |       /        \
# |           |     `//`      ` \\`
# |           | V-> A ---------> Y
# |           |     `\\`      `//`
# |           |------>\        /
# |                    \      /
# |                     \    /
# |                      \  /
# |                       C2

simulationInteraction <- function(sim_runs = 7500L, size = 15000L,
                                  distr_C1 = c("rnorm","rtruncnorm","runif","rbeta","rlnorm","rgamma","rchisq"),
                                  center_C1, var_C1, min_C1, max_C1,
                                  distr_C2, center_C2, var_C2, min_C2, max_C2,
                                  distr_V, center_V, var_V, min_V, max_V,
                                  Intr = c(0,0,0,0), #RR_C_A per 1 unit C in different strata of I
                                  mean_Y, var_Y,
                                  A_prob = 0.5, I_prob = 0.5,
                                  E_A_Y, E_C1_Y, E_C2_Y, RR_V_A,
                                  ...){
  # Define a function to simulate interactions. 
  # The function takes several parameters:
  # :sim_runs: integer. the number of simulation runs (default is 7500)
  # :size: integer. the size of the simulation (default is 15000)
  # :distr_C1, distr_C2, distr_V: character. the distributions for C1, C2, and V 
  # :center_C1, center_C2, center_V: numeric. the (approx.) centers for C1, C2, and V
  # :var_C1, var_C2, var_V: numeric. positive. the (approx.) variances for C1, C2, and V
  # :min_C1, max_C1, min_C2, max_C2, min_V, max_V: numeric. Min and max value for variables. Note: under some distributions there are by definition no max. and min.
  # :Intr: c("Intr_LevelM_C1","Intr_LevelM_C2","Intr_LevelN_C1","Intr_LevelN_C2"). Strength of effects of C-s on treatment for different level of the variable interacting with C-A effect.
  # :mean_Y, var_Y: the mean and variance for Y given C1=C2=A=0. Y is normally distributed here.
  # :A_prob, I_prob: marginal treatment probability without adjustment; proportion of I(level=="M").
  # :E_A_Y, E_C1_Y, E_C2_Y, RR_V_A: (Causal) effects for A_Y, C1_Y, C2_Y, V_A
  # ...: additional parameters. Current available option: formula_list: list of formulas (with only righthand part).
  
  # The function will return the follows:
  # :result: class "BiasSimulation" object. A list with following sublists:
  # :result$bias_table: a data frame with bias and root mean squared error of each formula under IPW and OR
  # :result$causal_structure: a list of current treatment/response/confounders/other variable notations.
  # :result$call: original function call
  # :result$args: argument list being passed to the function
  # :result$formulas: formula list used.
  # :result$simdata: list of results in each simulation run and A raw data sample.
  # :result$rawdatadescription: a look of raw data (one sample).
  
  #check arguments and store original call
  call <- match.call()
  optional_args <- list(...)
  arguments <- c(as.list(environment()),list(...))
  distr_options <- c("rnorm","rtruncnorm","runif","rbeta","rlnorm","rgamma","rchisq") # available distribution options
  distr_C1 <- match.arg(distr_C1,distr_options)
  distr_C2 <- match.arg(distr_C2,distr_options)
  distr_V <- match.arg(distr_V,distr_options)
  distr_Y <- "rnorm" #Y is by definition normally distributed
  
  generateContiVars <- function(distr, center, var, 
                           minimum, maximum,
                           size){ 
    # generate continuous variables using given characteristics.
    # One can add some constraints to the max. and min. of the data.
    variable_vec <- switch(distr,
                           "rnorm" = rnorm(n = size, mean = center, sd = sqrt(var)), # only center and var
                           "rtruncnorm" = rtruncnorm(n = size, a = minimum, b = maximum, mean = center, sd = sqrt(var)),
                           # center, var, max, min
                           "runif" = runif(n = size, min = minimum, max = maximum), # only max and min
                           "rbeta" = {
                             rbetaGen <- function(){ # max, min, center, var. will only be centered approximately.
                              total_range <- maximum - minimum 
                              scaled_mean <- (center - minimum) / total_range
                              scaled_var <- var / (total_range ^ 2)
                             
                              alpha <- ((1 - scaled_mean) / scaled_var - 1 / scaled_mean) * scaled_mean ^ 2
                              beta <- alpha * (1 / scaled_mean - 1)
                              beta_distribution <- rbeta(size, alpha, beta)
                              return(minimum + beta_distribution * total_range)
                             }
                             rbetaGen()
                           },
                           "rlnorm" = { # only center and var
                             rlnormGen <- function(){
                               sdlog <- sqrt(log(var/center^2 + 1))
                               meanlog <- log(center) - 0.5 * sdlog^2
                               lognormal_distribution <- rlnorm(size, meanlog, sdlog)
                               return(lognormal_distribution)
                             }
                             rlnormGen()
                           },
                           "rgamma" = rgamma(size, shape = (center^2) / var, scale = var / center), # only center and var
                           "rchisq" = rchisq(size,df = center), # only center
                           stop("???")
    )
    return(variable_vec)
  }
  
  ipwEstimator <- function(formula, data){
    # IP Weighting estimator. formula should have the format ~ C1 + C2 + ... 
      prob <- predict(glm(formula = update(formula, A ~ . ),family = binomial(),data = data),
              type = "response")
      overall_treat_prob <- sum(data$A)/NROW(data)
      IPW <- ifelse(data$A, overall_treat_prob/prob, (1-overall_treat_prob)/(1-prob))
      model <- glm(Y ~ A, family = gaussian(),weights = IPW, data = data)
      return(coef(model)[2])
  }
  
  outcomeRegEstimator <- function(formula, data){
    # Outcome regression approach. formula should have the format ~ C1 + C2 + ...
    model <- glm(formula = update(formula, Y ~ A + .),family = gaussian(),data = data)
    return(coef(model)[2])
  }

  buildDataSim <- function(){
    data_sim <- data.frame(
      C1 = generateContiVars(distr_C1, center_C1, var_C1, min_C1, max_C1, size),
      C2 = generateContiVars(distr_C2, center_C2, var_C2, min_C2, max_C2, size),
      V = generateContiVars(distr_V, center_V, var_V, min_V, max_V,size)
    )
    
    # yield I (variable interacting with C1/C2-A effect)
    data_sim$I <- rbinom(size,1,I_prob)
    
    # Convert RR to OR, Pr(A|C1=0,C2=0,Intr[1:4]=0,V=0) is used.
    OR_Intr <- Intr * (1 - A_prob) / (1 - Intr * A_prob)
    OR_V_A <- RR_V_A * (1 - A_prob) / (1 - RR_V_A * A_prob)
    
    # yield Treat (A) 
    LogOdds_A <- with(data_sim, log(A_prob / (1 - A_prob)) + 
                        (log(OR_Intr[1])*C1 + log(OR_Intr[2])*C2) * (I==0) +
                        (log(OR_Intr[3])*C1 + log(OR_Intr[4])*C2) * (I==1) +
                        log(OR_V_A*V)
    )
    
    post_prob_A <- exp(LogOdds_A)/(exp(LogOdds_A)+1)
    data_sim$A <- rbinom(size, 1, post_prob_A)
    
    
    # yield Y
    mu <- with(data_sim,mean_Y + E_A_Y * A + E_C1_Y * C1 + E_C2_Y * C2)
    data_sim$Y <- rnorm(size,mean = mu,sd = sqrt(var_Y))
    
    return(data_sim)
  }
  
  
  
  # define situations to be evaluated. You can either define yourself formula_list outside the function.
  if(exists("formula_list", where = sys.frame(which = 1L)) == FALSE){
    # If a formula list (with exact name formula_list) is defined outside the function and passed into the call, the part below will never run.
    formula_list <- c(
      ~ C1 + C2, # No interaction considered
      ~ (C1 + C2) * I - I, # Only interaction effects considered, without I itself
      ~ (C1 + C2) * I, # "Full Model"
      ~ C1 + C2 + V # adjusted for V. possible misspecified situation 1. IPW using this formula will estimate conditional effect w.r.t. V.
    )
  }
  
  
  raw_results_IPW <- raw_results_OR <- matrix(NA,nrow = sim_runs,ncol = length(formula_list))
  cat("\nSimulation running... Please wait... \n")
  pb <- txtProgressBar(min = 0, max = sim_runs, initial = 0, style = 3) # Progress bar

  for(i in 1:sim_runs){ # Run
    set.seed(2024+i)
    data_sim <- buildDataSim()
    for(f in 1:length(formula_list)){ # loop over all formulae
      raw_results_IPW[i,f] <- ipwEstimator(formula = formula_list[f][[1]],data = data_sim) #IPW estimates
      raw_results_OR[i,f] <- outcomeRegEstimator(formula = formula_list[f][[1]],data = data_sim) #Outcome Regression estimates
    }
    setTxtProgressBar(pb,i)
  }
  close(pb) # Close progress bar
  cat("\n Simulation done!\n")
  
  # prepare output
  results <- data.frame()
  for(i in 1:length(formula_list)){ # Calculate bias+RMSE and return main table.
    name = as.character(formula_list[i][[1]])[2]
    bias_IPW = mean(raw_results_IPW[,i]) - E_A_Y
    RMSE_IPW = sqrt(mean((raw_results_IPW[,i] - E_A_Y)^2))
    bias_OR = mean(raw_results_OR[,i]) - E_A_Y
    RMSE_OR = sqrt(mean((raw_results_OR[,i] - E_A_Y)^2))
    results <- rbind.data.frame(results, data.frame(name,bias_IPW,RMSE_IPW,bias_OR,RMSE_OR))
  }
  names(results) <- c("name","bias_IPW","RMSE_IPW","bias_OR","RMSE_OR")
  
  causal_structure <- list(
    treatment = "A",
    response = "Y",
    confounder = c("C1","C2"),
    interactions = c("I","I:RR(A|C1)","I:RR(A|C2)"),
    other_cause_tr = "V",
    other_cause_tr = NULL
  )
  
  returnlist <- list(bias_table = results, 
                     causal_structure = causal_structure,
                     call = call, 
                     args = arguments,
                     formulas = formula_list,
                     simdata = list(IPW = raw_results_IPW, OR = raw_results_OR,raw = buildDataSim()),
                     rawdatadescription = Hmisc::describe(buildDataSim())
                     )
  class(returnlist) <- "BiasSimulation"
  attr(returnlist,"simulation") <- c("runs" = sim_runs,"size" = size)
  return(returnlist)
}

summary.BiasSimulation <- function(x){
  stopifnot(class(x)=="BiasSimulation")
  cat("Simulation of a causal structure with confounding interaction:\n")
  cat("Simulation runs ", attr(x,"simulation")[1], "times; Sample size per run: ",attr(x,"simulation")[2],"\n")
  cat("Formula used: \n")
  for(i in 1:length(x$formulas)){
    cat(" ", as.character(x$formulas[i][[1]])," ;")
  }
  cat("\n\n")
  cat("Treatment: ", x$causal_structure$treatment[1], " Response: ", x$causal_structure$response[1],"\n")
  cat("True causal effects: ", x$args$E_A_Y, "\n\n")
  cat("Interaction: \n")
  cat("Pr(A|C1=c,I=0)/Pr(A|C1=c-1,I=0) = ", x$args$Intr[1], ";  Pr(A|C2=c,I=0)/Pr(A|C2=c-1,I=0) = ", x$args$Intr[2], "\n")
  cat("Pr(A|C1=c,I=1)/Pr(A|C1=c-1,I=1) = ", x$args$Intr[3], ";  Pr(A|C2=c,I=1)/Pr(A|C2=c-1,I=1) = ", x$args$Intr[4], "\n")
  cat("Noise (V): RR = ", x$args$RR_V_A,"\n\n\n")
  cat("Bias table: \n")
  defaultdigits <- getOption("digits")
  options(digits = 4)
  print(x$bias_table)
  options(digits = defaultdigits)
}

print.BiasSimulation <- function(x){
  NextMethod("print",x)
}

describeBiasSimulation <- function(x){
  x$rawdatadescription
}

makeBiasTable <- function(x, ...) UseMethod("makeBiasTable")

makeBiasTable.default <- function(x, ...){
  datarows0 <- c(x$args$RR_V_A,x$args$Intr)
  returndata <- data.frame()
  
  for(i in seq_len(NROW(x$bias_table))){
    datarows <- c(datarows0,x$bias_table[i,1],x$bias_table[i,2],x$bias_table[i,2]-x$bias_table[i,3],x$bias_table[i,2]+x$bias_table[i,3])
    returndata <- rbind.data.frame(returndata,datarows)
  }
  colnames(returndata) <- c("RR_V_A","Intr1","Intr2","Intr3","Intr4","Name","Bias","low","high")
  return(returndata)
}

makeBiasTable.BiasSimulation <- function(x, ...) NextMethod("makeBiasTable")


# other aux. functions

inflInteractionStrength <- function(Intr,steplength){ # define a function to update Intr
  current_C1 <- Intr[1]/Intr[3]
  current_C2 <- Intr[2]/Intr[4]
  Intr[1] <- (1+steplength)*Intr[1]
  Intr[2] <- (1+steplength)*Intr[2]
  return(Intr)
}
