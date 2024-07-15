# Dependencies: rlang, pryr, truncnorm, Hmisc, betareg, tmle, lmw, cobalt


simulationInteraction <- function(sim_runs = 1000L, size = 15000L,
                                  distr_C1 = c("rnorm","rtruncnorm","runif","rbeta","rlnorm","rgamma","rchisq"),
                                  center_C1, var_C1, min_C1, max_C1,
                                  distr_C2, center_C2, var_C2, min_C2, max_C2,
                                  distr_VA, center_VA, var_VA, min_VA, max_VA,
                                  distr_VY, center_VY, var_VY, min_VY, max_VY,
                                  IntrA = c(0,0,0,0), RR_IA_A = 1, E_IA_A = 0,
                                  IntrY = c(0,0,0,0), E_IY_Y = 0,
                                  IA_prob = 0.2, IY_prob = 0.5,
                                  mean_Y, var_Y,
                                  A_prob = 0.5, distr_A = NULL, mean_A = NA, var_A = NA, min_A = -Inf, max_A = Inf,
                                  E_A_Y, RR_VA_A, E_VA_A, E_VY_Y,
                                  continuous_A = FALSE,
                                  numerCond = ~ 1,
                                  bal_abs_scale = FALSE,
                                  use_DRestimator = TRUE,
                                  ...){
  # This is the main function.
  # The arguments for the function:
  # :sim_runs: integer. times of simulation runs
  # :size: integer. size of simulation dataset for each run of simulation
  # :distr_, center_, var_, min_, max_ for C1, C2, VA, VY, A, Y: double. distributions and parameters for continuous variables
  # :IntrA: vector of doubles of length 4. Effect of C1 and C2 on A in different strata of IA. IntrA[1] and [2] = (E_C1_A and E_C2_A) or (RR_C1_A and RR_C2_A) for IA=0; [3] and [4] for IA=1.
  # :IntrY: vector of doubles of length 4. Effect of C1 and C2 on Y in different strata of IY. Same as IntrA.
  # :IA_, IY_, A_prob: (prior) probability for binary variables
  # :RR_IA_A, RR_VA_A: double. Greater than 0. Risk ratio of exposures on a binary variable. per 1 unit increase (log scale) for continuous exposures.
  # :E_IA_A, E_IY_Y, E_VA_A, E_VY_Y: double.Effect of exposures on a continuous variable. per 1 unit increase for continuous exposures.
  # :E_A_Y: double. True causal effects.
  # :continuous_A: whether you want to explore a continuous exposure. Default set to FALSE.
  #               ** Although I write this down, I still never run it because the IP weights are too unstable and often result in huge bias (your arguments are right).
  # :numerCond: The formula that will be used if conditioning on some variables in IPTW. i.e., should be ~ V if the IPW is in the form Pr(A|V)/Pr(A|L,V) (and the formula for denominator should be ~ L + V)
  # :bal_abs_scale: whether you want to take only absolute value of each balance score (SMDs or Corr). This may influence the overall mean of balance scores. Default set to FALSE.
  # :use_DRestimator: whether you want the function to return estimates from doubly robust estimators (tmle::tmle)
  # ...: other arguments to be passed into the function call. Currently only "formula_list" is used.
  # ..formula_list: a list of formulas with only right-hand side (Leave treatment and outcome out). Will be used to fit models.
  
  # The function will return the follows:
  # :result: class "BiasSimulation" object. A list of length 8, with following sublists:
  # :result$bias_table: a data frame with bias and root mean squared error of each formula under IPW and OR
  # :result$balance: a list of matrices with balance scores for IPW and OR estimators. Will return the mean value and the standard error of balance scores.
  # :result$causal_structure: a list of current treatment/response/confounders/interactions/other variable notations.
  # :result$call: original function call
  # :result$args: argument list passed into the function call
  # :result$formulas: formula list used.
  # :result$simdata: list of results in each simulation run and A raw dataset sample.
  # :result$rawdatadescription: a look of raw data (one sample).
  
  #check arguments and store original call
  
  arguments <- c(as.list(environment()),list(...))
  optional_args <- list(...)
  call <- match.call()
  distr_options <- c("rnorm","rtruncnorm","runif","rbeta","rlnorm","rgamma","rchisq") # available distribution options for C1 and C2
  # rnorm: normal distribution; rtruncnorm: truncated normal distribution (using truncnorm package)
  # runif: uniform distribution; rbeta: beta distribution (scaled); rlnorm: log-normal distribution;
  # rgamma: Gamma distribution; rchisq: Chi-sq distribution.
  distr_C1 <- match.arg(distr_C1,distr_options)
  distr_C2 <- match.arg(distr_C2,distr_options)
  distr_VA <- match.arg(distr_VA,distr_options)
  distr_VY <- match.arg(distr_VY,distr_options)
  distr_Y <- "rnorm" #Y is by definition normally distributed
  
  if(continuous_A){ # simple check for continuous A.
    stopifnot(exists("distr_A",where = sys.frame(which = 1L)))
    distr_A <- match.arg(distr_A,c("rnorm","rtruncnorm","rgamma")) # currently only norm and gamma A is possible.
  }
  
  
  # define formula lists. You can either define yourself formula_list outside the function.
  if(exists("formula_list", where = sys.frame(which = 1L)) == FALSE){
    # If a formula list (with exact name formula_list) is defined outside the function and passed into the call, the part below will never run.
    formula_list <- c(
      ~ C1 + C2, # No interaction considered
      ~ (C1 + C2) * IA - IA, # Only interaction effects considered, without IA itself
      ~ (C1 + C2) * IA, # "Full Model"
      ~ C1 + C2 + VA # adjusted for VA. possible misspecified situation 1. 
    )
  }
  
  # set blank data frames and matrices
  raw_results_IPW <- raw_results_OR <- matrix(NA,nrow = sim_runs,ncol = length(formula_list))
  balance_IPW <- balance_OR <- vector("list", length = length(formula_list))
  raw_results_DR <- matrix(NA,nrow = sim_runs,ncol = length(formula_list))

  # set progress bar
  cat("\nSimulation running... Please wait... \n")
  pb <- txtProgressBar(min = 0, max = sim_runs, initial = 0, style = 3) # Progress bar

  
  for(i in 1:sim_runs){ # Loop over 1:sim_runs * I didn't use seq_len() because edge cases will not present...
    set.seed(2024+i) # set different seeds.
    data_sim <- buildDataSim(continuous_A = continuous_A) # build simulation dataset.
    for(f in 1:length(formula_list)){ # loop over all formulae
      if(!continuous_A){ # binary A simulation process. Estimates and Balance.
        IPWEsts <- ipwEstimator_binaryA(formula = formula_list[f][[1]],data = data_sim,numerCond = numerCond)
        OREsts <- outcomeRegEstimator(formula = formula_list[f][[1]],data = data_sim)
        raw_results_IPW[i,f] <- IPWEsts$est #IPW estimates
        raw_results_OR[i,f] <- OREsts #Outcome Regression estimates
        if(use_DRestimator){ # DR estimates
          DREsts <- drEstimator_binaryA(formula = formula_list[f][[1]], data = data_sim, return = "ATE")
          raw_results_DR[i,f] <- DREsts
        } else {
          raw_results_DR[i,f] <- NA
        }
        
        # balance scores
        
        BalanceCheckIPW <- balanceCheck(W = IPWEsts$IPW, data = data_sim, abs_scale = bal_abs_scale, continuous_A = continuous_A)
        BalanceScoreOR <- weightImputationOR(formula = formula_list[f][[1]], data = data_sim) # use lmw::lmw to impute weights for outcome regression estimates.
        BalanceCheckOR <- balanceCheck(W = BalanceScoreOR, data = data_sim, abs_scale = bal_abs_scale, continuous_A = continuous_A)
        if(i==1){
          balance_IPW[[f]] <- BalanceCheckIPW
          balance_OR[[f]] <- BalanceCheckOR
        } else {
          balance_IPW[[f]] <- cbind(balance_IPW[[f]],BalanceCheckIPW$balance)
          balance_OR[[f]] <- cbind(balance_OR[[f]],BalanceCheckOR$balance)
        }
        
      } else { # continuous A case. Only report Balance for IPW estimators. 
        IPWEsts <- ipwEstimator_continuousA(formula = formula_list[[f]],data = data_sim,numerCond = numerCond)
        OREsts <- outcomeRegEstimator(formula = formula_list[[f]],data = data_sim)
        raw_results_IPW[i,f] <- IPWEsts$est #IPW estimates
        raw_results_OR[i,f] <- OREsts #Outcome Regression estimates
        if(use_DRestimator){ # DR estimates
          DREsts <- drEstimator_continuousA(formula = formula_list[[f]], data = data_sim)
          raw_results_DR[i,f] <- DREsts # tbd.... now only returns NA
        } else {
          raw_results_DR[i,f] <- NA
        }
        # Check balance. OR only returns NA because there's no function to check balance of OR in continuous cases.
        BalanceCheckIPW <- balanceCheck(W = IPWEsts$IPW, data = data_sim, abs_scale = bal_abs_scale, continuous_A = continuous_A)
        if(i==1){
          balance_IPW[[f]] <- BalanceCheckIPW
          balance_OR[[f]] <- data.frame(BalanceCheckIPW$varname,rep(NA,NROW(BalanceCheckIPW)))
        } else {
          balance_IPW[[f]] <- cbind(balance_IPW[[f]],BalanceCheckIPW$balance)
          balance_OR[[f]] <- cbind(balance_OR[[f]],rep(NA,NROW(BalanceCheckIPW)))
        }
      }
    }
    setTxtProgressBar(pb,i) # update progress bar
  }
  close(pb) # Close progress bar
  cat("\n Simulation done!\n")
  
  # prepare output
  results <- data.frame()
  for(i in 1:length(formula_list)){ # Calculate bias+RMSE and return main table.
    name = as.character(formula_list[i][[1]])[2]
    bias_IPW = mean(raw_results_IPW[,i]) - E_A_Y # bias = mean - true_eff
    RMSE_IPW = sqrt(mean((raw_results_IPW[,i] - E_A_Y)^2)) # root mean squared error
    bias_OR = mean(raw_results_OR[,i]) - E_A_Y
    RMSE_OR = sqrt(mean((raw_results_OR[,i] - E_A_Y)^2))
    if(use_DRestimator){
      bias_DR <- mean(raw_results_DR[,i]) - E_A_Y
      RMSE_DR <- sqrt(mean((raw_results_DR[,i] - E_A_Y)^2))
    } else {
      bias_DR <- RMSE_DR <- NA
    }
    results <- rbind.data.frame(results, data.frame(name,bias_IPW,RMSE_IPW,bias_OR,RMSE_OR, bias_DR, RMSE_DR))
  }
  names(results) <- c("name","bias_IPW","RMSE_IPW","bias_OR","RMSE_OR", "bias_DR", "RMSE_DR")
  
  # preparing balance score output
  balance_list <- vector("list",length = length(formula_list))
  for(i in 1:length(formula_list)){
    balance_list[[i]] <- matrix(nrow = NROW(balance_IPW[[i]]),ncol = 4)
    names(balance_list)[i] <- as.character(formula_list[[i]])[2]
    rownames(balance_list[[i]]) <- balance_IPW[[i]][,1]
    colnames(balance_list[[i]]) <- c("balance_IPW","balanceSE_IPW","balance_OR","balanceSE_OR")
    balance_list[[i]][,1] <- rowSums(balance_IPW[[i]][,-1])/sim_runs # calculate "mean balance score"
    balance_list[[i]][,3] <- rowSums(balance_OR[[i]][,-1])/sim_runs
    for(j in 1:NROW(balance_IPW[[i]])){
      balance_list[[i]][j,2] <- sd(balance_IPW[[i]][j,-1])/sqrt(sim_runs) # calculate "standard error" for balance score.
      balance_list[[i]][j,4] <- sd(balance_OR[[i]][j,-1])/sqrt(sim_runs)
    }
  }
  
  
  
  # preparing other components
  confounder <- logical(2L)
  names(confounder) <- c("C1","C2")
  if(!continuous_A){ # indicator for the presence of confounding
    if(IntrA[1] != 1 && IntrA[3] != 1 && IntrY[1] && IntrY[3]){
      confounder[1] <- TRUE
    }
    if(IntrA[2] != 1 && IntrA[4] != 1 && IntrY[2] && IntrY[4]){
      confounder[2] <- TRUE
    }
  } else {
    if(IntrA[1] && IntrA[3] && IntrY[1] && IntrY[3]){
      confounder[1] <- TRUE
    }
    if(IntrA[2] && IntrA[4] && IntrY[2] && IntrY[4]){
      confounder[2] <- TRUE
    }
  }
  interactions <- matrix(NA,nrow = 2,ncol = 4)
  colnames(interactions) <- c("IA-C1A","IA-C2A","IY-C1Y","IY-C2Y")
  rownames(interactions) <- c("Additive","Multiplicative")
  if(!continuous_A){
    interactions[2,1] <- IntrA[3]/IntrA[1] # for binary A the effects are RR and we usually talk about modification on mult. scale.
    interactions[2,2] <- IntrA[4]/IntrA[2]
  } else {
    interactions[1,1] <- IntrA[3]-IntrA[1] # for continuous A the effects are E.
    interactions[1,2] <- IntrA[4]-IntrA[2]
  }
  
  interactions[1,3] <- IntrY[3]-IntrY[1]
  interactions[1,4] <- IntrY[4]-IntrY[2]
  
  causal_structure <- list(
    treatment = "A",
    response = "Y",
    confounder = confounder,
    interactions = interactions,
    noise_tr = ifelse(RR_VA_A!=1 || E_VA_A, "VA",NA),
    noise_resp = ifelse(E_VY_Y,"VY",NA)
  )
  
  returnlist <- list(bias_table = results, 
                     balance = balance_list,
                     causal_structure = causal_structure,
                     call = call, 
                     args = arguments,
                     formulas = formula_list,
                     simdata = list(IPW = raw_results_IPW, OR = raw_results_OR, DR = raw_results_DR, raw = buildDataSim()),
                     rawdatadescription = Hmisc::describe(buildDataSim())
                     )
  class(returnlist) <- "BiasSimulation"
  attr(returnlist,"simulation") <- c("runs" = sim_runs,"size" = size)
  attr(returnlist,"type_A") <- ifelse(continuous_A, "continuous", "binary")
  attr(returnlist$balance,"abs_scale") <- bal_abs_scale
  return(returnlist)
}



summary.BiasSimulation <- function(x){ # summary method for class "BiasSimulation"
  stopifnot(class(x)=="BiasSimulation")
  summary_obj <- list(bias_table = x$bias_table, balance = x$balance, causal_structure = x$causal_structure, 
                      formulas = x$formulas, attributes = attributes(x))
  attr(summary_obj,"true_effect") <- x$args$E_A_Y
  class(summary_obj) <- "summary.BiasSimulation"
  return(summary_obj)
}
  
print.summary.BiasSimulation <- function(x, printfull = FALSE){ # print method for summary.BiasSimulation
  cat("Simulation of a causal structure with confounding interaction:\n")
  cat("Simulation runs", x$attributes$simulation[1], "times; Sample size per run:", x$attributes$simulation[2],"\n")
  cat("Treatment:", x$attributes$type_A, x$causal_structure$treatment[1], "; Response:",x$causal_structure$response[1],"\n")
  cat("True causal effects:",attr(x,"true_effect"),"\n")
  cat("Interaction strengths:\n")
  print(x$causal_structure$interactions)
  cat("\n")
  cat("Formula used: \n")
  for(i in 1:length(x$formulas)){
    cat(as.character(x$formulas[i][[1]]),";")
  }
  cat("\n--------\n")
  cat("Bias table: \n")
  defaultdigits <- getOption("digits")
  options(digits = 4)
  print(x$bias_table)
  cat("\n")
  cat("--------\n")
  cat("Balance check results:\n")
  for(i in 1:length(x$balance)){
    cat("Formula",names(x$balance)[i],":\n")
    if(!printfull){
      print(x$balance[[i]][1:2,])
    } else {
      print(x$balance[[i]])
    }
    cat("\n")
  }
  cat("Note: for full balance table please use `print(summary(res), printfull = TRUE)`.\n")
  options(digits = defaultdigits)
}


print.BiasSimulation <- function(x){ # print method for class BiasSimulation
  print(x$bias_table)
  cat("\n--------\n")
  cat("Don't forget to save your simulation results the next time!")
}

describeBiasSimulation <- function(x){ # describe the simulation datasets.
  x$rawdatadescription
}

makeBiasTable <- function(x, ...){ # primarily written to combine biases of IPW estimates from different simulation runs together. not used now.
  datarows0 <- c(x$args$RR_VA_A,x$args$IntrA)
  returndata <- data.frame()
  
  for(i in seq_len(NROW(x$bias_table))){
    datarows <- c(datarows0,x$bias_table[i,1],x$bias_table[i,2],x$bias_table[i,2]-x$bias_table[i,3],x$bias_table[i,2]+x$bias_table[i,3])
    returndata <- rbind.data.frame(returndata,datarows)
  }
  colnames(returndata) <- c("RR_V_A","IntrA1","IntrA2","IntrA3","IntrA4","Name","Bias_IPW","low","high")
  return(returndata)
}


# other aux. functions

inflInteractionStrength <- function(Intr, step, multiplic = TRUE){ # define a function to update Intr
  if(multiplic){
    Intr <- log(Intr)
  }
  direction_C1 <- ifelse((Intr[3] - Intr[1]) >= 0,1,-1)
  direction_C2 <- ifelse((Intr[4] - Intr[2]) >= 0,1,-1)
  
  Intr[3] <- Intr[3] + 0.5 * step * direction_C1
  Intr[1] <- Intr[1] - 0.5 * step * direction_C1
  Intr[4] <- Intr[4] + 0.5 * step * direction_C2
  Intr[2] <- Intr[2] - 0.5 * step * direction_C2
  
  if(multiplic){
    Intr <- exp(Intr)
  }
  return(Intr)
}
