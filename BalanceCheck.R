# Dependency: cobalt, lmw
balanceCheck <- function(W, data, abs_scale = TRUE, continuous_A = FALSE, ...){ # Check IPW balance between treated and untreated group. Binary version.
  # abs_scale: present the absolute value of the distances.
  
  # preparing covariates to be checked.
  cov <- A ~ C1 + C2 + VA + VY + IA + IY # Interaction and noises added
  addl <- ~ C1:IA + C2:IA + C1:IY + C2:IY # checking interaction term balance (HOW TO INTERPRET?)
  if(E_VA_A == 0 && RR_VA_A == 1){
    update(cov, .~.-VA)
  }
  if(E_VY_Y == 0){
    update(cov, .~.-VY)
  }
  if(IntrA[1]==IntrA[3] && IntrA[2]==IntrA[4] && E_IA_A == 0 && RR_IA_A == 1){ #If Interaction betw IA and Effect C1_A/C2_A not present
    update(cov, .~.-IA)
    update(addl, .~. - C1:IA - C2:IA)
  }
  if(IntrY[1]==IntrY[3] && IntrY[2]==IntrY[4] && E_IY_Y == 0){ # If Interaction betw IY and Effect C1_Y/C2_Y not present
    update(cov, .~.-IY)
    update(addl, .~. - C1:IY - C2:IY)
  }
  
  # calculate standardized mean differences
  if(continuous_A == FALSE){
    baltab_SMD <- cobalt::bal.tab(x = cov, data = data, weights = W, int = FALSE, s.d.denom = "pooled",
                                  addl = addl, continuous = "std", binary = "std", abs = abs_scale)
    x <- data.frame(varname = rownames(baltab_SMD$Balance), balance = baltab_SMD$Balance$Diff.Adj)
  } else {
    baltab_Corr <- cobalt::bal.tab(x = cov, data = data, weights = W, int = FALSE, s.d.denom = "weighted",
                                   addl = addl, continuous = "std",binary = "std",abs = abs_scale)
    x <- data.frame(varname = rownames(baltab_Corr$Balance), balance = baltab_Corr$Balance$Corr.Adj)
  }
  
  
  attr(x,"continuous") <- continuous_A 
  return(x) # Output: A df with standardized diff. measures.
}

weightImputationOR <- function(formula, data, ...){
  formula <- update(formula,  . ~ . + A)
  LmwBalance <- lmw::lmw(formula = formula, treat = "A",data = data,
                         estimand = "ATE", method = "URI")
  return(LmwBalance$weights)
}