buildDataSim <- function(continuous_A = FALSE, ...){
  # continuous A and ,,, and different causal structures!!
  # output: a data.frame with additional class and attributes stating causal structures.
  # current variables to be simulated: C1, C2, VA, VY, A**, IA, IY. In some structures VY and IY will be set to zero.
  data_sim <- data.frame( # generate C1, C2, VA, VY
    C1 = generateContiVars(distr_C1, center_C1, var_C1, min_C1, max_C1, size),
    C2 = generateContiVars(distr_C2, center_C2, var_C2, min_C2, max_C2, size),
    VA = generateContiVars(distr_VA, center_VA, var_VA, min_VA, max_VA, size),
    VY = generateContiVars(distr_VY, center_VY, var_VY, min_VY, max_VY, size)
  )
  
  # generate IA and IY (variable interacting with C1/C2-A and C1/C2-Y effect)
  data_sim$IA <- rbinom(size,1,IA_prob)
  data_sim$IY <- rbinom(size,1,IY_prob)
  
  # yield A.
  if(continuous_A == FALSE){ # binary A data.
    # Convert RR to OR
    OR_IntrA <- IntrA * (1 - A_prob) / (1 - IntrA * A_prob)
    OR_VA_A <- RR_VA_A * (1 - A_prob) / (1 - RR_VA_A * A_prob)
    OR_IA_A <- RR_IA_A * (1 - A_prob) / (1 - RR_IA_A * A_prob)
    # Calculate Log odds. IntrA[1,2] are for IA=0, IntrA[3,4] are for IA = 1
    LogOdds_A <- with(data_sim, log(A_prob / (1 - A_prob)) + 
                        (log(OR_IntrA[1])*C1 + log(OR_IntrA[2])*C2) * (IA==0) +
                        (log(OR_IntrA[3])*C1 + log(OR_IntrA[4])*C2) * (IA==1) +
                        log(OR_VA_A)*VA + log(OR_IA_A)*IA
    )
    
    post_prob_A <- exp(LogOdds_A)/(exp(LogOdds_A)+1) # transform back 
    data_sim$A <- rbinom(size, 1, post_prob_A)
  } else if(continuous_A == TRUE) {# Continuous A.
    #IntrA : c(E_C1_A|IA=0,E_C2_A|IA=0,E_C1_A|IA=1,E_C2_A|IA=1)
    mu_A <- with(data_sim,
                 mean_A + (IntrA[1]*C1 + IntrA[2]*C2) * (IA==0) +
                          (IntrA[3]*C1 + IntrA[4]*C2) * (IA==1) + 
                          E_VA_A * VA + E_IA_A * IA
    )
    data_sim$A <- switch(distr_A, # Currently only norm, gamma, and chisq for A are supported.
                         "rnorm" = rnorm(size,mean = mu_A, sd = sqrt(var_A)),
                         "rlnorm" = {
                           rlnormgen <- function(){
                            sdlog = sqrt(log(var_A/mu_A^2+1))
                            meanlog = log(mu_A) - 0.5*var_A^2
                            rlnorm(size,meanlog = meanlog, sdlog = sdlog)
                           }
                           rlnormgen()
                         },
                         "rgamma" = rgamma(size, shape = (mu_A^2) / var_A, scale = var_A / mu_A),
                         "rchisq" = rchisq(size, df = mu_A),
                         stop("Specified distribution for A not supported yet")
                         )
  }
  
  
  
  # yield Y
  mu_Y <- with(data_sim,
                 mean_Y + E_A_Y * A + E_VY_Y * VY + E_IY_Y * IY + 
                 (IntrY[1]*C1 + IntrY[2]*C2) * (IY==0) +
                 (IntrY[3]*C1 + IntrY[4]*C2) * (IY==1)
              )
  data_sim$Y <- rnorm(size,mean = mu_Y,sd = sqrt(var_Y)) # to avoid bias from model misspecification Y can only be normally distributed currently.
  
  # preparing output
  class(data_sim) <- c(class(data_sim),"simdata")

  return(data_sim)
}