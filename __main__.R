# Library packages
library(rlang)
library(pryr)
library(truncnorm)
library(Hmisc)
library(cobalt)
library(lmw)
library(dbarts) # sometimes I need to manually load this package for tmle
library(SuperLearner)
library(tmle)
library(dplyr)
#library(ggplot2)

# Import function
# setwd(...)
source(".\\Estimators.R")
source(".\\BuildDataSim.R")
source(".\\GenerateContiVars.R")
source(".\\GetProbDensA.R")
source(".\\BalanceCheck.R")
source(".\\SimulationInteractionIPW.R") # main function and its summary() method.

## If data loaded, please use summary() for each result and each sub-list in result_list
## to see results.

###############################
# Assign parameters. 
#  Sometimes you will result in an A such that marginal Pr(A) is too small or large and NAs are produced.

# Runs and sample size per run.
sim_runs = 100L
size = 5000L

# Current available distribution options and parameters needed: 
# "rnorm","rtruncnorm","runif","rbeta","rlnorm","rgamma","rchisq"
# "rnorm": param: center, var
# "rtruncnorm": param: center, var, min, max
# "runif": param: min, max
# "rbeta": param: center (approx.), var, min, max
# "rlnorm": param: center, var
# "rgamma": param: center, var
# "rchisq": param: center

# A
mean_A = 2
var_A = 1
distr_A = "rnorm" # for A only normal, truncated normal, and gamma are available for now.
A_prob = 0.25
min_A = -Inf
max_A = Inf
E_A_Y = 3.00 
continuous_A = FALSE # if set continuous_A = FALSE, only A_prob will be used, and other parameters can be NAs.

# Y
mean_Y = 3.88
var_Y = 2.46

# Characteristics of C1
distr_C1 = "rtruncnorm"
center_C1 = 2
var_C1 = 3
min_C1 = 0
max_C1 = 10

# Characteristics of C2
distr_C2 = "runif"
center_C2 = NA
var_C2 = NA
min_C2 = 1.5
max_C2 = 4.5

# Characteristics of VA (`Noise`)
distr_VA = "rbeta"
center_VA = 3
var_VA = 2.5
min_VA = 0
max_VA = 6
RR_VA_A = 1
E_VA_A = 0

# VY
distr_VY = "rlnorm"
center_VY = 3.3
var_VY = 0.25
min_VY = NA
max_VY = NA
E_VY_Y = 0.13

# IA
IA_prob = 0.386
IntrA = c(0.68,2.05,0.68,2.05) # IntrA[1] = RR_C1_A for IA=0, IntrA[2] = RR_C2_A for IA=0; [3]and[4] are for IA = 1. In continuous case there're E_C1_A, ...
RR_IA_A = 1
E_IA_A = 0

# IY
IY_prob = 0.466
IntrY = c(2.22,-0.85,2.22,-0.85) # same as IntrA. c(E_C1_Y|IY=0, E_C2_Y|IY=0, ...)
E_IY_Y = 0



# models you want to use: 
# inside the function a formula_list is predefined as:
# c(~ C1+C2, ~ (C1+C2)*IA-IA, ~ (C1+C2)*IA, ~ C1+C2+V), i.e.,
# no interaction; only interaction terms; 3-variable model; mis-adjusted model.
# If you want to use your own models please assign a new formula_list here and modify 
#  the call (add "formula_list = your_formula_list")
formula_list <- c(
  ~ C1 + C2,
  ~ (C1 + C2) * IA - IA,
  ~ (C1 + C2) * IA
)
########################################

# 
call <- expr(
  simulationInteraction(sim_runs = sim_runs, size = size,
                        distr_C1 = distr_C1,
                        center_C1 = center_C1, var_C1 = var_C1, min_C1 = min_C1, max_C1 = max_C1,
                        distr_C2 = distr_C2, center_C2 = center_C2, var_C2 = var_C2, min_C2 = min_C2, max_C2 = max_C2,
                        distr_VA = distr_VA, center_VA = center_VA, var_VA = var_VA, min_VA = min_VA, max_VA = max_VA,
                        distr_VY = distr_VY, center_VY = center_VY, var_VY = var_VY, min_VY = min_VY, max_VY = max_VY,
                        IntrA = IntrA, RR_IA_A = RR_IA_A, E_IA_A = E_IA_A,#RR_C_A or E_C_A per 1 unit C in different strata of IA
                        IntrY = IntrY, E_IY_Y = E_IY_Y,#E_C_A in different strata of IY
                        IA_prob = IA_prob, IY_prob = IY_prob,
                        mean_Y = mean_Y, var_Y = var_Y,
                        A_prob = A_prob, distr_A = distr_A, mean_A = mean_A, var_A = var_A, min_A = min_A, max_A = max_A,
                        E_A_Y = E_A_Y, RR_VA_A = RR_VA_A, E_VA_A = E_VA_A, E_VY_Y = E_VY_Y,
                        continuous_A = continuous_A,
                        numerCond = ~ 1,
                        bal_abs_scale = FALSE,
                        use_DRestimator = TRUE,
                        formula_list = formula_list
  ) 
)



#0. Evaluate the call once.
result <- eval(call)
# Display the result.
summary(result)