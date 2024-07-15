library(pryr)
library(rlang)
library(truncnorm)
library(Hmisc)
library(cobalt)
library(lmw)

#library(betareg)


library(dbarts)
library(SuperLearner)
library(tmle)
#library(AIPW)

# A
mean_A = 2
var_A = 1
distr_A = "rnorm"
A_prob = 0.25
min_A = -Inf
max_A = Inf
E_A_Y = 2.5

# Y
mean_Y = 1.88
var_Y = 0.8

# C1
distr_C1 = "rtruncnorm"
center_C1 = 2
var_C1 = 3
min_C1 = 0
max_C1 = 10

# C2
distr_C2 = "rlnorm"
center_C2 = 2
var_C2 = 2
min_C2 = NA
max_C2 = NA

# VA
distr_VA = "rlnorm"
center_VA = 2.0
var_VA = 2.5
min_VA = NA
max_VA = NA
RR_VA_A = 1.22
E_VA_A = 0.4

# VY
distr_VY = "rnorm"
center_VY = 3.3
var_VY = 0.25
min_VY = NA
max_VY = NA
E_VY_Y = 0.13

# IA
IA_prob = 0.5
IntrA = c(1.35,2.05,1.5,0.77)
RR_IA_A = 0.89
E_IA_A = 0

# IY
IY_prob = 0.5
IntrY = c(2.22,-0.85,1.22,0.85)
E_IY_Y = 0

# size
size = 100000

#formula
formula_list <- c(
  ~ C1 + C2, 
  ~ (C1+C2)*IA,
  ~ (C1 + C2)* IY
)

source(".\\Estimators.R")
source(".\\BuildDataSim.R")
source(".\\GenerateContiVars.R")
source(".\\GetProbDensA.R")
source(".\\BalanceCheck.R")

data <- buildDataSim()
est1 <- ipwEstimator_binaryA(~ C1 + C2, data = data)
est2 <- ipwEstimator_binaryA(~ (C1 + C2) * IY, data = data)
est3 <- ipwEstimator_binaryA(~ C1 + C2, data = data[data$IY==1,])
est4 <- ipwEstimator_binaryA(~ (C1 + C2) * IY, data = data[data$IY==1,])
est5 <- ipwEstimator_binaryA(~ C1 + C2,data = data[data$IY==0,])
est6 <- ipwEstimator_binaryA(~(C1+C2)*IA,data = data[data$IY==0,])
est7 <- outcomeRegEstimator(~ C1 + C2, data = data)
est8 <- outcomeRegEstimator(~ (C1 + C2) * IY,data = data)
est9 <- outcomeRegEstimator(~ C1 + C2, data = data[data$IY==1,])
print(c(est1$est,est2$est,est3$est,est4$est,est5$est,est6$est,est7,est8,est9)-2.5)

data2 <- buildDataSim(continuous_A = TRUE)
est1_cont <- ipwEstimator_continuousA(~ (C1 + C2)*IA, data = data2,numerCond = ~ 1)
