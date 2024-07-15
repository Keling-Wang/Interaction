# Library packages
library(rlang)
library(pryr)
library(truncnorm)
library(Hmisc)
library(cobalt)
library(lmw)
library(dbarts)
library(SuperLearner)
library(tmle)
library(dplyr)
library(ggplot2)
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
# Assign parameters. Note: be careful with C1/C2/V distribution assignment and effects!
#  Sometimes you will result in an A such that marginal Pr(A) is too small or large and
#  followed by severe positivity violation in your simulation dataset and NAs produced.

# Runs and sample size per run.
sim_runs = 500L
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
distr_A = "rnorm"
A_prob = 0.25
min_A = -Inf
max_A = Inf
E_A_Y = 3.00

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
distr_C2 = "rgamma"
center_C2 = 4
var_C2 = 2
min_C2 = NA
max_C2 = NA

# Characteristics of VA (`Noise`)
distr_VA = "rbeta"
center_VA = 3
var_VA = 2.5
min_VA = 0
max_VA = 6
RR_VA_A = 1
E_VA_A = 0.4

# VY
distr_VY = "rlnorm"
center_VY = 3.3
var_VY = 0.25
min_VY = NA
max_VY = NA
E_VY_Y = 0.13

# IA
IA_prob = 0.386
IntrA = c(0.68,2.05,0.68,2.05)
RR_IA_A = 0.89
E_IA_A = 0

# IY
IY_prob = 0.466
IntrY = c(2.22,-0.85,2.22,-0.85)
E_IY_Y = -0.35



# models you want to use: 
# inside the function a formula_list is predefined as:
# c(~ C1+C2, ~ (C1+C2)*I-I, ~ (C1+C2)*I, ~ C1+C2+V), i.e.,
# no interaction; only interaction terms; 3-variable model; mis-adjusted model.
# If you want to use your own models please assign a new formula_list here and modify 
#  the call (add "formula_list = your_formula_list")
formula_list <- c(
  ~ C1 + C2,
  ~ (C1 + C2) * IA - IA,
  ~ (C1 + C2) * IA
)
########################################

# Function too long. Stored as a call.
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
                        continuous_A = FALSE,
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


#1. Inflate RR_C1_A and RR_C2_A by a same percentage. Pr(A|Cx,I=1)/Pr(A|Cx,I=0) remains the same.
#call <- result$call
steplength = 0.25 # Inflate by 25% each loop.
Intr_infl <- Intr <- c(0.35,0.15,0.25,0.185) # start from here
result_list_1 = list()
call$Intr <- expr(Intr_infl)

for(i in 1:8){
  Intr_infl <- Intr * ((1 + steplength)^i)
  res <- eval(call)
  summary(res)
  result_list_1 <- c(result_list_1,list(res))
}


#2. Inflate Pr(A|C,I=0)/Pr(A|C,I=1) (Interaction strength)
steplength = 0.10 # Inflate by 10% each loop
Intr <- c(0.35,0.15,0.25,0.185) #start from here


result_list_2 = list()
Intr_infl_st <- Intr
call$Intr <- expr(Intr_infl_st) # modify call
for(i in 1:8){
  Intr_infl_st <- inflInteractionStrength(Intr_infl_st,steplength)
  cat("Round ",i," ; Pr(A|C1,I=1)/Pr(A|C1,I=0) = ",round(Intr_infl_st[3]/Intr_infl_st[1],3),"\n")
  cat("          ","Pr(A|C2,I=1)/Pr(A|C2,I=0) = ",round(Intr_infl_st[2]/Intr_infl_st[4],3),"\n")
  res <- eval(call)
  summary(res)
  result_list_2 <- c(result_list_2,list(res))
}


#3. inflate RR_V_A (effect of V on A)
steplength_V = 0.25
steplength = 0.25
Intr_infl_st <- Intr <- c(1.45,1.35,1.25,1.185) #Start from Interaction I:C1-A = 0.86, I:C2-A = 0.88
RR_V_A_infl <- RR_V_A <- 1.00
call$Intr <- expr(Intr_infl_st)
call$RR_V_A <- expr(RR_V_A_infl)

result_list_3 <- list()
filename <- paste("result_", as.character(Sys.Date()),as.character(round(as.numeric(Sys.time()))),".txt",sep = "")

for(i in 1:5){
  Intr_infl_st <- inflInteractionStrength(Intr_infl_st,steplength)
  cat("\n\nRound ",i," ; Pr(A|C1,I=1)/Pr(A|C1,I=0) = ",round(Intr_infl_st[3]/Intr_infl_st[1],3),"\n")
  cat("          ","Pr(A|C2,I=1)/Pr(A|C2,I=0) = ",round(Intr_infl_st[4]/Intr_infl_st[2],3),"\n")
  for(j in -7:7){
    RR_V_A_infl <- RR_V_A * (1 + steplength_V)^j
    cat("Round",i,":",j,", RR_V_A = ",RR_V_A_infl)
    res <- eval(call)
    summary(res)
    write(capture.output(summary(res)), file = filename, append = TRUE)
    result_list_3 <- c(result_list_3,list(res))
  }
}

# Make plots
plot_data <- data.frame()
for(i in seq_len(length(result_list_3))){
  dt <- makeBiasTable(result_list_3[[i]])
  plot_data <- rbind.data.frame(plot_data,dt)
}
plot_data[,c(1:5,7:9)] <- as.data.frame(sapply(plot_data[,c(1:5,7:9)],as.numeric))
plot_data$IntrStrength <- plot_data$Intr3/plot_data$Intr1

# Create the line plot with error bars
ggplot(plot_data[plot_data$Name=="C1 + C2" | plot_data$Name=="(C1 + C2) * I" | plot_data$Name=="C1 + C2 + V",], aes(x = log(RR_V_A), y = Bias, colour = -log(IntrStrength), group = interaction(Intr1,Name))) +
  geom_line(aes(linetype = Name)) +
  geom_abline(slope = 0,intercept = 0,color = "red") +
  scale_color_continuous(type = "viridis") +
  #geom_errorbar(aes(ymin = low, ymax = high), width = 0.1) +
  labs(x = "ln(RR_V_A)") + 
  ggtitle(label = "Bias of IPW estimator (Only confounders) w.r.t. different noise and interaction strength") + 
  theme(plot.title = element_text(hjust = 0.5))

