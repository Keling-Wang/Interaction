source(".\\SimulationInteractionIPW.R")
call <- expr(simulationInteraction(sim_runs = 100L, size = 15000L,
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

res <- eval(call)
