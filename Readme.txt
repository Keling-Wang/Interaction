Files: 
 SimulationInteractionIPW.R: main function. contains main function and a summary() method for the output.
 __main__.R: interface. You can assign parameters and run simulations in this R script.

Mainly used formulars: 
 Simulation of A: 
	OR(C1 or C2 or Intr) = RR * (1 - A_prob0) / (1 - RR * A_prob0)
	logit(A) = logit(A_0) + sum(log(OR)s)
	A = rbinom(size, 1, expit(logit(A))
	
 Simulation of Y: 
	mean(Y) = mean_Y0 + E_A_Y * A + E_C1_Y * C1 + E_C2_Y * C2
	Y = rnorm(size, mean(Y), sqrt(var_Y))