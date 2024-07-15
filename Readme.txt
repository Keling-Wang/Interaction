Files: 
 __main__.R: interface. You can assign parameters and run simulations in this R script.
 BalanceCheck.R: calculate balance scores
   function: balanceCheck, weightImputationOR
 BuildDataSim.R: function for generating the simulation dataset. 
   function: buildDataSim
 Estimators.R: IPW, OR, and DR estimators
   function: ipwEstimator_binaryA, ipwEstimator_continuousA, outcomeRegEstimator, drEstimator_binaryA
 GenerateContiVars.R: generate continuous variables using given parameters and distributions.
   function: generateContiVars
 GetProbDensA.R: get probability density of treatment for IPW continuous A.
   function: getProbDensA
 SimulationInteractionIPW.R: main function. contains main function and a summary() method for the output.
   functions: simulationInteractionIPW, summary.BiasSimulation, print.summary.BiasSimulation, ...

Function call tree in simulationInteractionIPW: 

simulationInteractionIPW # main function
  |
  buildDataSim # build simulation dataset
  |  |
  |  generateContiVars # generate continuous variables
  |
  ipwEstimator_binaryA
  |
  outcomeRegEstimator
  |
  drEstimator_binaryA
  |  |
  |  tmle::tmle
  |
  balanceCheck
  |  |
  |  cobalt::bal_tab
  |
  weightImputationOR # impute (implied) weight for outcome regression
  |  |
  |  lmw::lmw
  |
  ipwEstimator_continuousA
  |  |
  |  getProbDensA # get probability density for A
  |
  drEstimator_continuousA
    |
    NA (not implemented)