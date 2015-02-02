# phyto-division-rate-model
Matrix model to estimate daily divsion rate of a phytoplankton population from diel change in cell size distributions 

Scripts are expecting user supplied inputs:
  Predefined cell size classes that cover range of cell volumes for given phytoplankton population
  Hourly observations of cell size distributions (counts of how many cells fall into each size bin for each hour)
  Radiation data that corresponds to the time frame for population observations
 
Transition matrix construction is performed by matrix_const.m
negloglike_calc.m calculates the negative log likelihood of a set of parameters given a day of observations assuming a Dirichlet-mulitnomial distirbution and a two-subpopulation model construction.
loglike_opt.m is a sample script that sets up an optimization using MATLAB's fmincon to find the set of parameters that maximize the likelihood function given a day of observations. 
sample_dat.mat contains simulated sample data for a Synechococccus population

Prerequisites
MATLAB

Full description of the model can be found in Hunter-Cevera et al. 2014 (PNAS), but also see Hunter-Cevera 2014 (MIT/WHOI thesis)

Notes on optimization: We use a strategy of many mulitple start points for the optimization to avoid finding only local minima. With 13 parameters in the two-component model version, we find sometimes (depending on the data) that there can be several combinations of parameters that result in very close likelihood values. 

