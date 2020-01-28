################################################################################################################
#
# robustness.R     28 January 2020
#
# Copyright (C) 2020 Naomi Radke
#
#
#
# This script measures the robustness of the Altherr thinning regime under model parameter, climate, wood price,
# carbon tax and discount rate uncertainty. Performance measures are NPV yield and NPV carbon (from objective fct.).
# The uncertainty scenarios have been created using LHS previously.
#
#
#
# robustness.R is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#
# This code is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU Lesser General Public
# License for more details (<http://www.gnu.org/licenses/>).
#
#
# Any questions? Naomi Radke (naomikradke@gmail.com) 
# ===================================================================================================
#
# Input:
#- clim_dat             Climate data at stand location
#         - stand_scen    stand index number
#         - year          projection year (5 year increment)
#         - DI            climate index 1
#         - GDD           climate index 2
#
#- stand_dat            Beech stand data
#         - stand_scen    stand index number
#         - tree          tree number
#         - species       tree species (integer; i.e, 41)
#         - dbh           diameter at breast height
#         - age           age of tree
#         - ni            number of trees per hectare
#
#
#- harvest_dat          Harvesting plan according to a certain thinning strategy
#         - stand_scen    stand index number
#         - year          projection year when harvesting 
#         - age           age of trees harvested during harvest year
#         - htype         type of harvest thinning (low, uniform, or high thinning)
#         - g             thinning value (if present)
#         - g_r           desired basal area
#
#- param_scen           Dataframe of growth model parameter values for every scenario
#         - X1            parameter 1
#         - X2            parameter 2
#         - ...            
#         - X22           parameter 22
#
#- clim_scen            Dataframe of drought index development over considered time period for every scenario
#         - scenario1     drought index for every 5-year step ]0,1] for scenario 1
#         - scenario2     drought index for every 5-year step ]0,1] for scenario 2
#         - ...
#         - scenario20    drought index for every 5-year step ]0,1] for scenario 20
#
#- pr_scen              Dataframe of timber price index development over considered time period for every scenario
#         - scenario1     price index for every 5-year step [%, where price in year 2010 = 100%] for scenario 1
#         - scenario2     price index for every 5-year step [%, where price in year 2010 = 100%] for scenario 2
#         - ...
#         - scenario10    price index for every 5-year step [%, where price in year 2010 = 100%] for scenario 10
#
#- dr_scen              Dataframe of discount rate development over considered time period for every scenario
#         - scenario1     discount rate for every 5-year step [-1,1] for scenario 1
#         - scenario2     discount rate for every 5-year step [-1,1] for scenario 2
#         - ...
#         - scenario10    discount rate for every 5-year step [-1,1] for scenario 10
#
# Output:
#- scen_obj              Dataframe of the objective values (NPV of timber yield and net carbon sequestration) for each scenario
#         - model         number of model scenario [1:84]
#         - DI            number of drought index scenario [1:20]
#         - DR            number of discount rate scenario [1:10]
#         - WP            number of timber price index scenario [1:10]
#         - NPV_yield     Net Present Value of timber yield from all harvest events in the considered time period
#         - Carb_seq      Net carbon sequestration in above-ground tree biomass and wood products in the considered time period
#         - robust        TRUE/FALSE minimum performance requirements for the two objectives are met
#
#         
################################################################################################################

#___PREPARATIONS___#

# Load required packages
  
  library(foreach) # for parallelizing
  library(doParallel) # for parallelizing
  library(ggplot2) # for density plots
  library(patchwork) # for multiple plotting
  library(plyr) # for counting TRUE and FALSE in robustness

# Load required data

  # initialization data to run forest growth model
  stand_dat <- read.csv(".csv", header=TRUE, na.strings='(null)',sep=";", dec=".")
  clim_dat <- read.csv(".csv", header=TRUE, na.strings='(null)',sep=";", dec=".")
  harvest_dat <- read.csv(".csv", header=TRUE, na.strings='(null)',sep=";", dec=".") # Altherr thinning
  
  # load input tables with scenarios
  param_scen <- read.csv("data/scenarios/params_acc_scen_disturbed.csv", row.names = 1)[,1:22] # acceptable model param scenarios
  clim_scen <- read.csv("data/scenarios/DI_scenarios.csv", row.names = 1) # nrow = nyears, ncol = 20 climate scenarios
  pr_scen <- read.csv("data/scenarios/PI_scen.csv", row.names = 1) # nrow = nyears, ncol = 10 scenarios price index
  dr_scen <- read.csv("data/scenarios/dr_scen_arma.csv", row.names = 1) # dr simulations with fitted ARMA

  
# Define vectors that contain the scenarios for each uncertainty
scen_model_params <- seq(1:nrow(param_scen))  # number of model parameter scenarios
scen_DI <- seq(1:ncol(clim_scen))             # number of drought index (climate) scenarios
scen_DR <- seq(1:ncol(dr_scen))               # number of discount rate scenarios
scen_WP <- seq(1:ncol(pr_scen))               # number of wood price scenarios


# Create scenarios through total enumeration of individual scenarios
scen <- expand.grid("model" = scen_model_params, "DI" = scen_DI, "DR" = scen_DR, "WP" = scen_WP)#, "CT" = scen_CT)

# Source the required functions
source("objective_function_parallel.R")  # source the objective function that calculates objective values for NPV timber and carbon

#(the following functions are part of the growth model and are not provided here)
source("growth_yield_stand_model_final.R") # source forest growth model
source("dominate_class_function.R") # source dominant class function needed for forest model

#___ROBUSTNESS MEASURE___#

# Define minimum performance requirements
min_NPV_yield = 0 # for NPV of timber yield
min_C_seq = 0     # for net carbon sequestration

# Run objective function for each scenario
obj_scen = obj.fun(scenarios = scen)

# Merge scenario data frame with objectives dataframe
scen_obj <- cbind(scen, obj_scen)
  
# Add a new column with TRUE/FALSE if min. requirements are fulfilled
  # define test function
  test_funct <- function(x){all(x["NPV_yield"] >= min_NPV_yield & x["Carb_seq"] >= min_C_seq)}

  # apply function and save test results in new column
  scen_obj$robust <- apply(scen_obj, 1, test_funct)
  
  # count robustness
  rob_perc <- count(scen_obj, vars = "robust")/nrow(scen_obj)*100 # % scenarios that fulfill min requirements = robustness index
  
  
  
  # write csv file
  currentDate <- Sys.Date()
  csvFileName <- paste("output/evaluation/robustness",currentDate,".csv",sep="")
  write.csv(scen_obj, file = csvFileName)
  
  
