######################################################################################################
#
# lGz_plot_lapply.R     28 January 2020
#
# Copyright (C) 2020 Naomi Radke
#
#
# 
# This script calculates the lGz (yearly total volume growth) for every model period and for every parameter scenario
# of the adapted forest growth model (original by Trasobares et al. (2016)). 
# It compares the observed lGz with the one produce by the model scenarios in order to find an acceptably fitting subset
# for use as scenarios that represent model parameter uncertainty.
#
#
#
# lGz_plot_lapply.R is free software: you can redistribute it and/or modify
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
#         - DI            stand drought index
#         - GDD           sum of degree-days (April-October)
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
#- harvest_dat          Harvesting plan according to the Altherr thinning regime
#         - stand_scen    stand index number
#         - year          projection year when harvesting 
#         - age           age of trees harvested during harvest year
#         - htype         type of harvest thinning (low, uniform, or high thinning)
#         - g             thinning value (if present)
#         - g_r           desired basal area
#
#- obs_harvest_dat      Observed harvest data for individual trees for each harvest event
#         - year          projection year when harvesting
#         - dbh           diameter at breast height (cm)
#         - height        height (m)
#         - n_harv        number of trees harvested at a certain dbh
#         - f             quotient for calculating tree height based on diameter and age
#         - age           age of tree (years)
#
#- obs_stand_aft        Observed stand data after each harvest event for individual trees
#         - year          projection year when harvesting
#         - age           age of tree (years)
#         - dbh           diameter at breast height (cm)
#         - height        height (m)
#         - n             number of trees standing at a certain dbh
#         
#- params               Matrix with nrow beech growth model parameter scenarios (disturbed optimal parameter values)
#         - X1            parameter 1
#         - X2            parameter 2
#         - ...            
#         - Xx           parameter x
#
#- DEoptim              Matrix of the best parameter set found by optimization
#         - V1            parameter 1
#         - V2            parameter 2
#         - ...            
#         - Vx           parameter x
#
#
# Output:
#- lGz_comparison      Dataframe of the stand's yearly volume growth at different stand ages for all params scenarios
#         - lGz         yearly stand volume growth at a specific stand age (m3/ha/yr)
#         - age         stand age (years)
#         - source      numbered model parameter scenario, best fit scenario or original scenario
#         - Sources     "observation" or "model parameter examples"
#
#- params_acc         Dataframe of the parameter scenarios and values of each growth model parameter
#         - X1            parameter 1
#         - X2            parameter 2
#         - ...            
#         - Xx           parameter x
#
#- lGz_comparison_acc  Dataframe of the stand's yearly volume growth at different stand ages for acceptable model scenarios, observation, original or best fit
#         - lGz         yearly stand volume growth at a specific stand age (m3/ha/yr)
#         - age         stand age (years)
#         - source      numbered model parameter scenario, best fit scenario or original scenario
#         - Sources     "observation" or "acceptable model parameter examples"
#               
############################################################################################################################

# clear work space
# clear work space
rm(list = ls())
graphics.off()

setwd("~/...")

# Load packages
library(ggplot2) # for plotting
library(dplyr) # to stack dataframes for plotting
library(GGally) # extension for ggplot to make the correlation pairs plot


# Load input data for model and observed data (not provided but can be exchanged with another forest growth model)
stand_dat <- read.csv(".csv", header=TRUE, na.strings='(null)',sep=";", dec=".")
clim_dat <- read.csv(".csv", header=TRUE, na.strings='(null)',sep=";", dec=".")
harvest_dat <- read.csv(".csv", header=TRUE, na.strings='(null)',sep=";", dec=".")
# observation
obs_harvest_dat <- read.csv(".csv", header=TRUE, na.strings='(null)',sep=";", dec=".") # trees harvested for each harvest year
obs_stand_aft <- read.csv(".csv", header=TRUE, na.strings='(null)',sep=";", dec=".") # stand data after harvest

# Load the data frame with parameter value scenarios
params <- read.csv("output/precalibration/params_DEoptim_disturbed.csv", row.names = 1)

# Load the best-fit parameter set created by DEoptim 
load("output/precalibration/DEoptim_out_5SE2019-11-29.RData") # best fit of original +-5SE by Trasobares et al 2016
DEoptim <- out_DEoptim$optim$bestmem

# Both the following forest growth model functions are not provided here
# Source the function to calculate the dominate class:
source("growth_model/dominate_class_function.R")

# Source the growth yield model (calibrated and new management module)
source("growth_model/growth_yield_stand_model_final.R")

# ----- lGZ MODEL DATA ---------------------------

# If the below function was run before (e.g. on a cluster), load the model_runs data and skip to the "best fit" section)
load("output/precalibration/model_runs.RData")

# Run the model for each parameter set and output a list of dataframes, each containing the lgZ of the model run
n.samples <- nrow(params)
model_runs <- lapply(c(1:n.samples), function(i){
  p <- as.numeric(params[i,]) # pick a parameter set
  
  model = forest.model(p, clim_dat, stand_dat, harvest_dat) # run the model
  
  # Separate the model output
  stand_scen_gp_table = model$stand
  harvest_stand_scen = model$harvest
  colnames(stand_scen_gp_table) <- c("Stand_scen","Year","Age","N","BA_stand","Dg","Ddom","Hdom","Volume_stand")
  colnames(harvest_stand_scen) <- c("Stand_scen","Year","Age","N","BA_harvested","VOL_h")
  
  # Subtracting "stand volume after thinning" at time = t-1 from "stand volume before next thinning" at time = t
  lGz_model_5 <- (stand_scen_gp_table[2,9] - stand_scen_gp_table[1,9])/5 # lgZ at age 63 (in 1969)
  lGz_model_10 <- (stand_scen_gp_table[4,9] - stand_scen_gp_table[3,9])/5
  lGz_model_15 <- (stand_scen_gp_table[6,9] - stand_scen_gp_table[5,9])/5
  lGz_model_20 <- (stand_scen_gp_table[8,9] - stand_scen_gp_table[7,9])/5
  lGz_model_25 <- (stand_scen_gp_table[10,9] - stand_scen_gp_table[9,9])/5
  lGz_model_30 <- (stand_scen_gp_table[12,9] - stand_scen_gp_table[11,9])/5
  lGz_model_35 <- (stand_scen_gp_table[14,9] - stand_scen_gp_table[13,9])/5
  lGz_model_40 <- (stand_scen_gp_table[16,9] - stand_scen_gp_table[15,9])/5
  lGz_model_45 <- (stand_scen_gp_table[18,9] - stand_scen_gp_table[17,9])/5
  
  lGz_model_5_45 <- c(lGz_model_5, lGz_model_10, lGz_model_15, lGz_model_20, lGz_model_25, lGz_model_30, lGz_model_35, lGz_model_40, lGz_model_45)
  age_model <- c(63, 68, 73, 78, 83, 88, 93, 98, 103)
  lGz_model <- data.frame(lGz= lGz_model_5_45, age= age_model, source = c(rep(paste(i, "model", sep = ""), length(age_model))))
  
})

# Run the model for the DEoptim parameter set and calculate the lGz, stored in a dataframe
best_fit = forest.model(p = as.numeric(DEoptim), clim_dat, stand_dat, harvest_dat)

# Separate the model output
stand_scen_gp_table = best_fit$stand
harvest_stand_scen = best_fit$harvest
colnames(stand_scen_gp_table) <- c("Stand_scen","Year","Age","N","BA_stand","Dg","Ddom","Hdom","Volume_stand")
colnames(harvest_stand_scen) <- c("Stand_scen","Year","Age","N","BA_harvested","VOL_h")

# Subtracting "stand volume after thinning" at time = t-1 from "stand volume before next thinning" at time = t
lGz_best_fit_5 <- (stand_scen_gp_table[2,9] - stand_scen_gp_table[1,9])/5 # lgZ at age 63 (in 1969)
lGz_best_fit_10 <- (stand_scen_gp_table[4,9] - stand_scen_gp_table[3,9])/5
lGz_best_fit_15 <- (stand_scen_gp_table[6,9] - stand_scen_gp_table[5,9])/5
lGz_best_fit_20 <- (stand_scen_gp_table[8,9] - stand_scen_gp_table[7,9])/5
lGz_best_fit_25 <- (stand_scen_gp_table[10,9] - stand_scen_gp_table[9,9])/5
lGz_best_fit_30 <- (stand_scen_gp_table[12,9] - stand_scen_gp_table[11,9])/5
lGz_best_fit_35 <- (stand_scen_gp_table[14,9] - stand_scen_gp_table[13,9])/5
lGz_best_fit_40 <- (stand_scen_gp_table[16,9] - stand_scen_gp_table[15,9])/5
lGz_best_fit_45 <- (stand_scen_gp_table[18,9] - stand_scen_gp_table[17,9])/5

lGz_best_fit_5_45 <- c(lGz_best_fit_5, lGz_best_fit_10, lGz_best_fit_15, lGz_best_fit_20, lGz_best_fit_25, lGz_best_fit_30, lGz_best_fit_35, lGz_best_fit_40, lGz_best_fit_45)
age_best_fit <- c(63, 68, 73, 78, 83, 88, 93, 98, 103)
lGz_best_fit <- data.frame(lGz= lGz_best_fit_5_45, age= age_best_fit, source = "model best_fit", Sources = "model best fit")

# Run the model with the original parameter set (Trasobares et al. 2016) and calculate the lGz, stored in a dataframe
p_orig <- c(-3.4731, -0.00081, 0.8485, -0.09665, -1.1036, -0.04753, 3.3119, 1.0857, # dbh growth
            2.044, 3.733, -0.306, 2.484, 4.164, -0.0372, # mortality
            -8.6828, -0.00227, 1.0648, -0.05216, -0.1545, 0.3371, 3.1526, 2.2973) # height growth)

original = forest.model(p_orig, clim_dat, stand_dat, harvest_dat)

# Separate the model output
stand_scen_gp_table = original$stand
harvest_stand_scen = original$harvest
colnames(stand_scen_gp_table) <- c("Stand_scen","Year","Age","N","BA_stand","Dg","Ddom","Hdom","Volume_stand")
colnames(harvest_stand_scen) <- c("Stand_scen","Year","Age","N","BA_harvested","VOL_h")

# Subtracting "stand volume after thinning" at time = t-1 from "stand volume before next thinning" at time = t
lGz_orig_5 <- (stand_scen_gp_table[2,9] - stand_scen_gp_table[1,9])/5 # lgZ at age 63 (in 1969)
lGz_orig_10 <- (stand_scen_gp_table[4,9] - stand_scen_gp_table[3,9])/5
lGz_orig_15 <- (stand_scen_gp_table[6,9] - stand_scen_gp_table[5,9])/5
lGz_orig_20 <- (stand_scen_gp_table[8,9] - stand_scen_gp_table[7,9])/5
lGz_orig_25 <- (stand_scen_gp_table[10,9] - stand_scen_gp_table[9,9])/5
lGz_orig_30 <- (stand_scen_gp_table[12,9] - stand_scen_gp_table[11,9])/5
lGz_orig_35 <- (stand_scen_gp_table[14,9] - stand_scen_gp_table[13,9])/5
lGz_orig_40 <- (stand_scen_gp_table[16,9] - stand_scen_gp_table[15,9])/5
lGz_orig_45 <- (stand_scen_gp_table[18,9] - stand_scen_gp_table[17,9])/5

lGz_orig_5_45 <- c(lGz_orig_5, lGz_orig_10, lGz_orig_15, lGz_orig_20, lGz_orig_25, lGz_orig_30, lGz_orig_35, lGz_orig_40, lGz_orig_45)
age_orig <- c(63, 68, 73, 78, 83, 88, 93, 98, 103)
lGz_orig <- data.frame(lGz= lGz_orig_5_45, age= age_orig, source = "model original parameters", Sources = "model original parameters")



# ----- lGZ OBSERVED DATA --------------------------


# First, calculate volume of each tree in order to lateron calculate stand volume

# Trees standing after harvest
obs_stand_aft$volume <- (pi/4)*(obs_stand_aft$dbh/100)^2*obs_stand_aft$height*0.5*obs_stand_aft$n

# Trees harvested
obs_harvest_dat$volume <- (pi/4)*(obs_harvest_dat$dbh/100)^2*obs_harvest_dat$height*0.5*obs_harvest_dat$n_harv

# Calculate stand volume after harvest and volume harvested of each observation year
obs_aft_volume <- aggregate(volume~year, data = obs_stand_aft, sum)
obs_harvest_volume <- aggregate(volume~year, data = obs_harvest_dat, sum)

# Calulate harvest volume of each observation year (some observation periods are 3 or 6 years apart, most are 5 years)

lGz_obs_3 <- (obs_aft_volume[2,2]+obs_harvest_volume[1,2]-obs_aft_volume[1,2])/3
lGz_obs_8 <- (obs_aft_volume[3,2]+obs_harvest_volume[2,2]-obs_aft_volume[2,2])/5
lGz_obs_13 <- (obs_aft_volume[4,2]+obs_harvest_volume[3,2]- obs_aft_volume[3,2])/5
lGz_obs_19 <- (obs_aft_volume[5,2]+obs_harvest_volume[4,2]- obs_aft_volume[4,2])/6
lGz_obs_25 <- (obs_aft_volume[6,2]+obs_harvest_volume[5,2]-obs_aft_volume[5,2])/6
lGz_obs_30 <- (obs_aft_volume[7,2]+obs_harvest_volume[6,2]-obs_aft_volume[6,2])/5
lGz_obs_35 <- (obs_aft_volume[8,2]+obs_harvest_volume[7,2]-obs_aft_volume[7,2])/5
lGz_obs_40 <- (obs_aft_volume[9,2]+obs_harvest_volume[8,2]-obs_aft_volume[8,2])/5
lGz_obs_45 <- (obs_aft_volume[10,2]+obs_harvest_volume[9,2]-obs_aft_volume[9,2])/5


lGz_obs_3_45 <- c(lGz_obs_3, lGz_obs_8, lGz_obs_13, lGz_obs_19, lGz_obs_25, lGz_obs_30, lGz_obs_35, lGz_obs_40, lGz_obs_45)
age_obs <- c(63, 66, 71, 76, 82, 88, 93, 98, 103) # stand age in years
lGz_obs <- data.frame(lGz= lGz_obs_3_45, age= age_obs, source = c(rep("observation (5-yr. average)", length(age_obs))),Sources = "observation")



# -----------------------------------------------------------------------


# Make a long dataframe out of model, best fit and observation lGz
lGz_models <- bind_rows(model_runs, .id = "model") # stack the models in the list into a single df
lGz_models$Sources <- "model parameter examples" # add colour column for ggplot
lGz_comparison <- rbind(lGz_models[,-1], lGz_obs, lGz_best_fit, lGz_orig) # finally, stack the models and the observed data into a single df

#lGz_comparison as a csv file as input to plot the lgz in precalibration_plot.R 
write.csv(lGz_comparison, "output/precalibration/lGz_comparison.csv") # load the model, observed and best fit runs


# -------------------- Models within acceptable range --------------------------

# acceptable is defined here as +- 2.5 m3/yr/ha divergence from observation
upper_range <- lGz_obs$lGz+2.5
lower_range <- lGz_obs$lGz-2.5
lower_range[5] <- lGz_obs$lGz[5]-6 # at age 83 the observed lGz differs heavily from all models, so we broaden the range here 
lower_range[7] <- lGz_obs$lGz[7]-3 # at age 88 larger diversion again

# test: is lGz of model i within the upper-lower range?
test_funct <- function(i){all(upper_range > model_runs[[i]]$lGz & model_runs[[i]]$lGz> lower_range) }

accepted <- lapply(c(1:n.samples), test_funct) # creates a list for the result of n.sample model runs
acc <- unlist(accepted) # change list into a vector
true <- which(acc == "TRUE") # a vector that shows which models are within acceptable range
params_acc <- params[true,] # subset parameter sets that are within acceptable range
param_min <- apply(params_acc, 2, min) # get the min value for each parameter (= for each column)
param_max <- apply(params_acc, 2, max) # and the max value


# Save global environment as R.Data for further use (not having to run the whole analysis again)
save.image(file = "precalibration/lGz_plot_apply_glob_environm.RData")

# Save acceptable parameter sets as csv as input for the robustness Sobol sensitivity analysis
write.csv(params_acc[,-23], "precalibration/Output/params_acc_scen_disturbed.csv") 

# The following calculations are usefule to prepare data for plotting with ggplot:

# Run the model_runs function again with only acceptable parameters and save
params <- params_acc # replace the parameter df with only the acceptable parameters
n.samples <- nrow(params)

model_runs <- lapply(c(1:n.samples), function(i){
  p <- as.numeric(params[i,]) # pick a parameter set
  
  model = forest.model(p, clim_dat, stand_dat, harvest_dat) # run the model
  
  # Separate the model output
  stand_scen_gp_table = model$stand
  harvest_stand_scen = model$harvest
  colnames(stand_scen_gp_table) <- c("Stand_scen","Year","Age","N","BA_stand","Dg","Ddom","Hdom","Volume_stand")
  colnames(harvest_stand_scen) <- c("Stand_scen","Year","Age","N","BA_harvested","VOL_h")
  
  # Subtracting "stand volume after thinning" at time = t-1 from "stand volume before next thinning" at time = t
  lGz_model_5 <- (stand_scen_gp_table[2,9] - stand_scen_gp_table[1,9])/5 # lgZ at age 63 (in 1969)
  lGz_model_10 <- (stand_scen_gp_table[4,9] - stand_scen_gp_table[3,9])/5
  lGz_model_15 <- (stand_scen_gp_table[6,9] - stand_scen_gp_table[5,9])/5
  lGz_model_20 <- (stand_scen_gp_table[8,9] - stand_scen_gp_table[7,9])/5
  lGz_model_25 <- (stand_scen_gp_table[10,9] - stand_scen_gp_table[9,9])/5
  lGz_model_30 <- (stand_scen_gp_table[12,9] - stand_scen_gp_table[11,9])/5
  lGz_model_35 <- (stand_scen_gp_table[14,9] - stand_scen_gp_table[13,9])/5
  lGz_model_40 <- (stand_scen_gp_table[16,9] - stand_scen_gp_table[15,9])/5
  lGz_model_45 <- (stand_scen_gp_table[18,9] - stand_scen_gp_table[17,9])/5
  
  lGz_model_5_45 <- c(lGz_model_5, lGz_model_10, lGz_model_15, lGz_model_20, lGz_model_25, lGz_model_30, lGz_model_35, lGz_model_40, lGz_model_45)
  age_model <- c(63, 68, 73, 78, 83, 88, 93, 98, 103)
  lGz_model <- data.frame(lGz= lGz_model_5_45, age= age_model, source = c(rep(paste(i, "model", sep = ""), length(age_model))))
  
})

lGz_models <- bind_rows(model_runs, .id = "model") # stack the models in the list into a single df
lGz_models$Sources <- "acceptable model parameter examples" # add colour column for ggplot
lGz_comparison_acc <- rbind(lGz_models[,-1], lGz_obs, lGz_best_fit, lGz_orig) # finally, stack the models and the observed data into a single df

write.csv(lGz_comparison_acc, "output/precalibration/lGz_comparison_acc.csv")
