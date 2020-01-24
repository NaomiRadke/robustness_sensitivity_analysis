######################################################################################################################
# 
# This script uses a global optimizing algorithm (DEoptim) to fit the parameters of the dbh andheight growth function 
# and survival function by the beech growth model of Trasobares et al. 2016 to the beech stand tested in this study.
#
# This code is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
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
#- stand_dat            Beech stand data
#         - stand_scen    stand index number
#         - tree          tree number
#         - species       tree species (integer; i.e, 41)
#         - dbh           diameter at breast height (cm)
#         - age           age of tree (years)
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
# - params_init         Matrix of samples of the original parameters +-5*standard error
#         - X1            parameter 1
#         - X2            parameter 2
#         - ...            
#         - X22           parameter 22
#
# Output:
# - out_DEoptim         List of 2 containing the results of the global optimization algorithm
#         - optim       List containing ptimal model parameter results
#               - bestmem   the best parameter set found
#               - bestval   best function evaluation value reached (in our case Root Mean Squared Error of yearly volume growth)
#               - nfeval    number of function evaluations
#               - iter      number of procedure iterations
#         - member      List containing all our settings and results of the optimization search process
#               - lower     lower bound of model parameter values we specified
#               - upper     upper bound of model parameter values we specified
#               - bestmemit best parameter set at each iteration
#               - bestvalit best function value (RMSE of yearly volume growth) at each iteration
#               - pop       the population generated at the last iteration
#               - storepop  list containing the intermediate populations
#
# - parameters_sample_all Matrix of n.sample Latin Hypercube Samples over +- 5% * RMSE of the optimal model parameter set 
#         - X1            parameter 1
#         - X2            parameter 2
#         - ...            
#         - X22           parameter 22
#
#
#
####################################################################################################################




# Clear the environment 
rm(list=ls())
graphics.off()

# Load the DEoptim R package
library(DEoptim)

# Load data for the forest.model 
# initialization data
stand_dat <- read.csv("Data/initial_state_trees_tending_wu_224.csv", header=TRUE, na.strings='(null)',sep=";", dec=".")
clim_dat <- read.csv("Data/climate_record_224.csv", header=TRUE, na.strings='(null)',sep=";", dec=".")
harvest_dat <- read.csv("Data/harvest_scenarios_scaled_224_wu.csv", header=TRUE, na.strings='(null)',sep=";", dec=".")

# observations
obs_harvest_dat <- read.csv("Data/individ_trees_harvest_224.csv", header=TRUE, na.strings='(null)',sep=";", dec=".") # trees harvested for each harvest year
obs_stand_aft <- read.csv("Data/indiv_dbh_h_understory_aft_harv_224.csv", header=TRUE, na.strings='(null)',sep=";", dec=".") # stand data after harvest


# For the observed data, calculate volume of each tree in order to lateron calculate stand volume
# First, calculate volume of each tree in order to lateron calculate stand volume
obs_stand_aft$volume <- (pi/4)*(obs_stand_aft$dbh/100)^2*obs_stand_aft$height*0.5*obs_stand_aft$n
# Trees harvested
obs_harvest_dat$volume <- (pi/4)*(obs_harvest_dat$dbh/100)^2*obs_harvest_dat$height*0.5*obs_harvest_dat$n_harv

# Calculate stand volume after harvest of each observation year
obs_aft_volume <- aggregate(volume~year, data = obs_stand_aft, sum)
obs_harvest_volume <- aggregate(volume~year, data = obs_harvest_dat, sum)



# Source the function to calculate the dominate class:
source("dominate_class_function.R")

# Source the model
source("fitting_vol_growth.R")


# load samples of the original parameters +-5*standard error
params_init <-  read.csv("Data/LHS_params_scen_large.csv", row.names = 1) 


# define bounds within which the algorithm searches (min and max value of each parameter)
bound.lower <- apply(params_init, 2, min) # the min value of each column 
bound.upper <- apply(params_init, 2, max) # the max value of each column

# define a function that calculates the RMSE of model vs observed data and gives a high RMSE if model does not function
model_RMSE = function(p, clim_dat, stand_dat, harvest_dat){
  
  tryCatch({
    forest.model(p, clim_dat, stand_dat, harvest_dat)
  }, error= function(e){
    return(1000) # return RMSE = 1000 in case an error occurs
  }) #END tryCatch
} # END model_RMSE function

# Now let the DEoptim algorithm find parameter values that lead to lowest RMSE
out_DEoptim <- DEoptim(model_RMSE, bound.lower, bound.upper, clim_dat = clim_dat, stand_dat = stand_dat, harvest_dat = harvest_dat, 
                       DEoptim.control(itermax = 500, trace = FALSE))
currentDate <- Sys.Date()
FileName <- paste("OutputDEoptim_out", currentDate, ".RData", sep = "")
save(out_DEoptim, file = FileName)


# Find a parameter values range with a low RMSE (RMSE+0.1*RMSE) to sample scenarios with an acceptable fit
RMSE_optim <- out_DEoptim$optim$bestval    #best RMSE
RMSE_max <- 1.05*RMSE_optim                 #max. allowable RMSE
iter_max <- which(out_DEoptim$member$bestvalit < RMSE_max)[1] # first iteration where RMSE < RMSE_max

# for each parameter select the max and min value at iteration (row) iter_max and after
val_max <- apply(out_DEoptim$member$bestmemit[iter_max:500,],2,max) # find max value of each parameter
val_min <- apply(out_DEoptim$member$bestmemit[iter_max:500,],2,min) # find max value of each parameter

# Conduct Latin Hypercube samples to be used to find parameter scenarios with an acceptable fit (see script: lGz_plot_lapply.R)

library(lhs) # Latin Hypercube Sampling package

# set bounds over which to sample
bound.upper <- val_max
bound.lower <- val_min

# some settings
n.parameter<- 22
set.seed(1234)
n.sample <- 100000

# do the sampling
parameters.sample <- data.frame(randomLHS(n.sample, n.parameter))

# match sample from [0,1] range to actual parameter values
parameters.sample.all <- parameters.sample # initialize

for (p in 1:n.parameter) {
  parameters.sample.all[,p] <- bound.lower[p] + (bound.upper[p]-bound.lower[p])*(parameters.sample[,p])
}

write.csv(parameters.sample.all, "precalibration/Output/params_DEoptim_disturbed.csv")
