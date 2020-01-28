################################################################################################################
#
# Sobol_driver_carb.R     28 January 2020
#
# Copyright (C) 2020 Naomi Radke
#
#
#
# This script performs a Sobol sensitivity analysis to rank the impact of the different uncertain parameters 
# on the net carbon sequestration.
# It sources the objective function which in turn sources the forest growth model.
#
#
#
# Sobol_driver_carb.R is free software: you can redistribute it and/or modify
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
#         - dbh           diameter at breast height (dbh)
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
# Output:
#- Sobol_out_carb       Output by the SobolSalt function which also provides the sensitivity indices. See ?sobolSalt for documentation
#
#- sens_carb            Dataframe with 1st and total order sensitivity indices for all uncertain factors
#         - sens_type     "1st" for 1st order sensitivity or "total" for total order sensitivity
#         - params        "model" for model parameters or "drought index" as drought index as uncertain factor
#         - sens_index    1st or total order sensitivity indices [0,1], the higher the stronger their contribution to total variance
#         - CI_low        lower confidence interval value for sensitivity indices
#         - CI_up         upper confidence interval value for sensitivity indices
#
#
####################################################################################################################
# set working directory
#setwd("~/...")

# load required packages
library(ggplot2) # for plotting
library(doParallel) # for going parallel
library(foreach) # for going parallel
library(sensitivity) # for Sobol sensitivity analysis

# load data to run forest model
stand_dat <- read.csv(".csv", header=TRUE, na.strings='(null)',sep=";", dec=".")
clim_dat <- read.csv(".csv", header=TRUE, na.strings='(null)',sep=";", dec=".")
harvest_dat <- read.csv(".csv", header=TRUE, na.strings='(null)',sep=";", dec=".") # Altherr thinning

# load input tables with scenarios
param_scen <- read.csv("data/scenarios/params_acc_scen_disturbed.csv", row.names = 1) # model parameter scenarios
clim_scen <- read.csv("data/scenarios/DI_scenarios.csv", row.names = 1) # climate scenarios, nrow = nyears, ncol = 20 scenarios

# Some settings
n.boot <- 700                                 # number of bootstraps for Sobol analysis
n.param <- 2                                  # number of parameters


# create n.scen scenarios using LHS
set.seed(2017) # set seed for reproducability
lhs_scen_carb <- data.frame(randomLHS(n.scen, n.param))


# Source the required functions
source("objective_function_carb_vect.R")  # source the objective function that calculates objective values for NPV timber and carbon
#(the following functions are part of the growth model and are not provided here)
source("growth_yield_stand_model_final.R") # source forest growth model
source("dominate_class_function.R") # source dominant class function needed for forest model



###--------- Calculate Sobol sensitivity indices (1st order and total) -------------------------------------------

# vector with objects to be exported/ required for running parallel
export.names <- c("obj.fun", "f.d_dom", "forest.model", "Carbon_fct", 
                  "clim_dat", "stand_dat", "harvest_dat", 
                  "param_scen", "clim_scen")




##---Sobol NPV carbon----
# Create function to feed into Sobol that outputs an NPV timber value
Carbon_fct <- function(dataframe.in){ # dataframe.in will be x1 and x2 of the Sobol
  
  # match the [0,1] scenarios to their actual scenario number
  scen <- dataframe.in # initialize a dataframe for the scenarios fed into the objective function
  
  scen[,1] <- cut(dataframe.in[,1], breaks = seq(0,1,1/nrow(param_scen)), labels = FALSE) # match model param scenarios
  scen[,2] <- cut(dataframe.in[,2], breaks = seq(0,1,1/ncol(clim_scen)), labels = FALSE) # match drought index scenarios
  
  
  # initialize the output
  nr <- nrow(scen)
  output <- rep(0,nr)
  
  
  # start parallelization
  cores = detectCores()-1   # Use one core less than available to not overload the PC
  cl <- makeCluster(cores)
  print(paste('Starting cluster with ',cores,' cores', sep=''))
  registerDoParallel(cl)
  
  # Run model in parallel
  print(paste('Starting ',nr,' model projections ...',sep=''))
  final.output <- foreach(i=1:nr, .combine = c,# save output in a vector
                          #.export= ls(globalenv()),
                          .export = export.names,
                          .inorder = FALSE) %dopar% {
                            model.out <- obj.fun(as.numeric(scen[i,]))  # run the model
                            output[i] <- as.numeric(model.out) # store the model output of each run in the output vector                           
                          }# END foreach
  
  
  
  print(paste(' ... done.'))
  
  stopCluster(cl)
  
  # in case the Sobol function requires a zero-mean output, output should be centered at 0
  output.avg <- mean(final.output)
  final.output <- final.output - output.avg
  
  return(final.output)   # return output to Sobol analysis
  
} # END function NPV_carbon_fct

# Define the two samples X1 and X2 for Sobol sensitivity analysis

size = 80000 # how many samples to draw for X1 and X2 each
set.seed(135)
X1 = lhs_scen_carb[sample(c(1:nrow(lhs_scen_carb)), size = size),]
set.seed(246)
X2 = lhs_scen_carb[sample(c(1:nrow(lhs_scen_carb)), size = size),] 



# Run Sobol for carbon seq
Sobol_out_carb <- sobolSalt(model = Carbon_fct,
                            X1,
                            X2,
                            scheme = "A",
                            nboot = n.boot)


# Save output
currentDate <- Sys.Date()
FileName <- paste("output/evaluation/Sobol_out_carb", currentDate, ".RData", sep = "")
save(Sobol_out_carb, file = FileName)



#plot(sobol_out, choice = 1)

# Test: Did we have enough samples? If largest width of CI is less than 10% of the highest total sensitivity index (Wong 
# & Keller 2017). This shows that the CI width converges (doesn't get smaller with an increased number of samples/boots)


params <- c("model", "drought index")
sens_1st_tot_carb <- data.frame(cbind(params,
                                 "sens_1st" = Sobol_out_carb$S[,1], #sensitivity of each parameter = Var(param_i)/Var(model output)
                                 "CI_low_1st" = Sobol_out_carb$S[,4],
                                 "CI_up_1st" = Sobol_out_carb$S[,5],
                                 "sens_tot" = Sobol_out_carb$T[,1], #total sensitivity of each parameter = 1- (Var_model - Var_param_i)/Var_model
                                 "CI_low_tot" = Sobol_out_carb$T[,4],
                                 "CI_up_tot" = Sobol_out_carb$T[,5]))

max_CI_width <- 0.1* max(as.numeric(paste(sens_1st_tot_carb$sens_tot))) # the largest CI width should be 10% of the largest total sensitivity index
larg_CI_width <- max(sens_1st_tot$CI_up_tot - sens_1st_tot$CI_low_tot) # the actual largest CI with of total sens.
larg_CI_width <= max_CI_width                       # test: is largest CI small enough? TRUE or FALSE.. if FALSE, increase 
# n.sample and n.boot      

# Make a nicer sensitivity plot than what soboSalt provides

# First, create a dataframe for use with ggplot
sens_carb <- data.frame(sens_type = rep(c("1st", "total"), each = 2),
                        
                        params = rep(params[1:2], 2),
                        sens_index = c(as.numeric(paste(sens_1st_tot_carb$sens_1st))[1:2], as.numeric(paste(sens_1st_tot_carb$sens_tot))[1:2]),
                        CI_low = c(as.numeric(paste(sens_1st_tot_carb$CI_low_1st))[1:2], as.numeric(paste(sens_1st_tot_carb$CI_low_tot))[1:2]),
                        CI_up = c(as.numeric(paste(sens_1st_tot_carb$CI_up_1st))[1:2], as.numeric(paste(sens_1st_tot_carb$CI_up_tot))[1:2]))

sens_carb$params <- factor(sens_carb$params, levels = unique(sens_carb$params)) # turn column into a factor to get a correct order when plotting

# create bar plot
sobol_carb <- ggplot(data=sens_carb, aes(x= params, y= sens_index, fill = sens_type))+
  geom_bar(stat="identity", position = position_dodge(), width = 0.5)+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width = 0.2, position = position_dodge(0.5))+
  theme_bw()+
  scale_fill_manual(values = c("grey", "grey50"))+
  labs(y = "Sobol sensitivity indices", x = "", fill = "sensitivity type", title = "Net carbon sequestration")+
  scale_y_continuous(breaks=seq(0,1, 0.2), limits = c(0,1))

