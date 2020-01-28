############################################################################################################################
#
# objective_function_yield_vect.R     28 January 2020
#
# Copyright (C) 2020 Naomi Radke
#
#
#
# This function calculates the Net Present Value of timber yield based on the forest growth 
# model output data (stand volume growth (m3/ha) and harvested volume (m3/ha) in 5-year intervals).
# It can be sourced for the Sobol Sensitivity analysis.
#
#
#
# objective_function_yield_vect.R is free software: you can redistribute it and/or modify
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
#- stand_dat            Beech stand data
#         - stand_scen    stand index number
#         - tree          tree number
#         - species       tree species (integer; i.e, 41)
#         - dbh           diameter at breast height (cm)
#         - age           age of tree
#         - ni            number of trees per hectare
#
#
#- harvest_dat          Harvesting data
#         - stand_scen    stand index number
#         - year          projection year when harvesting (at year 30, 60, 90, 100, and 110)
#         - age           age of trees harvested during harvest year
#         - htype         type of harvest thinning (low, uniform, or high thinning)
#         - g             thinning value (if present)
#         - g_r           desired basal area
#
#- scenarios              matrix with nrow scenarios
#         - p             model parameters
#         - DI            stand drought index
#         - dr            discount rate
#         - pr            timber price
#
# Output:
# - objectives            dataframe with one column and one row for each scenario
#         - NPV_yield     NPV of harvested trees (yield) over the simulated period (Euro/ha)
#  
#
#
####################################################################################################################  




obj.fun <- function(scenario){  
  
  # Some settings 
  n.mod.params <- 22 # number of model parameters
  
  objectives <- data.frame(NPV_yield = numeric(0)) # create an empty df for output of each scenario
  ind.model <-  1  # column in scenarios that correspond to model parameters
  ind.clim <- 2        # column that corresponds to rcp scenarios
  ind.dr <- 3             # column that corresponds to discount rate scenarios
  ind.pr <- 4
  wp_init <- 116 # price index from 2016

  # Run the forest growth model for each scenario
 
    
    # Map model scenario number to respective model parameter values
    scen_p <- as.numeric(scenario[ind.model]) # get the model parameter scenario number
    p <- as.numeric(param_scen[scen_p,]) # pick the parameters value of that model parameter scenario
    
    # Map timber price scenario number to respective time series data
    scen_wp <- as.numeric(scenario[ind.pr])
    wp <- as.numeric(pr_scen[, scen_wp])
    
    # Add the timber price as column to the harvest input table so Revenue can be calculated within model
    harvest_dat$pr_exp <- wp
    
    # Map climate scenario number to respective time series data
    scen_clim <- as.numeric(scenario[ind.clim]) # get the climate scenario number
    clim_dat$DI <- as.numeric(clim_scen[, scen_clim]) # pick the time series values for the climate scenario
    clim_dat$GDD <- 1
    
    # Map discount rate scenario number to respective time series data
    scen_dr <- as.numeric(scenario[ind.dr])
    dr <- as.numeric(dr_scen[,scen_dr])
   
  
    
    
    model.out <- forest.model(p, clim_dat, stand_dat, harvest_dat) # output of the model run for scenario i
    
    # Separate the harvest output for further use
    harvest_out <- model.out$harvest  # subset the harvest output
    colnames(harvest_out) <- c("Stand_scen","Year","Age","N","BA_harvested","VOL_h","Rev_sum", "TC", "u")
    
    stand_out <- model.out$stand # subset the stand output
    colnames(stand_out) <- c("Stand_scen","Year","Age","N","BA_stand","Dg","Ddom","Hdom","Volume_stand")
    
    stand_final <- rev(model.out$individual)[[2]] # subset the final stand state (before growth) for further use
    

    #NPV YIELD
    Rev_end <- ((0.8+(rev(wp)[1]-14.88)*0.0368)*(stand_final$dbh-20))*stand_final$n # Revenue at the final stand state
    
    Rev_end[Rev_end<0] <- 0 # replace negative revenues by 0 
    NPV_end <- sum(Rev_end)/(1+rev(dr)[1])^rev(harvest_out$Year)[1]
    NPV_end_fc <- (sum(Rev_end) - 175)/(1+rev(dr)[1])^rev(harvest_out$Year)[1] # deduct 175 Euro/ha/yr fix cost for the fix cost version
    Rev_init <- ((0.8+(wp_init-14.88)*0.0368)*(stand_dat$dbh-20))*stand_dat$ni
    Rev_init[Rev_init<0] <- 0
    NPV_init <- sum(Rev_init)  
    NPV_init_fc <- NPV_init - 175
    
    
    NPV_yield_fc <- NPV_end_fc - NPV_init_fc + sum((harvest_out$Rev_sum-175)/((1+dr)^harvest_out$Year)) #
    objectives <- NPV_yield_fc
    
  
  return(objectives)
  
} # END obj.fun