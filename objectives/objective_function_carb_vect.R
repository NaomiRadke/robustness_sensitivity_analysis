############################################################################################################################
#
# This function outputs the carbon storage in above-ground biomass and wood products and can be sourced for the
# Sobol Sensitivity analysis.
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
#- scenarios            Matrix with nrow scenarios
#         - p             model parameters
#         - DI            stand drought index
#         - dr            discount rate 
#         - pr            timber price index
#
# Output:
# - objectives          Dataframe with one column and one row for each scenario
#         - Net_C_seq     Net carbon sequestration in above-ground tree biomass and wood products (tons C/ha) over
#                         considered time period      
#
#
####################################################################################################################  




obj.fun <- function(scenario){  
  
  # Some settings 
  n.mod.params <- 22 # number of model parameters
  
  objectives <- data.frame(Net_C_seq = numeric(0)) # create an empty df for output of each scenario
  ind.model <-  1  # column in scenarios that correspond to model parameters
  ind.clim <- 2        # column that corresponds to rcp scenarios


  # Run the forest growth model for each scenario
 
    
    # Map model scenario number to respective model parameter values
    scen_p <- as.numeric(scenario[ind.model]) # get the model parameter scenario number
    p <- as.numeric(param_scen[scen_p,]) # pick the parameters value of that model parameter scenario
    
    
    # Map climate scenario number to respective time series data
    scen_clim <- as.numeric(scenario[ind.clim]) # get the climate scenario number
    clim_dat$DI <- as.numeric(clim_scen[, scen_clim]) # pick the time series values for the climate scenario
    clim_dat$GDD <- 1
    
    
    model.out <- forest.model(p, clim_dat, stand_dat, harvest_dat) # output of the model run for scenario i
    
    # Separate the harvest output for further use
    harvest_out <- model.out$harvest  # subset the harvest output
    colnames(harvest_out) <- c("Stand_scen","Year","Age","N","BA_harvested","VOL_h","Rev_sum", "TC", "u")
    
    stand_out <- model.out$stand # subset the stand output
    colnames(stand_out) <- c("Stand_scen","Year","Age","N","BA_stand","Dg","Ddom","Hdom","Volume_stand")
    
    stand_final <- rev(model.out$individual)[[2]] # subset the final stand state (before growth) for further use
    
   
    #NET CARBON BALANCE
    # Carbon balance = C-growth of stand - deadwood decomposition mass - periodic emissions from wood product pool
    
    Vol_growth <- 0               # initialize Volume growth vector
    
    for(t in 1:nrow(stand_out)){
      Vol_growth[t] <- stand_out[t+1,"Volume_stand"] - stand_out[t, "Volume_stand"]
    }
    Vol_growth <- Vol_growth[seq(1, length(Vol_growth)-1, 2)] # keep only every 2nd value as that is the growth
    
    # turn volume growth into carbon growth
    stem_dens <- 0.5583 # stem density of beech stem (Pistorius et al. 2006)
    vef <- 1.2 # volume expansion factor to total aboveground biomass (Pistorius et al. 2006)
    cf <- 0.5 # carbon factor (carbon in dry biomass) (IPCC GPG)
    C_growth <- Vol_growth*stem_dens*cf*vef
    
    # Periodic emissions of wood products
    
    # Create data frame with nrow = model period, ncol = periodic emissions of each harvest event
    period_em <- data.frame(matrix(nrow = nrow(harvest_out), ncol = nrow(harvest_out)))
    
    # loop through each harvest year and each model period to fill the matrix period_em
    for(year in 1:nrow(period_em)){
      TCy <- harvest_out$TC[year] # total carbon content of harvested trees
      uy <- harvest_out$u[year] # periodic emissions of trees harvested in that year
      for(period in 1:nrow(period_em)){
        if(TCy > uy){
          period_em[(period+(year-1)),year] <- uy
          TCy <- TCy-uy
        } else {
          period_em[(period+(year-1)),year] <- TCy
          TCy = 0
          break
        }  
      }
    }
    
    
    C_em_wp <- rowSums(period_em[c(1:nrow(harvest_out)),], na.rm = TRUE)  # calculate the total emissions from wood products for each period
    
    # Calculate carbon balance for each year (growth - periodic emissions from wood products) and multiply by the discounted C tax
    Carb_seq <- sum(C_growth - C_em_wp)
    Net_c_seq <- Carb_seq#+C_end - C_init
    
    objectives <- Net_c_seq
  
  # Calculate robustness (= % of scenarios in which minimum requirements regarding objectives are met)
  
  return(objectives)
  
} # END obj.fun