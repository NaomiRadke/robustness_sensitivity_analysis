############################################################################################################################
#
# This function outputs the NPV of timber yield and the net carbon storage and can be sourced for 
# robustness analysis.
#
# This code is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#
# Any questions? Naomi Radke (naomikradke@gmail.com)
# 
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
#         - GDDI          sum of degree-days (April-October)
#         - dr            discount rate (4 scenarios)
#         - pr            timber price
#         - ct            carbon tax            
#
# Output:
# - output                dataframe with a row for each scenario and columns for the 2 objectives
#         - NPV_yield     NPV of harvested trees (yield) over the simulated period (Euro/ha)
#         - Carb_seq      net carbon storage in trees standing and wood products over the simulated period (ton C/ha) 
#
#########################################################################################################



obj.fun <- function(scenarios){  
  
  # Some settings 
  n.mod.params <- 22 # number of model parameters
  n_scen <- nrow(scenarios) # number of scenarios
  objectives <- data.frame(NPV_yield = numeric(0), Carb_seq = numeric(0)) # create an empty df for output of each scenario
  ind.model <-  1  # column in scenarios that correspond to model parameters
  ind.clim <- 2        # column that corresponds to rcp scenarios
  ind.dr <- 3             # column that corresponds to discount rate scenarios
  ind.pr <- 4
  wp_init <- 116 # price index from 2016
  
  # Go parallel for running n.scen scenarios
  
  cores = detectCores()-1   # Use one core less than available to not overload the PC
  cl <- makeCluster(cores)
  print(paste('Starting cluster with ',cores,' cores', sep=''))
  registerDoParallel(cl)
  
  # Run the forest growth model for each scenario
 output <-  foreach(i = 1:n_scen, .combine = rbind,# save output in a vector
                    .export= ls(globalenv()),
                    .inorder = FALSE) %dopar% {
                      scen_p <- as.numeric(scenarios[i,ind.model]) # get the model parameter scenario number
                      p <- as.numeric(param_scen[scen_p,]) # pick the parameters value of that model parameter scenario
                      
                      # Map timber price scenario number to respective time series data
                      scen_wp <- as.numeric(scenarios[i, ind.pr])
                      wp <- as.numeric(pr_scen[, scen_wp])
                      
                      # Add the timber price as column to the harvest input table so Revenue can be calculated within model
                      harvest_dat$pr_exp <- wp
                      
                      # Map climate scenario number to respective time series data
                      scen_clim <- as.numeric(scenarios[i, ind.clim]) # get the climate scenario number
                      clim_dat$DI <- as.numeric(clim_scen[, scen_clim]) # pick the time series values for the climate scenario
                      clim_dat$GDD <- 1
                      
                      # Map discount rate scenario number to respective time series data
                      scen_dr <- as.numeric(scenarios[i, ind.dr])
                      dr <- as.numeric(dr_scen[,scen_dr])
                      
                      
                      model.out <- forest.model(p, clim_dat, stand_dat, harvest_dat) # output of the model run for scenario i
                      
                      # Separate the harvest output for further use
                      harvest_out <- model.out$harvest  # subset the harvest output
                      colnames(harvest_out) <- c("Stand_scen","Year","Age","N","BA_harvested","VOL_h","Rev_sum", "TC", "u")
                      
                      stand_out <- model.out$stand # subset the stand output
                      colnames(stand_out) <- c("Stand_scen","Year","Age","N","BA_stand","Dg","Ddom","Hdom","Volume_stand")
                      
                      stand_final <- rev(model.out$individual)[[2]] # subset the final stand state (before growth) for further use
                      
                      # NPV TIMBER YIELD
 
                      Rev_end <- ((0.8+(rev(wp)[1]-14.88)*0.0368)*(stand_final$dbh-20))*stand_final$n # revenue at final stand state
       
                      Rev_end[Rev_end<0] <- 0 # replace negative revenues by 0 
                      NPV_end <- sum(Rev_end)/(1+rev(dr)[1])^rev(harvest_out$Year)[1] # Net Present Value
                      NPV_end_fc <- (sum(Rev_end) - 175)/(1+rev(dr)[1])^rev(harvest_out$Year)[1] # deduct 175 Euro/ha/yr fix cost for the fix cost version
                      Rev_init <- ((0.8+(wp_init-14.88)*0.0368)*(stand_dat$dbh-20))*stand_dat$ni
                      Rev_init[Rev_init<0] <- 0
                      NPV_init <- sum(Rev_init)  
                      NPV_init_fc <- NPV_init - 175
                      
                      NPV_yield_fc <- NPV_end_fc - NPV_init_fc + sum((harvest_out$Rev_sum-175)/((1+dr)^harvest_out$Year)) #
                      objectives[1,1] <- NPV_yield_fc #add value to the objectives dataframe for this scenario
                      
                      #NPV CARBON BALANCE
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
                      
                      # Calculate carbon balance for each year (growth - periodic emissions from wood products) 
                      Carb_seq <- sum(C_growth - C_em_wp)
                      
                      objectives[1,2] <- Carb_seq # add value to the objectives dataframe
                      return(objectives)
                       
                    } # END foreach scenarios
 
 return(output)

} # END obj.fun