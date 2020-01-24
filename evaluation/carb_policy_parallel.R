#########################################################################################################################
#
# This script calculates how different levels of a carbon tax influence the economic robustness of the BAU strategy. 
# The economic robustness is the sum of the NPV of timber yield and the NPV of carbon sequestration.
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
#- obj_val              Dataframe of the objective values (NPV of timber yield and net carbon sequestration) for each scenario
#         - model         number of model scenario [1:84]
#         - DI            number of drought index scenario [1:20]
#         - DR            number of discount rate scenario [1:10]
#         - WP            number of timber price index scenario [1:10]
#         - NPV_yield     Net Present Value of timber yield from all harvest events in the considered time period
#         - Carb_seq      Net carbon sequestration in above-ground tree biomass and wood products in the considered time period
#         - robust        TRUE/FALSE minimum performance requirements for the two objectives are met
#
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
#- output               Vector of discounted carbon sequestration for every scenario
#         - disc_carb_seq discounted carbon sequestration over considered time period (ton C/ha)
#         
########################################################################################################################

# Load required packages
library(ggplot2) # for plotting
library(RColorBrewer) # to change plot colours
library(plot3D) # for 3D-plotting robustness
library(foreach) # for parallelizing
library(doParallel) # for
library(scales) # for converting to percentages

# Load required data

# the objective values (NPV timber yield and net carbon sequestration) for each scenario of the scenario ensemble
obj_val <- read.csv("evaluation/Output/robustness_bestRMSE272019-12-06.csv", row.names = 1) 

# model initialization data
stand_dat <- read.csv("Data/initial_state_trees_tending_wu_224.csv", header=TRUE, na.strings='(null)',sep=";", dec=".")
clim_dat <- read.csv("Data/climate_future_224.csv", header=TRUE, na.strings='(null)',sep=";", dec=".")
harvest_dat <- read.csv("Data/harvest_scenarios_Altherr.csv", header=TRUE, na.strings='(null)',sep=";", dec=".") # Altherr thinning

# Load scenarios
param_scen <- read.csv("Data/scenarios/params_acc_scen_disturbed.csv", row.names = 1)[,1:22] #
clim_scen <- read.csv("Data/scenarios/DI_scenarios.csv", row.names = 1) # nrow = nyears, ncol = 20 scenarios
pr_scen <- read.csv("Data/scenarios/PI_scen.csv", row.names = 1) # nrow = nyears, ncol = 10 scenarios price index
dr_scen <- read.csv("Data/scenarios/dr_scen_arma.csv", row.names = 1) # dr simulations with fitted ARMA

# Load models
source("growth_yield_stand_model_final.R") # source forest growth model
source("dominate_class_function.R") # source dominant class function needed for forest model


# Calculate the discounted net carb seq for each year and every scenario, where each scenario is a list


# create a list of discounted net carb seq for each scenario
n.samples <- nrow(obj_val)
carb_scen_disc <- data.frame(disc_carb_seq = numeric(0)) # create an empty df for output of each scenario

# Go parallel for running n.samples scenarios
cores = detectCores()-1   # Use one core less than available to not overload the PC
cl <- makeCluster(cores)
print(paste('Starting cluster with ',cores,' cores', sep=''))
registerDoParallel(cl)

# Run the forest growth model for each scenario
output <-  foreach(i = 1:n.samples, .combine = rbind,# save output in a vector
                   .export= ls(globalenv()),
                   .inorder = FALSE) %dopar% {

  
  # Some settings 
  ind.model <-  1  # column in scenarios that correspond to model parameters
  ind.clim <- 2        # column that corresponds to rcp scenarios
  ind.dr <- 3
  
  # Run the forest growth model for each scenario
  
  # Map model scenario number to respective model parameter values
  scen_p <- as.numeric(obj_val[i,ind.model]) # get the model parameter scenario number
  p <- as.numeric(param_scen[scen_p,]) # pick the parameters value of that model parameter scenario
  
  # Map climate scenario number to respective time series data
  scen_clim <- as.numeric(obj_val[i,ind.clim]) # get the climate scenario number
  clim_dat$DI <- as.numeric(clim_scen[, scen_clim]) # pick the time series values for the climate scenario
  clim_dat$GDD <- 1
  
  # Map discount rate scenario number to respective time series data
   scen_dr <- as.numeric(obj_val[i,ind.dr])
   dr <- as.numeric(dr_scen[,scen_dr])
  
  
  model.out <- forest.model(p, clim_dat, stand_dat, harvest_dat) # output of the model run for scenario i
  
  # Separate the harvest output for further use
  harvest_out <- model.out$harvest  # subset the harvest output
  colnames(harvest_out) <- c("Stand_scen","Year","Age","N","BA_harvested","VOL_h","Rev_sum", "TC", "u")
  
  stand_out <- model.out$stand # subset the stand output
  colnames(stand_out) <- c("Stand_scen","Year","Age","N","BA_stand","Dg","Ddom","Hdom","Volume_stand")
  
  stand_final <- rev(model.out$individual)[[2]] # subset the final stand state (before growth) for further use
  
  
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
  
  # Calculate carbon balance for each year (growth - periodic emissions from wood products) and multiply by the discounted C tax
  
  Net_carb_seq <- sum((C_growth - C_em_wp)/((1+dr)^harvest_out$Year)) #net carbon sequestration for the scenario
  
  carb_scen_disc[1,1] <- Net_carb_seq
  return(carb_scen_disc)
  
} # end lapply.. output: list with one value for each scenario

save(output,file = "evaluation/Output/carb_scen_disc")




# calculate NPV of carbon seq as a function of carbon tax

Carb_tax <- c(seq(0,120,0.5)) # carbon tax scenarios
ctax_scen <- data.frame(matrix(nrow = nrow(obj_val), ncol = length(Carb_tax)))

for (i in 1:length(Carb_tax)){
  
ctax_scen[,i] <- output$disc_carb_seq*Carb_tax[i] 
  
}

# calculate the net revenue (NPV yield + NPV carb) for each C-tax scenario in a new data frame
net_scen <- ctax_scen + obj_val$NPV_yield

# fraction of the number of negative net revenues for each C-tax scenario (column) 
robust <- colSums(net_scen>0)/nrow(net_scen)
robust <- percent(robust)


# Plot robustness vs carbon tax
robust_Ctax <- cbind.data.frame(Carb_tax, robust) # make a dataframe of the Ctax scenarios and respective robustness


ggplot(data = robust_Ctax, aes(x=Carb_tax, y = robust))+
  geom_point(aes(col = robust), size = 4)+
  labs(y= "Robustness (% of scenarios with an NPV > 0)", x = "Carbon Policy (Euro/ton C)")+
  scale_color_gradient(low = "red", high = "green", labels = scales:: percent_format(accuracy = 1))+
  theme_bw()+ 
  theme(text = element_text(size=14), legend.position = c(0.8,0.2))+
  labs(color = "Robustness")+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))
  
