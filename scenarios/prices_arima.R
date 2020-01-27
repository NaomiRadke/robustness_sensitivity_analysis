######################################################################################################
# This script fits an ARIMA model to a timber price index time series in order to simulate future paths.
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
#- WP                   Time series data of timber price index (averaged over products and species) in Germany (1977-2019)
#         - year          year (yyyy)
#         - price_index_all_wood - price index (%, where price in 2010 = 100%)
#         
# Output:
#- sim                  Yearly time series data of 10 price index simulations from the fitted ARIMA model (2020-2080)
#- PI_plot              Yearly price index
#         - observation (1977-2019)
#         - simulations (2020-2080)      
#- sim_df_5yr           Matrix of 5-yearly time series data (rows) for the 10 timber price simulations (columns) (2020-2080)      
#
#
######################################################################################################

# Set working directory
setwd("Analysis")

# Load required packages
library(forecast) # for fitting arima and forecasting
library(tseries)  # for turning data into time series data
library(ggplot2)  # for plotting
library(dplyr)
library(zoo)

# Load price data as time series
WP <- read.csv("data/scenarios/price_index_destatis_1977-18.csv", sep = ";")# read the data from csv
WP_ts <- ts(WP$price_index, start = 1977)                  # turn it into time series
plot.ts(WP_ts)                                             # check for stationarity 

#### MANUAL AND AUTOMATED ARIMA ######################################

# STEP 1: Test for stationarity/ Identify d

# If non-stationary, use diff to de-trend the data
#WP_ts_diff1 <- diff(WP_ts, differences = 1)
#plot(WP_ts_diff1)


# Conduct formal test for stationarity

  # e.g. with the Augmented Dickey-Fuller unit root test (p-val should be below 0.05/5%)
#  print(adf.test(WP_ts_diff1))


# STEP 2: Identifying p and q
  
  # making a correlogram and partial correlogram
#  acf(WP_ts_diff1)    # shows that q is 1
#  pacf(WP_ts_diff1)   # shows that p is 0
  
#  WP_ts_arima <- Arima(WP_ts, order = c(1,1,0), include.mean = FALSE) # This could be our arima model
  
  
  # Let's see what the automated arima function gives us 
  auto <- auto.arima(WP_ts, stepwise = FALSE, approximation = FALSE, seasonal = FALSE)  # let the auto.arima function do that for us
                                                                      # stepwise+approximation set to FALSE to 
                                                                      # consider all possible models in the search
  summary(auto) # shows the fitted arima model
  
# STEP 3: Check the residuals
  checkresiduals(auto)  # check for normality / ACF = whether autorcorrelations are within threshold
  
  
# STEP 4: Forecast  

  # for the manual arima
#  WP_forecast <- forecast(WP_ts_arima, h = 60)
#  summary(WP_forecast)
  
  # for the automated arima
  for_auto_arima <- forecast(auto, h = 60)
 
  # plot the best fit arima -> will add uncertain forecast scenario lines to the graph lateron
   g <- autoplot(for_auto_arima)
  


####### CREATING  SCENARIOS #########################
  
  # Some settings
  n.yrs = 60 # number of forecasted years
  n.sim = 10 # number of sampled future trajectories
  
  # Simulating n.sim future trajectories 
  sim <- ts(matrix(0, nrow = n.yrs, ncol = n.sim), start = end(WP_ts)[1]+1) # empty ts matrix to be filled with simulations
  
  # simulate and set seed for every run
  for(i in 1:n.sim){
    sim[,i] <- simulate(auto, nsim = n.yrs, seed = 100+i)
  }
  
  # save ts object
  save(sim, file = "output/scenarios/timb_pric_simulations.RData")
  
  # Plot the recoreded time series and plausible future trajectories  
  PI_plot <- autoplot(WP_ts)+
    xlim(c(1970,2080))+
    ylim(c(0,200))+
    autolayer(sim, colour = FALSE, col="blue")+
    labs(y = "price index (price in 2010=100%)", x = "time", title = "c)")+
    theme_bw()+
    theme(legend.position = "none")
    
  save(PI_plot, file ="PI_plot.RData")
  

  # change the ts object into a dataframe
  sim_df <- as.data.frame(sim)
  
  # calculate 5-yr averages for PI (moving averages) and keep only every 5th column
  sim_df_5yr <- rollapply(sim_df, 5, FUN = mean, by =5, align = "left")
  
  # change order of rows according to more or less declining price scenarios
  sim_df_5yr <- sim_df_5yr[,c("Series 7", "Series 6", "Series 9", "Series 10", "Series 8","Series 3", "Series 4","Series 2", "Series 5","Series 1")]
  
  
  # save the index scenarios as csv file 
  write.csv(sim_df_5yr, file = "output/scenarios/PI_scen.csv")
 