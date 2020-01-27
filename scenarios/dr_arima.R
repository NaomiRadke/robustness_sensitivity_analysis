#################################################################################################################
#
# Create discount rate scenarios by fitting an ARIMA model
# Discount rate data: real interest rates on households' deposits (time deposits) 
# by Deutsche Bundesbank (https://www.bundesbank.de/en/statistics/money-and-capital-markets/interest-rates-and-yields/real-interest-rates-on-households-deposits) 
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
#- dr                   Time series data of real interest rates on houselhold deposits in Germany by Deutsche Bundesbank (1967-2019)
#         - date          date (dd.mm.yyyy)
#         - dr            real interest rates on households' deposits (time deposits)
#         
# Output:
#- sim                  Yearly time series data of 10 discount rate simulations from the fitted ARIMA model (2020-2080)
#- dr_plot              Yearly discount rate
#         - observation (1967-2019)
#         - simulations (2020-2080)      
#- sim_df_5yr           Matrix of 5-yearly time series data (rows) for the 10 discount rate simulations (columns) (2020-2080)      
#
#################################################################################################################

# load packages
library(zoo) # to calculate a moving average (using rollapply())
library(forecast) # for fitting arima and forecasting
library(ggplot2)  # for plotting
library(tseries) # for formatting as time series object


dr <- read.csv("data/scenarios/real_ir_deposits_bbk.csv", header = TRUE, sep = ";")

# change data into time series format
dr$date <- as.Date(dr$date, "%d.%m.%Y") # change date column to date format
dr$Year <- format(dr$date, "%Y") # add a year column in date format
dr$dr <- as.numeric(as.character(dr$dr)) # change column from factor to numeric format
dr_yearly <- aggregate(dr ~ Year, dr, mean) # calculate the yearly mean
dr_yearly$Year <- as.numeric(dr_yearly$Year) # change Year column to numeric so I can better convert it to time series
dr_orig_ts <- ts(dr_yearly$dr, start = 1967) # create time series object from original data

#STEP 1: Examine data
# check data
plot(dr_orig_ts)



# STEP 2: check stationarity with the sgumented Dickey-Fuller Test
adf.test(dr_orig_ts, alternative = "stationary")
# if p-value > 0.05 we cannot reject the H0 hypothesis, meaning our data is non-stationary and that differencing is required

# STEP 3: check autocorellations and choose order of ARIMA
acf(dr_orig_ts)
pacf(dr_orig_ts)

dr_orig_ts_d1 <- diff(dr_orig_ts, differences = 1) # differencing with level 1 (because only at lag 1 above line in PACF)
plot(dr_orig_ts_d1)

#--> spikes at particular lags of the differenced series can help inform choice of p and q

acf(dr_orig_ts_d1) # acf of differenced data --> q could be 3
pacf(dr_orig_ts_d1)# --> p could be 2

# --> in the next step evaluate the model parameter choice

# STEP 4: Evaluate and iterate ARIMA
fit_arma <- Arima(dr_orig_ts, order = c(2,0,3)) # fit an ARMA model 

tsdisplay(residuals(fit_arma)) # check residuals


# STEP 5: Forecast
fcast_arma <- forecast(fit_arma, h=60) # forecast 60 years into the future
plot(fcast_arma)



# STEP 6: FINALLY simulate scenarios
# Some settings
n.yrs = 60 # number of forecasted years
n.sim = 10 # number of sampled future trajectories

# Simulating n.sim future trajectories (from adding bootstrapped residuals to forecast)
sim <- ts(matrix(0, nrow = n.yrs, ncol = n.sim), start = end(dr_orig_ts)[1]+1) # empty ts matrix to be filled with simulations

# simulate and set seed for reproducability
for(i in 1:n.sim){
  sim[,i] <- simulate(fit_arma, nsim = n.yrs, seed = 100+i)
}

# save ts object
save(sim, file = "scenarios/discount_rate/Output/dr_simulations_fitarma.RData")

# Plot the recoreded time series and plausible future trajectories  
dr_plot <- autoplot(dr_orig_ts)+
  #ylim(c(0,))+
  autolayer(sim, colour = FALSE, col="blue")+
  labs(y = "real discount rate (%)", x = "time", title = "b)")+
  theme_bw()+
  theme(legend.position = "none")

save(dr_plot, file ="output/scenarios/dr_arma_plot.RData")

# change the ts object into a dataframe
sim_df <- as.data.frame(sim)

# calculate 5-yr averages for PI (moving averages) and keep only every 5th column
sim_df_5yr <- rollapply(sim_df, 5, FUN = mean, by =5, align = "left")
sim_df_5yr <- sim_df_5yr/100 # change from percentage to real number

# save the scenarios as csv file 
write.csv(sim_df_5yr, file = "output/scenarios/dr_scen_arma.csv")

