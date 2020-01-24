######################################################################################################################
#
# This script calculates the yearly and 5-yearly drought index DI and degree day index GDDI for the years 2018-2088 
# for RCPs 2.6, 4.5, 6.0 and 8.5 for each of the following CMIP5 models:
# 1.) GFDL-ESM2M 
# 2.) HadGEM2-ES 
# 3.) IPSL-CM5A-LR 
# 4.) MIROC-ESM-CHEM 
# 5.) NorESM1-M 
#
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
#- all_csv              List of climate hindcasts and forecasts of the above mentioned Global Climate Models and RCPs at stand location
#         - date          date (yyyy-mm-dd)
#         - prec          daily precipitation sum (mm)
#         - temp          daily temperature average (degree Celsius)
#
#- ts_obj               Data frame with observed climatic time series data at the forest stand cell (1950-2016)
#         - date          date(yyyy-mm-dd)
#         - temp          daily temperature average (degree Celsius)
#         - prec          daily precipitation sum (mm)
#         
#
# Output:               
#- DI                   Average Drought Index for every fith year (2018-2088) and every RCP and Global Climate Model
#         - DI            drought index ]0,1]
#
#- bg_DI                Average Drought Index over all RCP and Global Climate Models for every fifth year (2018-2088)
#         - DI            drought index ]0,1]
#
#- GDDI                 Average Degree Day Index for every fith year (2018-2088) and every RCP and Global Climate Model
#         - GDDI          degree day index ]0,1]
#
# ... and several plot data
# - DI_plot             Time series plot of yearly average drought index (April-October)
#          - observation (1950-2016)
#          - hindcasts (1950 -2018)
#          - forecasts (2018-2080)
#
#- temp_yearly          Average yearly mean temperature (degree Celsius, April-October)
#          - observation (1950-2016)
#          - hindcasts (1950 -2018)
#          - forecasts (2018-2080)
#
#- precip_yearly        Yearly precipitation sum (mm, April-October)
#          - observation (1950-2016)
#          - hindcasts (1950 -2018)
#          - forecasts (2018-2080)
#
######################################################################################################################

#------ General settings---#
setwd("C:/....") # set working directory
model_names <- c("gfdl_2.6", "gfdl_4.5", "gfdl_6.0", "gfdl_8.5", "hadgem_2.6", "hadgem_4.5", "hadgem_6.0", "hadgem_8.5", "ipsl_2.6",
            "ipsl_4.5", "ipsl_6.0", "ipsl_8.5", "miroc_2.6", "miroc_4.5", "miroc_6.0", "miroc_8.5", "noresm_2.6", "noresm_4.5",
            "noresm_6.0", "noresm_8.5") # vector with the model+rcp names to loop through

file_names <- as.vector(c(list.files("data/climate/ISIMIP_Data_Jared_stand_cell")[1:20]))

#-------Load packages -----#
library(lubridate) # to round daily data to monthly using function floor_date
library(dplyr) # to group data by month/ year etc. 
library(SPEI) # to calculate potential evapotranspiration
library(zoo) # to calculate a moving average (using rollapply())
library(ggplot2) # for plotting

#------ Load data ---------#
# Read all climate csv files and save in the list all_csv

all_csv <- lapply(file_names, function(i){read.csv(paste("data/climate/ISIMIP_Data_Jared_stand_cell/",i, sep = ""), 
                                                   header = TRUE, 
                                                   sep=",",
                                                   col.names = c("date", "prec", "temp", "Na", "Na"), 
                                                   stringsAsFactors=FALSE)[,1:3]}) # keep only columns 1 to 3
#-----
# Read the observed climate data
ts_obs <- read.csv("data/climate/climate_daily_record.csv", row.names = 1)[,c(1,3,2)] #read csv and change the column order

#-----------

# Set name of each list element
names(all_csv) <- model_names

#Add the observed data to the list to include it in the further process
all_csv[["observation"]] <- ts_obs

#Change the date columns from integer to date format
all_csv <- lapply(all_csv, transform, date = as.Date(date))



#------ Calculate the monthly, yearly (April-Oct) and 6-yearly average for data recorded and models------

#--MONTHLY----


  # create a monthly list
list_monthly <- list()

  # create column that is 1 if temperature is >= 5°C and otherwise 0
all_csv <- lapply(all_csv, function(x){x$is_degree_day <- ifelse(x$temp >5,1,0);return(x)})

 # add a column that gives the degree days (dd), multiplying is_degree_day with temp
all_csv <- lapply(all_csv, function(x){x$dd <- x$is_degree_day * x$temp;return(x)})

# summarize temp (mean), prec (sum) and dd (sum) monthly
for (i in 1:21){ list_monthly[[i]] <- all_csv[[i]] %>%  
  group_by(month = floor_date(date, "month")) %>% 
  summarise(temp = mean(temp), 
            prec = sum(prec), 
            sdd = sum(dd))}

# calculate monthly potential evapotranspiration using the SPEI package
list_monthly <- lapply(list_monthly, function(x) {x$PET <- as.numeric(thornthwaite(x$temp, 49.25, na.rm=FALSE));return(x)}) # Change time series (output of thornthwaite) into vector

# Add the fitting model name to each dataframe of the list
names(list_monthly) <- append(model_names, "observation")

#--SEASONALLY----

# create a season list (sum of April - October for each year)
list_season <- list()

# sum up prec, sdd and PET for each season
for (i in 1:21){ list_season[[i]] <- list_monthly[[i]] %>%  
group_by(year = floor_date(month, "year")) %>%
  filter(month(month) > 3 & month(month) < 11) %>% summarise( prec = sum(prec),
                                                              sdd = sum(sdd),
                                                              PET = sum(PET),
                                                              temp = mean(temp))}

# Add the fitting model name to each dataframe of the list
names(list_season) <- append(model_names, "observation")

# calculate DI and GDDI
#DI
a <- 0.78 
WHC_m <- 150 #soil water holding capacity (for 1m of soil depth), originally 210 (Trasobares et al. 2016)
DI = function(x){x$DI <- (a*ifelse(x$prec/x$PET>1,1, x$prec/x$PET))+((1-a)*(WHC_m/210));return(x)}
list_season <- lapply(list_season, DI)

#GDDI
thr = 3000 # original:1900 -> threshold value.. when sum of degree days is lower it has an effect on dbh growth
GDDI = function(x){x$GDDI <- ifelse(x$sdd/thr >1,1,x$sdd/thr) ; return(x)}
list_season <- lapply(list_season, GDDI)

#--6-YEAR-AVERAGE----
list_averages <- list()

#lapply(list_season, function(x) {rollapply(x$prec, 6, mean)})

for(i in 1:21){ list_averages[[i]] <- list_season[[i]] %>%
  mutate(prec_6_av = rollapply(data = prec, width = 6, FUN = mean, fill = NA, na.rm = T, align = "right"),
         GDDI_6_av = rollapply(data = GDDI, width = 6, FUN = mean, fill = NA, na.rm = T, align = "right"),
         DI_6_av = rollapply(data = DI, width = 6, FUN = mean, fill = NA, na.rm = T, align = "right"))
  
}
# Add the fitting model name to each dataframe of the list
names(list_averages) <- append(model_names, "observation")

# Save the "DI" values of every model in a dataframe
DI <- sapply(list_averages, function(x) x[["DI"]]) # extract the DI column from every dataframe in the list
DI <- as.data.frame(DI[-21]) # turn it into a dataframe but removing the observed data
DI$year <- list_averages$gfdl_2.6$year # add a column with years
DI <- DI %>% filter(year >= "2019-01-01" & year <= "2079-01-01")%>% 
  filter(row_number() %% 5 ==1)%>% # delete all years before 2024 and keep only every 5th row starting with the 1st row
  select(- year) # delet the year column since we don't need it anymore

model_order <- c("gfdl_2.6", "hadgem_2.6", "ipsl_2.6", "miroc_2.6", "noresm_2.6", 
                 "gfdl_4.5","hadgem_4.5", "ipsl_4.5", "miroc_4.5","noresm_4.5",
                 "gfdl_6.0", "hadgem_6.0", "ipsl_6.0", "miroc_6.0", "noresm_6.0",
                 "gfdl_8.5", "hadgem_8.5", "ipsl_8.5", "miroc_8.5", "noresm_8.5")

DI <- DI[, model_order]

write.csv(DI, "data/climate/Output/DI_scenarios.csv")

# "best guess" DI scenario (average of all scenarios)
bg_DI <- rowMeans(DI)

write.csv(bg_DI, "data/climate/Output/bg_DI_scen.csv")

# Save the "GDDI" values of every model in a dataframe
GDDI <- sapply(list_averages, function(x) x[["GDDI"]]) # extract the GDDI column from every dataframe in the list
GDDI <- as.data.frame(GDDI[-21]) # turn it into a dataframe but removing the observed data
GDDI$year <- list_averages$gfdl_2.6$year # add a column with years
GDDI <- GDDI %>% filter(year >= "2019-01-01" & year <= "2079-01-01")%>% 
  filter(row_number() %% 5 ==1)%>% # delete all years before 2024 and keep only every 5th row starting with the 1st row
  select(- year) # delet the year column since we don't need it anymore

model_order <- c("gfdl_2.6", "hadgem_2.6", "ipsl_2.6", "miroc_2.6", "noresm_2.6", 
                 "gfdl_4.5","hadgem_4.5", "ipsl_4.5", "miroc_4.5","noresm_4.5",
                 "gfdl_6.0", "hadgem_6.0", "ipsl_6.0", "miroc_6.0", "noresm_6.0",
                 "gfdl_8.5", "hadgem_8.5", "ipsl_8.5", "miroc_8.5", "noresm_8.5")

GDDI <- GDDI[, model_order]

write.csv(GDDI, "data/climate/Output/GDDI_scenarios.csv")

#------ Plot the monthly, yearly (April-Oct) and 6-yearly average for data recorded and models------


# Yearly (=seasonally) for DI and GDDI
df <- bind_rows(list_season, .id = "model") # stack the models in the list into a single df
df <- df %>% filter(year <= "2080-01-01")

# DI
DI_plot <- ggplot(df, aes(x=year,y=DI, group = model)) + 
  ylim(c(0.3,1))+
  scale_x_date(breaks=as.Date(c("1960-01-01", "2000-01-01", "2040-01-01", "2080-01-01")), date_labels = "%Y", limits = as.Date(c("1950-01-01", "2080-01-01")))+
  geom_line(colour = "blue")+
  geom_line(data = subset(df, model == "observation"), colour = "black")+
  labs(y="drought index (]0,1])", x= "time", title = "a)")+
  theme_bw()+
  theme(legend.position = "none")

save(DI_plot, file = "scenarios/climate_scenarios/DI_plot.RData")

ggplot(df[seq(1,nrow(df), 5),], aes(year,DI, colour = model)) + 
  ylim(0,1)+
  geom_line()

# GDDI
ggplot(df, aes(year,GDDI, colour = model)) + 
  ylim(0,1)+
  geom_line()

# plot every 5th year
ggplot(df[seq(1,nrow(df), 5),], aes(year,GDDI, colour = model)) + 
  geom_line()

# temperature - and save as tiff
tiff("plots/temp.tiff",  width = 7, height = 4.5, units = "in" ,  res = 600)
ggplot(df, aes(x=year, y=temp, group = model))+
  scale_x_date(breaks=as.Date(c("1960-01-01", "2000-01-01", "2040-01-01", "2080-01-01")), date_labels = "%Y", limits = as.Date(c("1950-01-01", "2080-01-01")))+
  geom_line(colour = "blue")+
  geom_line(data = subset(df, model == "observation"), colour = "black")+
  labs(y="average seasonal temperature (April-October, °Celsius)", x= "time")+
  theme_bw()+
  theme(legend.position = "none")
dev.off()

# precipitation
tiff("plots/prec.tiff",  width = 7, height = 4.5, units = "in" ,  res = 600)
 ggplot(df, aes(x=year, y=prec, group = model))+
  scale_x_date(breaks=as.Date(c("1960-01-01", "2000-01-01", "2040-01-01", "2080-01-01")), date_labels = "%Y", limits = as.Date(c("1950-01-01", "2080-01-01")))+
  geom_line(colour = "blue")+
  geom_line(data = subset(df, model == "observation"), colour = "black")+
  labs(y="average seasonal precipitation (April-October, mm)", x= "time")+
  theme_bw()+
  theme(legend.position = "none")
dev.off()

