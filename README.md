# Robustness and sensitivity analysis

This repository contains the analysis scripts, input- and output data on which the results of the manuscript "Identifying decision-relevant uncertainties for dynamic adaptive forest management" are based.
The analysis scripts are written in the R Programming language.

It contains the following scripts/data in the following sub-folders:
### /precalibration
- **DEoptim_calibration_parallel.R**: Runs a global opimization algorithm (DEoptim) over a defined value range for each growth model parameter
- **lGz_plot_lapply.R**: Calculates the periodic stand volume growth (m3/ha/yr) for a sample of model parameter combinations to find an acceptable subset

### /objectives
- **objective_function_parallel.R**: calculates the 2 objectives Net Present Value (NPV) of timber yield and net carbon sequestration using forest growth model output (stand volume growth and volume harvested for every period). Sourced in robustness.R
- **objective_function_yield_vect.R**: calculates the NPV of timber yield using forest growth model output (stand volume growth and volume harvested for every period). Sourced in Sobol_driver_yield.R
- **objective_function_carb_vect.R**: calculates the net carbon sequestration using forest growth model output (stand volume growth and volume harvested for every period). Sourced in Sobol_driver_carb.R

### /scenarios 
- **climate_scenarios.R**: time-variable climate scenarios
- **prices_scenarios.R**: time-variable timber price index scenarios
- **dr_arima.R**: time-variable discount rate scenario

### /evaluation
- **robustness.R**: calculates the objectives' performance under each scenario and Y/N if minimum performance requirements are met
- **carb_policy_parallel.R**: calculates the effect of introducing different levels of carbon taxes on the strategy's robustness
- **Sobol_driver_yield.R**: calculates the 1st and total order sensitivities i.e. the relative effect of every uncertainty on NPV of timber yield
- **Sobol_driver_carb.R**: calculates the 1st and total order sensitivities i.e. the relative effect of every uncertainty on net carbon sequestration

### /data
--> contains the data files used for precalibration and scenario creation. NOTE that the model parameter scenarios are based on the precalibration results and the data for the model parameter scenario creation can be found in the /output/precalibration folder.

### /data/precalibration
- LHS_params_scen_large.csv			


### /data/scenarios
- **price_index_destatis_1976-18.csv**: yearly data (1976-2018) for the price index of all wood harvested in Germany
- **real_ir_deposits_bbk.csv**: monthly data on the real interest rate on household deposits by the Deutsche Bundesbank (1967-2019) 
- **readme.txt**: information on the price index data
- **ISIMIP_Data_stand_cell folder**: contains daily precipitation (mm) and temperature (°C) data (1950-2016) for the studied stand for 5 different global climate models/ 4 RCPs (.csv format + readme.txt)

### /output           
--> contains the output files of the growth model precalibration, scenario creation and evaluation (robustness, sensitivity, carbon policy) scripts

### /output/precalibration
- **DEoptim_out_5SE2019-11-29.RData**: RData containing, among others, the optimal model parameter set found by the DEoptim algorithm
- **lGz_comparison.csv**: matrix containing the periodic volume growth (m3/ha/yr) for observed, best fit and a set of model parameter scenarios
- **lGz_comparison_acc.csv**: matrix containing the periodic volume growth (m3/ha/yr) for observed, best fit and a set of acceptably performing model parameter scenarios
- **model_runs.RData**: list for each model parameter scenario containing periodic volume growth (m3/ha/yr) for each time period 
- **params_DEoptim_disturbed.csv**: sample of model parameters that diverge maximum 0.05* Root Mean Squared Error (in periodic stand volume growth) from the optimal model parameter set (found when running DEoptim script)

### /output/scenarios                            
- **DI_scenarios.csv**: climate (drought index) scenarios
- **dr_scen_arma.csv**: discount rate scenarios
- **dr_simulations_fitarma.RData**: time series of each discount rate simulation
- **params_acc_scen_disturbed.csv**: model parameter scenarios
- **PI_scen.csv**: timber price index scenarios 
- **timb_pric_simulations.RData**: time series of each timber price index simulation 

### /output/evaluation
- **carb_scen_disc.RData**: vector containing discounted net carbon sequestration for every scenario
- **robustness_2019_12_06.csv**: matrix containing the objective values (NPV timber yield and net carbon sequestration)  as well as Y/N are minimum performance requirements fulfilled for every scenario
- **Sobol_out_carb_2019-12-06.RData**: RData containing (among others) the 1st order and total sensitivity indices + confidence intervals for the effect on net carbon sequestration
- **Sobol_out_NPVyield_2019-12-07.RData**: RData containing (among others) the 1st order and total sensitivity indices + confidence intervals for the effect on NPV of timber yield


Note that in each of these scripts, some edits will be necessary:
- file names with date-stamps in the names
- folders that are pointed to for in- or output data

This repository DOES NOT contain:
- the forest growth model we used to actually run the analysis. It is a climate-sensitive, empirical beech growth model created by Trasobares et al. (2016) and was adapted for use in this study. 
- the data used to:
	- calibrate the model: individual tree long-term stand data (diameters and numbers at different stand ages)
	- input data to the model (initial stand data, climate projection data, harvest pattern data)
	
	
We share these scripts to make the individual analysis steps transparent. We hope readers find them useful and can transfer them to their individual contexts of decision-making under deep uncertainty.


## WORKFLOW

![Analysis flow and matching scripts.](analysis_flow.jpg)

### --> precalibration of the growth model parameters
Run DEoptim_calibration_parallel.R to find a best fit parameter set for the forest growth and yield model. We use the global optimization algorithm DEoptim() from the DEoptim R package to run the optimization between the range of the original model parameter values (Trasobares et al. 2016) +/- 5 * standard error. This script outputs a large sample of perturbations of this optimal parameter set.

Run lgz_plot_lapply.R which will calculate the stand's periodic volume growth for the sample created in DEoptim_calibration_parallel.R, for the observed data, for the best fit model parameter values and the original parameter values (by Trasobares et al. 2016). It subsets those samples that show an acceptable fit to the observed data. This subset represents the model parameter uncertainty and will be used to evaluate robustness and sensitivity of a thinning strategy to the ensemble of uncertainties.

### --> scenario creation
Run climate_scenarios.R, prices_scenarios.R and dr_arima.R to create the exogeneous uncertainty scenarios, namely climate, timber price and discount rate scenarios. They are all time-varying scenarios.

### --> objectives
The objective functions are wrappers around the beech growth model. Using the model output (stand's volume growth and harvested volume at different points in time) they calculate the two objectives: Net Present Value of timber yield and net carbon storage. The objective_function_parallel.R is sourced in the robustness.R and the objective_function_yield_vect.R and objective_function_carb_vect.R in Sobol_driver_yield.R and Sobol_driver_carb.R script, respectively.

### ---> evaluation of robustness and sensitivities
Run robustness.R which i) creates an uncertainty ensemble by enumerating all model parameter, climate, price and discount rate scenarios, ii) calculates the objective values for each uncertainty scenario using the objective_function_parallel.R and iii) derives the robustness of the current management strategy (Altherr thinning) which is the fraction of scenarios in which minimum performance requirements are fulfilled.

Run carb_policy_parallel.R which runs the forest growth model under each scenario to calculate the discounted net carbon storage in above-ground biomass and wood products. It then calculates how different carbon prices influence the level of robustness (i.e. the fraction of scenarios in which performance requirements are met) by setting the "carbon income" off against negative Net Present Values from timber yield.

Run Sobol_driver_yield.R and Sobol_driver_carb.R to calculate the 1st order and total Sobol' sensitivity indices, which indicate the relative impact of each uncertain factor on the objectives' variance. Sobol_driver_yield.R calculates these indices for the variance in Net Present Value of timber yield and Sobol_driver_carb.R for the variance in net carbon sequestration. The Sobol scripts source the objective_function_yield_vect.R and objective_function_carb_vect.R, respectively.




**Reference beech growth model:** *Trasobares, A., Zingg, A., Walthert, L., Bigler, C. (2016) A climate-sensitive empirical growth and yield model for forest management planning of even-aged beech stands. Eur J Forest Res 135(2), 263-282. DOI: 10.1007/s10342-015-0934-7.*

**Questions?** Naomi Radke (naomikradke@gmail.com)

