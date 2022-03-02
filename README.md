# INDITEK

A global model of diversification (#genera My^-1) of marine invertebrates in the Phanerozoic (from 541 Ma to present).

This work was funded by national research grant CGL2017-91489-EXP (INDITEK project) from the Spanish government.

The model is written in MATLAB 2013b in a 1.7 GHz 2-core Intel Core i7, and tested with MATLAB 2021a in a MacOS 2.3 GHz 8-Core Intel Core i9, and with MATLAB 2020b on Windows with a 2.5 GHz Intel i5-3210M, and on Linux Debian with a 2.6 GHz Intel Core i9-9980HK processor.

## How to run the model:

The main module of INDITEK is **inditek_main.m**. To generate maps of the spatial distributions of diversity and time trajectories of global diversity based on the logistic and exponential diversification models, run inditek_main. In inditek_main you can manually change model parameters and extinction pattern to apply and compare the fossil time series to the output. 

To run the model, you need to have in the same folder the functions (.m), the folder toolbox and the folder data. The folder toolbox contains public functions used by INDITEK. The folder data contains:
- Point_ages_xyz.mat=floor age data from the plate-tectonic/paleo-elevation model.
- Point_foodtemp.mat=food-temp data from the cGenie earth-system model.
- landShelfOceanMask.mat= 0-2 mask to distingish land-shelf-ocean grids.
- POSoceanRidge.mat= 1 mask locating troughs to trace transects from peaks of diversity.
- LonDeg.mat= degrees of longitud according to the latitude with distance equivalent to 1 degree at the equator. This is used to search for nearest neighbour (NN) coastal points in a square to account for immigration.
- rhoExt.xlsx= mass extinction patterns to input in the model.
- FossilTimeSeries.xlsx= fossil time series to compare with the model time series.

### Model functions & outputs (outputs correspond to figure 1 & 2 in the manuscript):
The module inditek_main runs the following sequence of functions with their corresponding outputs:

- **inditek_rhonet.m**: calculates diversification rate (rho) and effective carrying capacity (Keff) ---> *particleRhoExpLog.mat*
- **inditek_alphadiv.m**: computes diversity in the model particles ---> *INDITEKexponential_alpha.mat & INDITEKlogistic_alpha.mat*
- **inditek_gridding_alphabyK.m**: interpolates logistic diversity/Keff in 0.5ºx0.5º grids ---> *INDITEKalphabyk_grid.mat
- **inditek_gridding_alphadiv.m**: interpolates diversity in 0.5ºx0.5º grids ---> *INDITEKexponential_grid.mat & INDITEKlogistic_grid.mat*
- **inditek_peaks_troughs.m**: detects diversity peaks & troughs ---> *INDITEKexponential_PT.mat & INDITEKlogistic_PT.mat*
- **inditek_bresenham.m**: traces diversity along transects ---> *INDITEKexponential_transects.mat & INDITEKlogistic_transects.mat*
- **inditek_diversity_transects.m**: diversity integration along transects ---> *INDITEKexponential_transectsDiv.mat & INDITEKlogistic_transectdiv.mat*
- **inditek_gammadiversity.m**: computes global diversity from transect diversities ---> *INDITEKexponential_global.mat & INDITEKlogistic_global.mat*
- **inditek_plotgamma.m**: Compares the model diversity time series with the fossil time series ---> *INDITEK_ModelvsFossil.jpeg*

### Loop to calibrate the logistic model Kmin-Kmax (output corresponds to extended data figure 6 in the manuscript):
To calibrate the Kmin and Kmax parameters of the logistic model, use the function **inditek_LoopKminKmax.m**. This function runs simulations of pair-wise Kmin and Kmax combinations in a geometric sequence of base 2, from 2 to 256 genera. Then, the function tests the effect of changing the Kmin and Kmax values on the Lin’s concordance correlation coefficient (CCC) calculated with **inditek_Model_Fossil_CCC.m** between the normalized diversities generated by the model and those estimated from the fossil record. The outputs of this function are *inditek_LogisticModel_KminKmax.mat* & *inditek_KminKmax_CCC.jpeg*.

To get the outputs of the calibrated logistic model, set the Kmin-Kmax combination in **inditek_main.m** and relaunch it (make sure you have moved the .mat files from the first launch in a separate folder). Rename the INDITEKlogistic_global.mat output file as INDITEKcaliblogistic_global.mat , add the INDITEKlogistic_global.mat file from the first run to the main folder, and plot results with **inditek_figures** .

### Figures 1 to 3:
The script **inditek_figures.m** plots the main results of the manuscript: (1) the time series of the logistic, exponential and calibrated model; (2) the calibrated model maps; (3) The diversity/carrying capacity maps of the calibrated model. It needs to be in the same folder of INDITEK*.mat files generated with inditek_main.m. For figure 2 and 3 you further need to download the free MATLAB package m_map (https://www.eoas.ubc.ca/~rich/map.html) and add it to the toolbox folder.

Explanations of all the functions are written inside them. For any further doubt, do not hesitate to contact me: carmencomas@gmail.com or cgcomas@icm.csic.es
