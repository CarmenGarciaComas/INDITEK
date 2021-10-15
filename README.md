# INDITEK

A global model of diversification (#genera My^-1) of marine benthic invertebrates in the Phanerozoic (from 541 Ma to present).

The model is written in MATLAB 2010b in a 1.7 GHz 2-core Intel Core i7, and tested with MATLAB 2021a in a MacOS 2.3 GHz 8-Core Intel Core i9, with MATLAB 2020b on Windows with a 2.5 GHz Intel i5-3210M, and on Linux Debian with a 2.6 GHz Intel Xeon E5-2640 processor.


## How to run the model:

The main module of INDITEK is **inditek_main**. To obtain the exponential and logistic models, run inditek_main. In inditek_main you can manually change model parameters and extinction pattern to apply and compare the fossil time series to the output. 

To run the model, you need to have in the same folder the functions (.m), the folder toolbox and the folder data. The folder toolbox contains public functions used by INDITEK. The folder data contain floor age data from the plate-tectonic/paleo-elevation model and food-temp data from the cGenie earth-system model.

### Model functions & outputs (outputs correspond to figure 1 & 2 in the manuscript):
The module inditek_main runs the following sequence of functions with their corresponding outputs:

- **inditek_alphadiv**: computes diversity in the model particules ---> *INDITEKexponential_alpha.mat & INDITEKlogistic_alpha.mat*
- **inditek_gridding_alphadiv**: interpolates diversity in 0.5x0.5º grids ---> *INDITEKexponential_grid.mat & INDITEKlogistic_grid.mat*
- **inditek_peaks_troughs**: detects diversity peaks & troughs ---> *INDITEKexponential_PT.mat & INDITEKlogistic_PT.mat*
- **inditek_bresenham**: transect diversity from peaks to troughs ---> *INDITEKexponential_transects.mat & INDITEKlogistic_transects.mat*
- **inditek_diversity_transects**: diversity integration along transects ---> *INDITEKexponential_transectsDiv.mat & INDITEKlogistic_transectdiv.mat*
- **inditek_gammadiversity**: computes global diversity from transect diversities ---> *INDITEKexponential_global.mat & INDITEKlogistic_global.mat*
- **inditek_plotgamma**: Compares the model diversity time series with the fossil time series ---> *INDITEK_ModelvsFossil.jpeg*

### Loop to calibrate the logistic model Kmin-Kmax (output corresponds to figure 3 in the manuscript):
The script to calculate the Lin’s concordance correlation coefficient (CCC) between the logistic model time series and the fossil time series with a Kmin-Kmax combinations in the geometric sequence of base 2, from 2 to 256 genera, is **inditek_LoopKminKmax**. Kmin-Kmax are the range of carrying capacities in the global maps according to food availability.

### Figure 1 to 3:
The script to plot the main manuscript figures 1-3 is **inditek_figures**. This needs to be in the same folder of INDITEK*.mat files generated with inditek_main.m & inditek_loopKminKmax.m
