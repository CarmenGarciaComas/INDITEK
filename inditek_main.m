close all
clear all

addpath(genpath('~/toolbox'))

disp('Run this script in a folder without .mat files from old runs. If it is OK press any key to continue')
pause

% folder with public matlab tools used in our code:

addpath(genpath('~/toolbox/')) 

% ----------------------------------- PARAMETERS --------------------------

%CHOOSE model parameters:

Q10 = 1.75; %n.u.
kfood = 0.5; %[POC mol * m-2 yr-1]
rhomin = 0.001; %[MA-1]  
rhomax = 0.035; %[MA-1]
Kmax=16;% Carrying capacity of #genera at maximum food availability
Kmin=4; % Carrying capacity of #genera at minimum food availability

% CHOOSE extinction pattern:

% 1. Zaffos curve
% 2. Alroy curve
% 3. Sepkoski curve

ext_pattern=2;

% ----------------------------------- load data----------------------------
% Sea floor age:
load data/Point_ages_xyz
% Food & temperature:
load data/Point_foodtemp

% ----------------------------------- MODEL -------------------------------

% NET DIVERSIFICATION RATE:
[rho_ocean,rho_shelf,K_ocean,K_shelf,Rho_explain,Point_timeslices]=inditek_rhonet(rhomin,rhomax,kfood,Q10,Kmax,Kmin,food_ocean,food_shelf,temp_ocean,temp_shelf,ext_pattern,shelf_lonlatAge,ocean_lonlatAge,Point_timeslices); 

% ALPHA DIVERSITY (particles):
%logistic model
label_cases = {'logistic'}; 
inditek_alphadiv(label_cases,rho_ocean,rho_shelf,K_ocean,K_shelf,Point_timeslices,ext_pattern,shelf_lonlatAge,ocean_lonlatAge);
%exponential model
K_ocean=inf*ones(size(K_ocean));K_shelf=inf*ones(size(K_shelf)); % overwrite K_ocean & K_shelf as infinite for the exponential model to run
label_cases = {'exponential'}; 
inditek_alphadiv(label_cases,rho_ocean,rho_shelf,K_ocean,K_shelf,Point_timeslices,ext_pattern,shelf_lonlatAge,ocean_lonlatAge)
clear rho_ocean rho_shelf K_ocean K_shelf

% DIVERSITY BY GRID:
label_cases = {'exponential','logistic'}; 
warning('off')
inditek_gridding_alphadiv(label_cases,ext_pattern)
% GAMMA DIVERSITY (global):
inditek_peaks_troughs(label_cases,ext_pattern)
inditek_bresenham(ext_pattern)
inditek_diversity_transects(ext_pattern)
inditek_gammadiversity(ext_pattern)

% plot global diversity time series
inditek_plotgamma(ext_pattern)

return
