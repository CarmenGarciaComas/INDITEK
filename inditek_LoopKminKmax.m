


close all
clear all

disp('CHECK that in the folder CODE-RUNS there are no .mat files(from old runs; these should be oved to another folder). If it is OK press any key to continue')
pause

% folder with public matlab tools used in our code:

addpath(genpath('~/toolbox/'))

% ----------------------------------- PARAMETERS --------------------------

%CHOOSE model parameters:

Q10 = 1.75; %n.u.
kfood = 0.5; %[POC mol * m-2 yr-1]
rhomin = 0.001; %[MA-1]
rhomax = 0.035; %[MA-1]

% CHOOSE extinction pattern:

% 1. Zaffos curve
% 2. Alroy curve
% 3. Sepkoski curve

ext_pattern=3;

% ----------------------------------- load data----------------------------
% Sea floor age:
load data/Point_ages_xyz
% Food & temperature:
load data/Point_foodtemp

N=[2,4,8,16,32,64,128,256];
n=0;
CCC_5=NaN(length(N),length(N));

for Kmin=N(1:2)% Carrying capacity of #genera at maximum food availability
    f=find(Kmin==N);
    n=n+1
    m=0;
    for pos=f+1%:length(N)   
        if n==length(N)
            Kmax=N(f);
            m=n
        else
            Kmax=N(pos); % Carrying capacity of #genera at minimum food availability
            m=pos
        end  
        % NET DIVERSIFICATION RATE:
        [rho_ocean,rho_shelf,K_ocean,K_shelf,Rho_explain,Point_timeslices]=inditek_rhonet(rhomin,rhomax,kfood,Q10,Kmax,Kmin,food_ocean,food_shelf,temp_ocean,temp_shelf,ext_pattern,shelf_lonlatAge,ocean_lonlatAge,Point_timeslices); 
        % ALPHA DIVERSITY (particles):
        %logistic model
        label_cases = {'logistic'};
        inditek_alphadiv(label_cases,rho_ocean,rho_shelf,K_ocean,K_shelf,Point_timeslices,ext_pattern,shelf_lonlatAge,ocean_lonlatAge);
        clear rho_ocean rho_shelf K_ocean K_shelf
        % DIVERSITY BY GRID:
        warning('off')
        inditek_gridding_alphadiv(label_cases,ext_pattern)
        % GAMMA DIVERSITY (global):
        inditek_peaks_troughs(label_cases,ext_pattern)
        inditek_bresenham(ext_pattern)
        inditek_diversity_transects(ext_pattern)
        inditek_gammadiversity(ext_pattern)
        % Lin's Coefficient of concordance (CCC)
        CCC= inditek_Model_Fossil_CCC(Q10,kfood,rhomin,rhomax,Kmax,Kmin,ext_pattern);      
        CCC_5(m,n)=CCC;
        delete('*.mat')
    end
end

Param=[Q10;kfood;rhomin;rhomax];
Paramlegend={'Q10';'kfood';'rhomin';'rhomax'};
CCC=flip(CCC_5,1);

save inditek_LogisticModel_KminKmax CCC Param* ext_pattern

figure
imAlpha=ones(size(CCC));
imAlpha(isnan(CCC))=0;
imagesc(CCC,'AlphaData',imAlpha);
imagesc(1:size(CCC,1),1:size(CCC,2),CCC,'AlphaData',imAlpha)
set(gca,'Fontname','timesnewroman','fontsize',12)
a=xlabel('Kmin (# genera)');
set(a,'Fontname','timesnewroman','fontsize',14);
a=ylabel('Kmax (# genera)');
set(a,'Fontname','timesnewroman','fontsize',14);
a=title(['CCC'])
set(a,'Fontname','timesnewroman','fontsize',16);
set(gca, 'ytick',1:size(CCC,1),'yticklabel',{'2','4','8','16','32','64','128','256'},'xtick',1:size(CCC,1),'xticklabel',{'2','4','8','16','32','64','128','256'});
colorbar
colormap(jet(100))
caxis([0,1])
print -dpng -r200 inditek_KminKmax_CCC

return
