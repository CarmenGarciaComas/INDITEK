
% NET DIVERSIFICATION RATE (Rho)
% inditek_rhonet produces net diversification rate (Rho) with Qlim = Qfood*Qtemp
% Qfood given by the growth saturation curve & Qtemp given by the Eppley curve

function [rho_ocean,rho_shelf,K_ocean,K_shelf,Rho_explain,Point_timeslices] = inditek_rhonet(rhomin,rhomax,kfood,Q10,Kmax,Kmin,food_ocean,food_shelf,temp_ocean,temp_shelf,ext_pattern,shelf_lonlatAge,ocean_lonlatAge,Point_timeslices)

disp('** inditek_rhonet.m **')

% Rho in moments of extinction to substitute rho in those time slices where extinction has taken place
[D,T] = xlsread('data/rhoExt.xlsx');
rhoExt=D(:,2:end);
clear D T
rhoExt=rhoExt(:,ext_pattern);

rhoExtocean=repmat(rhoExt,1,size(ocean_lonlatAge,1));
rhoExtocean=rhoExtocean';
rhoExtshelf=repmat(rhoExt,1,size(shelf_lonlatAge,1));
rhoExtshelf=rhoExtshelf';

% Upper (M) and lower (m) limits of food as the 99 & 1 of the whole time series
a=[food_ocean(isnan(food_ocean)==0);food_shelf(isnan(food_shelf)==0)];
Mfood=quantile(a,0.99);
mfood=quantile(a,0.01);

% Effective carrying capacity: max N of genera supported according to food;
K_shelf=Kmax - (Kmax-Kmin)*((Mfood-food_shelf)./(Mfood-mfood));
K_ocean=Kmax - (Kmax-Kmin)*((Mfood-food_ocean)./(Mfood-mfood));
% bounded between Kmax & Kmin
K_shelf(K_shelf>Kmax)=Kmax;
K_shelf(K_shelf<Kmin)=Kmin;
K_ocean(K_ocean>Kmax)=Kmax;
K_ocean(K_ocean<Kmin)=Kmin;

% Rho without extinctions
        
rho_ocean1=NaN(size(food_ocean));
rho_shelf1=NaN(size(food_shelf));
    
    for i=1:size(food_ocean,2)
        % Upper (M) and lower (m) limits of temperature as the 99 & 1 percentiles in each time frame
        a=[temp_ocean(:,i);temp_shelf(:,i)];
        Mtemp=quantile(a(isnan(a)==0),0.99);
        mtemp=quantile(a(isnan(a)==0),0.01);
        %ocean particles 
        Qfood = food_ocean(:,i) ./ (kfood + food_ocean(:,i)); % food limitation
        
        EppleyCurve = Q10.^((temp_ocean(:,i) - mtemp)/10);
        EppleyCurve_max = Q10^((Mtemp - mtemp) / 10); %normalization factor for the Eppley curve.
        EppleyCurve = EppleyCurve / EppleyCurve_max;
        Qtemp = EppleyCurve; % temperature limitation
        
        Qtemp(Qtemp>1)=1;
        Qfood(Qfood>1)=1;
        Qlim = (Qfood .* Qtemp); % interactive food-temp. limitation
        
        rho_ocean1(:,i) = rhomax - (rhomax-rhomin) * (1.0 - Qlim);
        %shelf particles
        Qfood = food_shelf(:,i) ./ (kfood + food_shelf(:,i)); % food limitation
        
        EppleyCurve = Q10.^((temp_shelf(:,i) - mtemp)/10);
        EppleyCurve_max = Q10^((Mtemp - mtemp) / 10); %normalization factor for the Eppley curve.
        EppleyCurve = EppleyCurve / EppleyCurve_max;
        Qtemp = EppleyCurve; % temperature limitation
        
        Qtemp(Qtemp>1)=1;
        Qfood(Qfood>1)=1;
        Qlim = (Qfood .* Qtemp); % interactive food-temp. limitation
        
        rho_shelf1(:,i) = rhomax - (rhomax-rhomin) * (1.0 - Qlim);
    end   

rho_ocean1(rho_ocean1>rhomax)=rhomax;
rho_ocean1(rho_ocean1<rhomin)=rhomin;
rho_shelf1(rho_shelf1>rhomax)=rhomax;
rho_shelf1(rho_shelf1<rhomin)=rhomin;

rho_ocean=rhoExtocean;
rho_shelf=rhoExtshelf;
all_timeslices = fliplr([Point_timeslices(end):1:Point_timeslices(1)]);
posPT=NaN(length(Point_timeslices),1);

% save rho for the 82 time frames at their corresponding position (posPT) in the -541MA:-1MA:0MA frames
for i=1:length(Point_timeslices)
    a=find(all_timeslices == Point_timeslices(i));  
    f=find(rho_ocean(:,a)<0);
    if isempty(f) %time frames of the 82 without extinction
        rho_ocean(:,a)=rho_ocean1(:,i);
        rho_shelf(:,a)=rho_shelf1(:,i);
        posPT(i)=a;
    end
end
posPT(isnan(posPT))=[];
% find the gaps without extinction between -541MA to 0MA to fill with copies of next point_timeslice (posPT)
% inditek_alphadiv will then run during the gap with the configuration of the next time frame to obtain local diversities at that time frame
f=find(rho_ocean(1,:)==0.01); 

for i=1:length(posPT)-1
    rho_ocean(:,f(f>posPT(i) & f<posPT(i+1)))=repmat(rho_ocean(:,posPT(i+1)),1,length(f(f>posPT(i) & f<posPT(i+1))));
    rho_shelf(:,f(f>posPT(i) & f<posPT(i+1)))=repmat(rho_shelf(:,posPT(i+1)),1,length(f(f>posPT(i) & f<posPT(i+1))));
end

Rho_explain=['Model with extinctions[rho=' num2str(rhomin) '-' num2str(rhomax) '; food saturation curve, kc=,' num2str(kfood) ';temperature Eppley curve, Q10=' num2str(Q10) ']'];

save particleRhoExpLog rho_ocean rho_shelf K_ocean K_shelf Rho_explain Point_timeslices -v7.3
return


