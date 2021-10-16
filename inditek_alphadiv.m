function [] = inditek_alphadiv(label_cases,rho_ocean,rho_shelf,K_ocean,K_shelf,Point_timeslices,ext_pattern,shelf_lonlatAge,ocean_lonlatAge)

disp('** inditek_alphadiv.m **')

% New ocean crust Do=1
% New coastal particles or inundated land Do= nearest neighbour (NN) D
% Dt=Do*exp^(rho*t);
% Dt+1=rho*Dt*(1-Dt/K); % exponential model= K=inf; logistic model = K (dependent on food)
% Because time steps are >1MA and I want to calculate to the 1MA step, I add deltaAge=time elapsed (this can be different from time frame n & n+1 difference because some particules
% were created sometime within the XMA time)

pt=Point_timeslices+1;
pt=fliplr(pt);

Do=1; % initialise diversity at time -541 MA with #1 genus area^-^1 (and for new ocean crust formed in the mid-ocean ridges

load data/LonDeg.mat % longitude width search according to the latitude for a Lat by Lon window to search for nearest neighbour (NN) particles to initiate D at continental margins just submerged
latWindow=2.5; lonWindow=2.5; % to shorten computation time, we search for NN particles in a window of length equivalent to 2.5 by 2.5 degrees at the equator
olab=0; %label for the few newly submerged coastal particles that have no nearby submerged coastal particles to receive inmigrated diversity (we get ocean NN diversity instead)

D_ocean=NaN(size(ocean_lonlatAge,1),542);
D_shelf=NaN(size(shelf_lonlatAge,1),542);

% positions to work at 1MA pace in >1MA solved pace (about 5 MA)
count=0; % time frame resolved (MA) (there are 82 out of 542MA defined by the Point_timeslices)
step=0;  % 82 time frames (steps in the loop)

ts2=Point_timeslices(1)+1; % next Point_timeslice (to fill the gap between both at each loop)

%ini:time frame from which we start accumulating diversity, dependent on which fossil diversity we are using to fit extinctions and compare pattern
if ext_pattern==3
    ini=4;
else
    ini=2;
end


for ts=Point_timeslices%current Point_timeslice
    step=step+1;
    count=count+(ts2-ts);
    
    ageO=ocean_lonlatAge(:,step,3);
    posO=find(isnan(ageO)==0); % ocean particles "active"
    
    ageS=shelf_lonlatAge(:,step,3);
    posS=find(isnan(ageS)==0 & ageS>0); %   % shelf particles "active" ; 0s correspond to land above sea water (so these are considered inactive D=NaN)
    
   
    
    if ts>Point_timeslices(ini) %initialise diversity
        D_ocean(posO,count)=Do;
        D_shelf(posS,count)=Do;  
    else
        
        %######################## New ocean particules (start D=Do)
        
        pos2O=find(ageO(posO)==0); % New particle in the ocean (newly formed crust have age=0 & thus Do=1)
        if length(unique(rho_ocean(:,count)))==1 && unique(rho_ocean(:,count))<0   % extinction period
            D_ocean(posO(pos2O),count)=Do+(rho_ocean(posO(pos2O),count))*Do;  % starting diversity in a moment of extinction is affected by extinction
        else
            D_ocean(posO(pos2O),count)=Do;
        end
        
        
        deltaAgeO=ageO(posO)-ocean_lonlatAge(posO,step-1,3); %particle time active within the time frame
        
        %######################## Ocean particules did not exist in time-1 and now start accumulating diversity
        
        pos2O=find(isnan(deltaAgeO)==1 & ageO(posO)>0); % particle did not exist in time count2 & it has appeared sometime within the time gap (count2:count)
        if isempty(pos2O)==0
            
            %They may appear at any moment of the time gap (nma)
            
            nma=unique(ageO(posO(pos2O)));
            for k=1:length(nma)
                pos=find(ageO(posO(pos2O))==nma(k));
                D_ocean(posO(pos2O(pos)),count-nma(k)+1)=Do+rho_ocean(posO(pos2O(pos)),count-nma(k)+1)*Do;  %accumulating diversity from the moment of appearance in the gap+1
                for gap=count-nma(k)+2:count
                    if length(unique(rho_ocean(:,gap)))==1 && unique(rho_ocean(:,gap))<0   % extinction period
                        D_ocean(posO(pos2O(pos)),gap)=D_ocean(posO(pos2O(pos)),gap-1)+D_ocean(posO(pos2O(pos)),gap-1).*rho_ocean(posO(pos2O(pos)),gap);
                    else
                        D_ocean(posO(pos2O(pos)),gap)=D_ocean(posO(pos2O(pos)),gap-1)+D_ocean(posO(pos2O(pos)),gap-1).*rho_ocean(posO(pos2O(pos)),gap).*(1-(D_ocean(posO(pos2O(pos)),gap-1)./K_ocean(posO(pos2O(pos)),step)));
                    end
                end
            end
        end
    
        deltaAgeS=ageS(posS)-shelf_lonlatAge(posS,step-1,3);% a vector of 5MA & some NaN if the plate did not exist at time count-1
        
        %######################## Continental shelf particules did not exist or were not inundated in time-1 and now start accumulating diversity
        
        pos2S=find(isnan(deltaAgeS)==1 & ageS(posS)<=ts2-ts); % particles in time t-1 that did not exist
        pos2S=[pos2S;find(shelf_lonlatAge(posS,step-1,3)==0 & ageS(posS)<=ts2-ts)];  % add also particles in time t-1 that were above land
        
        if isempty(pos2S)==0
            for k=1:length(pos2S)
                % find particles with diversity to 'import' to initialise diversity in the new shelf particles
                f=find(abs(abs(shelf_lonlatAge(posS(pos2S(k)),step,2))-LonDeg(:,1))==min(abs(abs(shelf_lonlatAge(posS(pos2S(k)),step,2))-LonDeg(:,1))));
                lon=lonWindow*LonDeg(f,2);
                lim=find(shelf_lonlatAge(posS,step,1)<=shelf_lonlatAge(posS(pos2S(k)),step,1)+lon & shelf_lonlatAge(posS,step,1)>=shelf_lonlatAge(posS(pos2S(k)),step,1)-lon & shelf_lonlatAge(posS,step,2)<=shelf_lonlatAge(posS(pos2S(k)),step,2)+latWindow  & shelf_lonlatAge(posS,step,2)>=shelf_lonlatAge(posS(pos2S(k)),step,2)-latWindow );
                f=find(D_shelf(posS(lim),count2)>0);
                lim=lim(f);
                if isempty(f)
                    lim=find(ocean_lonlatAge(posO,step,1)<=shelf_lonlatAge(posS(pos2S(k)),step,1)+lon & ocean_lonlatAge(posO,step,1)>=shelf_lonlatAge(posS(pos2S(k)),step,1)-lon & ocean_lonlatAge(posO,step,2)<=shelf_lonlatAge(posS(pos2S(k)),step,2)+latWindow  & ocean_lonlatAge(posO,step,2)>=shelf_lonlatAge(posS(pos2S(k)),step,2)-latWindow );
                    f=find(D_ocean(posO(lim),count2)>0);
                    lim=lim(f);
                    olab=1; % if no close shelf particle with diversity accumulated (inside the 2x2deg. window), we assume inmigration from nearest ocean particle in the 2x2deg window.
                end
                if isempty(f)
                    d=Do; %if no particle with accumulated D in the fixed area (Lat*~Lon window), I set d to Do (1) (initiate diversity as in ocean ridges)
                else
                    dist=NaN(length(lim),1);
                    if olab==1; % set initial diversity (d) to nearest ocean particle in the LatxLon window
                        olab=0; %reset tag for next time frame
                        for l=1:length(lim)
                            dist(l)=lldistkm([shelf_lonlatAge(posS(pos2S(k)),step,2),shelf_lonlatAge(posS(pos2S(k)),step,1)],[ocean_lonlatAge(posO(lim(l)),step,2),ocean_lonlatAge(posO(lim(l)),step,1)]);
                        end
                        lim(dist==0)=[];
                        dist(dist==0)=[];
                        f=find(dist==min(dist));
                        lim=lim(f);
                        d=nanmean(D_ocean(posO(lim),count2));
                        
                    else    % set initial diversity (d) to nearest shelf particle in the latxLon window
                        for l=1:length(lim)
                            dist(l)=lldistkm([shelf_lonlatAge(posS(pos2S(k)),step,2),shelf_lonlatAge(posS(pos2S(k)),step,1)],[shelf_lonlatAge(posS(lim(l)),step,2),shelf_lonlatAge(posS(lim(l)),step,1)]);
                        end
                        lim(dist==0)=[];
                        dist(dist==0)=[];
                        f=find(dist==min(dist));
                        lim=lim(f);
                        d=nanmean(D_shelf(posS(lim),count2));
                    end
                end
                % accumulate diversity during the time gap 
                if d>K_shelf(posS(pos2S(k)),step)
                    d=K_shelf(posS(pos2S(k)),step); %force local extinction if imported diversity is greater than K (greater than hat the location can support according to food availability).
                end
                if length(unique(rho_shelf(:,count2+1)))==1 && unique(rho_shelf(:,count2+1))<0 % extinction period
                    D_shelf(posS(pos2S(k)),count2+1)=d+rho_shelf(posS(pos2S(k)),count2+1)*d; 
                else
                    D_shelf(posS(pos2S(k)),count2+1)=d+rho_shelf(posS(pos2S(k)),count2+1)*d.*(1-(d/K_shelf(posS(pos2S(k)),step)));
                end
                
                for gap=count2+2:count
                    if length(unique(rho_shelf(:,gap)))==1 && unique(rho_shelf(:,gap))<0 % extinction period
                        D_shelf(posS(pos2S(k)),gap)=D_shelf(posS(pos2S(k)),gap-1)+D_shelf(posS(pos2S(k)),gap-1).*rho_shelf(posS(pos2S(k)),gap);
                    else
                        D_shelf(posS(pos2S(k)),gap)=D_shelf(posS(pos2S(k)),gap-1)+D_shelf(posS(pos2S(k)),gap-1).*rho_shelf(posS(pos2S(k)),gap).*(1-(D_shelf(posS(pos2S(k)),gap-1)./K_shelf(posS(pos2S(k)),step))).*((ts2-ts)./length(count2+2:count));
                    end
                end
                if  D_shelf(posS(pos2S(k)),gap)>K_shelf(posS(pos2S(k)),step)
                    D_shelf(posS(pos2S(k)),gap)=K_shelf(posS(pos2S(k)),step);
                end
            end
        end
        
        %######################## Special case of continental shelf particules that did not exist in time-1 and were artificially added in the Gplate model to fill gaps with age of nearest neighbour continental-shelf particles (thus we start diversity with nearest continental shelf particles)
        
        pos2S=find(isnan(deltaAgeS)==1 & ageS(posS)>ts2-ts); % particle in time t-1 did not exist or was above land and suddenly they have age greater than the time gap (ts2-ts)
        pos2S=[pos2S;find(shelf_lonlatAge(posS,step-1,3)==0 & ageS(posS)>ts2-ts)];
         
        if isempty(pos2S)==0
            for k=1:length(pos2S)
                % I search the nearest neighbour by iteratively increasing the searching window to shorten the computation time
                lim=find(shelf_lonlatAge(posS,step,1)<=shelf_lonlatAge(posS(pos2S(k)),step,1)+5 & shelf_lonlatAge(posS,step,1)>=shelf_lonlatAge(posS(pos2S(k)),step,1)-5 & shelf_lonlatAge(posS,step,2)<=shelf_lonlatAge(posS(pos2S(k)),step,2)+5  & shelf_lonlatAge(posS,step,2)>=shelf_lonlatAge(posS(pos2S(k)),step,2)-5 );
                f=find(D_shelf(posS(lim),count2)>0);
                if isempty(f) 
                    lim=find(shelf_lonlatAge(posS,step,1)<=shelf_lonlatAge(posS(pos2S(k)),step,1)+10 & shelf_lonlatAge(posS,step,1)>=shelf_lonlatAge(posS(pos2S(k)),step,1)-10 & shelf_lonlatAge(posS,step,2)<=shelf_lonlatAge(posS(pos2S(k)),step,2)+10  & shelf_lonlatAge(posS,step,2)>=shelf_lonlatAge(posS(pos2S(k)),step,2)-10 );
                    f=find(D_shelf(posS(lim),count2)>0);
                end
                if isempty(f)
                    lim=find(shelf_lonlatAge(posS,step,1)<=shelf_lonlatAge(posS(pos2S(k)),step,1)+15 & shelf_lonlatAge(posS,step,1)>=shelf_lonlatAge(posS(pos2S(k)),step,1)-15 & shelf_lonlatAge(posS,step,2)<=shelf_lonlatAge(posS(pos2S(k)),step,2)+15  & shelf_lonlatAge(posS,step,2)>=shelf_lonlatAge(posS(pos2S(k)),step,2)-15 );
                    f=find(D_shelf(posS(lim),count2)>0);
                end
                if isempty(f)
                    lim=find(shelf_lonlatAge(posS,step,1)<=shelf_lonlatAge(posS(pos2S(k)),step,1)+30 & shelf_lonlatAge(posS,step,1)>=shelf_lonlatAge(posS(pos2S(k)),step,1)-30 & shelf_lonlatAge(posS,step,2)<=shelf_lonlatAge(posS(pos2S(k)),step,2)+30  & shelf_lonlatAge(posS,step,2)>=shelf_lonlatAge(posS(pos2S(k)),step,2)-30 );
                    f=find(D_shelf(posS(lim),count2)>0);
                end
                lim=lim(f);
                dist=NaN(length(lim),1);
                for l=1:length(lim)
                    dist(l)=lldistkm([shelf_lonlatAge(posS(pos2S(k)),step,2),shelf_lonlatAge(posS(pos2S(k)),step,1)],[shelf_lonlatAge(posS(lim(l)),step,2),shelf_lonlatAge(posS(lim(l)),step,1)]);
                end
                lim(dist==0)=[];
                dist(dist==0)=[];
                f=find(dist==min(dist));
                lim1=lim(f);
                d=nanmean(D_shelf(posS(lim1),count2));
                if d>K_shelf(posS(pos2S(k)),step)
                    d=K_shelf(posS(pos2S(k)),step); %force local extinction if imported diversity is greater than K.
                end    
                % accumulate diversity during the time gap 
                if length(unique(rho_shelf(:,count2+1)))==1 && unique(rho_shelf(:,count2+1))<0 % extinction period 
                    D_shelf(posS(pos2S(k)),count2+1)=d+rho_shelf(posS(pos2S(k)),count2+1)*d;
                else 
                    D_shelf(posS(pos2S(k)),count2+1)=d+rho_shelf(posS(pos2S(k)),count2+1)*d.*(1-(d./K_shelf(posS(pos2S(k)),step)));
                end   
                for gap=count2+2:count
                    if length(unique(rho_shelf(:,gap)))==1 && unique(rho_shelf(:,gap))<0
                        D_shelf(posS(pos2S(k)),gap)=D_shelf(posS(pos2S(k)),gap-1)+D_shelf(posS(pos2S(k)),gap-1).*rho_shelf(posS(pos2S(k)),gap);
                    else
                        D_shelf(posS(pos2S(k)),gap)=D_shelf(posS(pos2S(k)),gap-1)+D_shelf(posS(pos2S(k)),gap-1).*rho_shelf(posS(pos2S(k)),gap).*(1-(D_shelf(posS(pos2S(k)),gap-1)./K_shelf(posS(pos2S(k)),step)));
                    end
                end
                if D_shelf(posS(pos2S(k)),gap)>K_shelf(posS(pos2S(k)),step)
                    D_shelf(posS(pos2S(k)),gap)=K_shelf(posS(pos2S(k)),step); %force local extinction if imported diversity is greater than K.
                end
            end
        end
        
        %######################## Particules continuing to accumulate diversity
        
        pos2O=find(deltaAgeO>0); %existing ocean particles continue to accumulate diversity
        pos2S=find(deltaAgeS>0 & deltaAgeS<=ts2-ts & shelf_lonlatAge(posS,step-1,3)~=0); %exisiting shelf particles with normal behaviour continue to accumulate diversity
        
        for gap=count2+1:count 
            if D_ocean(posO(pos2O),gap-1)>K_ocean(posO(pos2O),step)
                    d=K_ocean(posO(pos2O),step); %force local extinction if D>K
                else
                    d=D_ocean(posO(pos2O),gap-1);
            end
            if length(unique(rho_ocean(:,gap)))==1 && unique(rho_ocean(:,gap))<0 % extinction period
                D_ocean(posO(pos2O),gap)=d+d.*rho_ocean(posO(pos2O),gap).*(deltaAgeO(pos2O)./length(count2+1:count));
            else   
                D_ocean(posO(pos2O),gap)=d+d.*rho_ocean(posO(pos2O),gap).*(1-(d./K_ocean(posO(pos2O),step))).*(deltaAgeO(pos2O)./length(count2+1:count)); 
            end
            if D_shelf(posS(pos2S),gap-1)>K_shelf(posS(pos2S),step)
                d=K_shelf(posS(pos2S),step); %force local extinction if D> K.
            else
                d=D_shelf(posS(pos2S),gap-1);
            end
            if length(unique(rho_shelf(:,gap)))==1 && unique(rho_shelf(:,gap))<0 % extinction period
                D_shelf(posS(pos2S),gap)=d+d.*rho_shelf(posS(pos2S),gap).*(deltaAgeS(pos2S)./length(count2+1:count)) ;  
            else
                D_shelf(posS(pos2S),gap)=d+d.*rho_shelf(posS(pos2S),gap).*(1-(d./K_shelf(posS(pos2S),step))).*(deltaAgeS(pos2S)./length(count2+1:count)) ;  
            end
            
            if D_shelf(posS(pos2S),gap)>K_shelf(posS(pos2S),step)
                D_shelf(posS(pos2S),gap)=K_shelf(posS(pos2S),step);
            end
            if D_ocean(posO(pos2O),gap)>K_ocean(posO(pos2O),step)
                D_ocean(posO(pos2O),gap)=K_ocean(posO(pos2O),step);
            end
            
        end
    end
    ts2=ts;
    count2=count;
end

%flip to order from point time slice 1 (0MA) to 542 (541MA) 
D_ocean=flip(D_ocean,2);
D_shelf=flip(D_shelf,2);
%get the 82 Point time slices for which the model is resolved
D_ocean=D_ocean(:,pt);
D_shelf=D_shelf(:,pt);
%flip back once the point time slices are compiled
D_ocean=flip(D_ocean,2);
D_shelf=flip(D_shelf,2);

clear d pos*

eval(['save INDITEK' char(label_cases)  '_alpha D_ocean D_shelf Point_timeslices shelf_lonlatAge ocean_lonlatAge -v7.3']);

return

