function [] = inditek_gammadiversity(ext_pattern)

disp('** inditek_gammadiversity.m **')

% calculate global (gamma) diversity by iteratively adding a % of diversity of each transect according to their minimum peak distance to to peaks with greater transect diversity (zig-zag integration)
if ext_pattern==3
    ini=3;
else
    ini=1;
end
list=dir('*transectsDiv.mat');
eval(['load ' list(1).name ' Point_timeslices']);

for i=1:length(list)
    i
    eval(['load ' list(i).name ' Point_timeslices LAT LON Dtransects']);
    
    gammaD=NaN(length(Point_timeslices),1);
    gammaD(1:ini)=1; %initial conditions
    
    for j=ini+1:length(Point_timeslices)
        
        lat=LAT{j};lon=LON{j};
        lat=lat(1,:); %peaks lat
        lon=lon(1,:); %peaks lon
        
        D=Dtransects{j}';
        
        distP=NaN(size(D,1),size(D,1)); % distance between peaks
        
        [d,s]=sort(D,'descend'); %sort transect diversity from greatest to lowest d= sorted transect diversity
        lat2=lat(s);lon2=lon(s); %reorder lat-lon peak position
        for l=1:length(D)
            for m=1:length(D)
                distP(l,m)=lldistkm([lat2(l),lon2(l)],[lat2(m),lon2(m)]);
            end
        end
        
        distP=triu(distP,+1);
        tf = tril(true(size(distP,1),size(distP,2)));
        distP(tf)=NaN; %keep only the upper part of the peak distance matrix and NaN the lower part
        
        J= 0.9356.*exp(-0.0024.*distP)+0.0645; % Jaccard according to fig8 Miller 2009 and assuming Dist as mean of great circle range
        
        Dg=d(1);
        for l=2:length(d)
            f=find(distP(1:l-1,l)==nanmin(distP(1:l-1,l))); %find the minimum distance between the peak to be integrated and the rest of peaks already integratedalready integrated
            if length(f)>1
                f=f(1); %if there are more than 1 peak with equal distance, we pick the one with the greatest transect diversity
            end
            v=(J(f,l).*(1+(max(d(l-1:l)))/min(d(l-1:l))))./(1+J(f,l)) ; % overlap index translated from Jaccard (J)
            if v>1 % often because sorted peak diversities may be very different in D and NN distance too short such that shared % is very high and thus translation from J to V would be >1 which is imposible!
                v=(J(f,l).*2)./(1+J(f,l)); %in these cases we assume d(l-1)=d(l) to avoid d(l)*(1-v)<1 and thus subtracting diversity instead of adding it
            end
            Dg=nansum([Dg,(d(l)*(1-v))]);
        end
        gammaD(j)=Dg; 
    end
    clear D Dg 
    a=list(i).name;
    f=strfind(a,'_');
    a=a(1:f(1));
    eval(['save ' a 'global gammaD Point_timeslices ext_pattern']);
end

return


