function [] = inditek_diversity_transects(ext_pattern)

disp('** inditek_diversity_transects.m **')

% To calculate diversity along transects cut at 555km length
if ext_pattern==3
    ini=3;
else
    ini=1;
end
list=dir('*transects.mat');

eval(['load ' list(1).name ' Point_timeslices']);


for i=1:length(list)
    
    eval(['load ' list(i).name]);
    
    Dtransects=cell(length(Point_timeslices),1); %diversity along transect   
    for j=ini+1:length(Point_timeslices)
        
        div=DIV{j}; %diversity along transects (columns) from peak to trough
        lat=LAT{j};lon=LON{j}; % position of grids crossed by transects from peak to trough
        
        dist=NaN(size(div));
        V=dist; % overlap index
        V(1,:)=0; %to initialise 100% of shared diversity with Div at the peak D(1) (D*(1-0)=D)
        divV=dist;
        for k=1:size(div,2)
            
            L=max(find(isnan(div(:,k))==0)); % to perform loop until the length of each transect
            
            if L>1 % for a few peaks I did not find a near trough without passing most through land and thus there is not transect to go through
                
                for l=2:L %to 'cut' transects at 555km length
                    dist(l,k)=lldistkm([lat(1,k),lon(1,k)],[lat(l,k),lon(l,k)]);
                    if dist(l,k)>=555
                        break
                    end
                end
                
                L=max(find(isnan(dist(:,k))==0)); % transects are now of 555 km or less
                dist=NaN(size(div));
                for l=2:L
                    dist(l,k)=lldistkm([lat(l-1,k),lon(l-1,k)],[lat(l,k),lon(l,k)]);   %pair-wise distance of contiguous diversity along the transect alpha(l-1), alpha(l)
                    J= 0.9356.*exp(-0.0024.*dist(l,k))+0.0645; % Jaccard (J) according to fig8 Miller 2009 and assuming Dist as mean of great circle range
                    V(l,k)=(J.*(1+(max(div(l-1:l,k))/min(div(l-1:l,k)))))/(1+J); % from J to overlap coeff. (V)
                    if V(l,k)>1 % case in which contiguous diversities have such a big diversity difference that V>1 which is absurd (but can happen because J to V is generalised to a single shape, then we assume that alpha(l-1) & alpha(l)are the same to avoid substracting diversity instead of adding it (1-V<0)
                        V(l,k)=(J.*2)/(1+J); %we assume max(alpha(l-1,l))/min(alpha(l-1,l))=1
                    end
                end
            end
        end
        divV=div.*(1-V); %matrix with proportions of alpha(l) different from the previous diversities in the cut transects
        divV(isnan(divV)==1)=0;
        divV=sum(divV);
        Dtransects{j}=divV; % diversities added along the cut transects
    end
    
    a=list(i).name;
    a=a(1:end-4);
    
    eval(['save '  a 'Div Point_timeslices X Y LAT LON DIV Dtransects']);
end

