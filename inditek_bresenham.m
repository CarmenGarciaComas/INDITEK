function [] = inditek_bresenham(ext_pattern)

disp('** inditek_bresenham.m **')

% 1# with the output from extrema2 (peaks) and location of ocean ridges (troughs), compute the lines from P to the nearest T with the bresenham algorithm
% 2# trace the grids that a line connecting P & T will pass
% 3# record the diversity of the grids along it

load data/LonDeg.mat % degrees to search for a 1 degree lon at equator
load data/landShelfOceanMask

landShelfOceanMask=flip(landShelfOceanMask,3);
if ext_pattern==3
    ini=3;
else
    ini=1;
end
list=dir('*_PT.mat');

for i=1:length(list)
    
    eval(['load ' list(i).name]);
    
    % 1. find nearest trough for each peak that does not have more than 20% of grids with land
    % 2. trace the grids that a line connecting P & T will pass
    % 3. record the diversity of the grids along it
    
    LAT=cell(size(PosMAX,1),1);
    LON=LAT;
    DIV=LAT;   
   
    for j=ini+1:size(PosMAX,1)
        
        lsoM=landShelfOceanMask(:,:,j);
        lsoM=lsoM';
        llP=[Y(PosMAX{j}),X(PosMAX{j})]; % peak positions
        llT=[Y(PosMIN{j}),X(PosMIN{j})]; % trough positions (ocean ridges identified as ocean grids with D<=1)
        [Pr,Pc] = ind2sub(size(Y),PosMAX{j}); %convert peak positions to row-col position for the bresenham imput
        [Tr,Tc] = ind2sub(size(Y),PosMIN{j}); %convert trough positions to row-col position for the bresenham imput
        
        d=D(:,:,j);
        
        % delete the few peaks <-75 latitude (too complicated to find troughs nearby in peaks inside the Antartica land mass, inbetween land)
        f=find(llP(:,1)<-75);
        llP(f,:)=[];
        Pr(f)=[];Pc(f)=[];
        
        lat=NaN(360,size(llP,1)); %maximum length a transect could have: case of longest longitudinal line at the poles (once corrected if direction is wrong)
        lon=lat;
        div=lat;
        
        for k=1:length(llP)
            
            f=find(abs(abs(llP(k,1))-LonDeg(:,1))==min(abs(abs(llP(k,1))-LonDeg(:,1)))); % which lat are we in? to determine size of search window
            lon2=20*LonDeg(f,2);
            
            if lon2>360 %the whole circunference (cases near the poles)
                lon2=360;
            end
            %search in a rectangle of area equivalent to 40Lon*20Lat deg at the equator
            %search window of longer longitude especially when window increases (following the distribution of peaks and troughs in the maps)
            lim=find(llT(:,2)<=llP(k,2)+lon2 & llT(:,2)>=llP(k,2)-lon2 & llT(:,1)<=llP(k,1)+10 & llT(:,1)>=llP(k,1)-10);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isempty(lim)==0
                dist=NaN(length(lim),1);
                posdel=zeros(length(lim),1);
                posdel2=zeros(length(lim),1);
                L=posdel;
                for l=1:length(lim)
                    dist(l)=lldistkm(llP(k,:),llT(lim(l),:));
                    [x,y]=bresenham(Pr(k),Pc(k),Tr(lim(l)),Tc(lim(l))); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                    pos=sub2ind(size(X),x,y);
                    f2=find(lsoM(pos)==2); %land grids crossed by the line
                    ll=length(x);
                    if length(f2)/length(x)>=0.2 % if 20% of grids are land
                        posdel(l)=1; %delete
                    end
                    % Compute the line on the other direction
                    if Pc(k)<Tc(lim(l))
                        cP=Pc(k)+(size(Y,2)-Tc(lim(l)));
                        cT=1;
                        [x,y]=bresenham(Pr(k),cP,Tr(lim(l)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                        y=y-(size(Y,2)-Tc(lim(l)));
                        y(y<=0)=size(Y,2)+y(y<=0);
                        pos=sub2ind(size(X),x,y);
                        f2=find(lsoM(pos)==2);
                        if length(f2)/length(x)>=0.2; % if 20% of grids are land
                            posdel2(l)=1;
                        end
                        
                        if length(x)<ll
                            L(l)=1; % transects shorter than the left to right direction
                        end
                     elseif Pc(k)>Tc(lim(l))
                        cT=Tc(lim(l))+(size(Y,2)-Pc(k));
                        cP=1;
                        [x,y]=bresenham(Pr(k),cP,Tr(lim(l)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                        y=y-(size(Y,2)-Pc(k));
                        y(y<=0)=size(Y,2)+y(y<=0);
                        pos=sub2ind(size(X),x,y);
                        f2=find(lsoM(pos)==2);
                        if length(f2)/length(x)>=0.2; % if 20% of grids are land
                            posdel2(l)=1;
                        end
                        
                        if length(x)<ll
                            L(l)=1; % transects shorter than the left to right direction
                        end
                    end
                end
                
                f1=find(L==0);
                f=find(posdel(f1)==1);
                f1=f1(f);
                
                f2=find(L==1);
                f=find(posdel2(f2)==1);
                f2=f2(f);
                
                L(unique([f1;f2]))=[];
                lim(unique([f1;f2]))=[];
                dist(unique([f1;f2]))=[];
                
                if isempty(dist)==0
                    % shortest distance and find closest with no big land mass in between
                    f=find(dist==min(dist));
                    if length(f)>1
                        f=f(1);
                    end
                    if L(f)==1
                        if Pc(k)<Tc(lim(f))
                            cP=Pc(k)+(size(Y,2)-Tc(lim(f)));
                            cT=1;
                            [x,y]=bresenham(Pr(k),cP,Tr(lim(f)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                            y=y-(size(Y,2)-Tc(lim(f)));
                            y(y<=0)=size(Y,2)+y(y<=0);
                        elseif Pc(k)>Tc(lim(f))
                            cT=Tc(lim(f))+(size(Y,2)-Pc(k));
                            cP=1;
                            [x,y]=bresenham(Pr(k),cP,Tr(lim(f)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                            y=y-(size(Y,2)-Pc(k));
                            y(y<=0)=size(Y,2)+y(y<=0);
                        end
                    elseif L(f)==0
                        [x,y]=bresenham(Pr(k),Pc(k),Tr(lim(f)),Tc(lim(f))); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                    end
                    
                    Y(x,y);lat(1:length(x),k)=ans(:,1);
                    X(x,y);lon(1:length(x),k)=ans(1,:)';
                    d(x,y); div(1:length(x),k)=diag(ans)';
                    
                else %in case there is no T with transect not crossing >20% land grids, I increase the search window and do the same procedure
                    
                     f=find(abs(abs(llP(k,1))-LonDeg(:,1))==min(abs(abs(llP(k,1))-LonDeg(:,1)))); % which lon are we in? ~ size of search window
                     lon2=100*LonDeg(f,2); %search in a greater rectangle of 100Lon*80Lat deg at the equator
                     if lon2>360 %the whole circunference (cases near the poles)
                        lon2=360;
                     end
                    lim=find(llT(:,2)<=llP(k,2)+lon2 & llT(:,2)>=llP(k,2)-lon2 & llT(:,1)<=llP(k,1)+40 & llT(:,1)>=llP(k,1)-40);
                    %%%%%%%%%%%%%%%%
                    if isempty(lim)==0
                        dist=NaN(length(lim),1);
                        posdel=zeros(length(lim),1);
                        posdel2=zeros(length(lim),1);
                        L=posdel;
                        for l=1:length(lim)
                            dist(l)=lldistkm(llP(k,:),llT(lim(l),:));
                            [x,y]=bresenham(Pr(k),Pc(k),Tr(lim(l)),Tc(lim(l))); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                            pos=sub2ind(size(X),x,y);
                            f2=find(lsoM(pos)==2);
                            ll=length(x);
                            if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                posdel(l)=1;
                            end
                            % Compute the line on the other direction
                            if Pc(k)<Tc(lim(l))
                                cP=Pc(k)+(size(Y,2)-Tc(lim(l)));
                                cT=1;
                                [x,y]=bresenham(Pr(k),cP,Tr(lim(l)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                y=y-(size(Y,2)-Tc(lim(l)));
                                y(y<=0)=size(Y,2)+y(y<=0);
                                pos=sub2ind(size(X),x,y);
                                f2=find(lsoM(pos)==2);
                                if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                    posdel2(l)=1;
                                end

                                if length(x)<ll
                                    L(l)=1; % transects shorter than the left to right direction
                                end
                            elseif Pc(k)>Tc(lim(l))
                                cT=Tc(lim(l))+(size(Y,2)-Pc(k));
                                cP=1;
                                [x,y]=bresenham(Pr(k),cP,Tr(lim(l)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                y=y-(size(Y,2)-Pc(k));
                                y(y<=0)=size(Y,2)+y(y<=0);
                                pos=sub2ind(size(X),x,y);
                                f2=find(lsoM(pos)==2);
                                if length(f2)/length(x)>=0.2 % if 20% of grids are land
                                    posdel2(l)=1;
                                end

                                if length(x)<ll
                                    L(l)=1; % transects shorter than the left to right direction
                                end
                            end
                        end

                        f1=find(L==0);
                        f=find(posdel(f1)==1);
                        f1=f1(f);

                        f2=find(L==1);
                        f=find(posdel2(f2)==1);
                        f2=f2(f);

                        L(unique([f1;f2]))=[];
                        lim(unique([f1;f2]))=[];
                        dist(unique([f1;f2]))=[];

                        if isempty(dist)==0
                            % shortest distance and find closest with no big land mass in between
                            f=find(dist==min(dist));
                            if length(f)>1
                                 f=f(1);
                            end
                            if L(f)==1
                                if Pc(k)<Tc(lim(f))
                                    cP=Pc(k)+(size(Y,2)-Tc(lim(f)));
                                    cT=1;
                                    [x,y]=bresenham(Pr(k),cP,Tr(lim(f)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                    y=y-(size(Y,2)-Tc(lim(f)));
                                    y(y<=0)=size(Y,2)+y(y<=0);
                                elseif Pc(k)>Tc(lim(f))
                                    cT=Tc(lim(f))+(size(Y,2)-Pc(k));
                                    cP=1;
                                    [x,y]=bresenham(Pr(k),cP,Tr(lim(f)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                    y=y-(size(Y,2)-Pc(k));
                                    y(y<=0)=size(Y,2)+y(y<=0);
                                end
                            elseif L(f)==0
                                [x,y]=bresenham(Pr(k),Pc(k),Tr(lim(f)),Tc(lim(f))); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                            end

                            Y(x,y);lat(1:length(x),k)=ans(:,1);
                            X(x,y);lon(1:length(x),k)=ans(1,:)';
                            d(x,y); div(1:length(x),k)=diag(ans)';
                        else
                            %weird peaks with no trough found in the big window without crossing more than 20% land (set manually)
                           lat(1,k)=llP(k,1);
                           lon(1,k)=llP(k,2);
                           div(1,k)=d(Pr(k),Pc(k));
                        end
                    end    
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
            else % if I do not find troughs in the area equivalent to 40Lon*20Lat at the equator rectangle, I increase the rectangle to 100x40
                %continue
                f=find(abs(abs(llP(k,1))-LonDeg(:,1))==min(abs(abs(llP(k,1))-LonDeg(:,1)))); % which lon are we in? ~ size of search window
                lon2=50*LonDeg(f,2); %search in a rectangle of 100Lon*40Lat deg at the equator
                if lon2>360 %the whole circunference (cases near the poles)
                    lon2=360;
                end
                lim=find(llT(:,2)<=llP(k,2)+lon2 & llT(:,2)>=llP(k,2)-lon2 & llT(:,1)<=llP(k,1)+20 & llT(:,1)>=llP(k,1)-20);
                %%%%%%%%%%%%%%%%%%%%%%%
                if isempty(lim)==0
                    dist=NaN(length(lim),1);
                    posdel=zeros(length(lim),1);
                    posdel2=zeros(length(lim),1);
                    L=posdel;
                    for l=1:length(lim)
                        dist(l)=lldistkm(llP(k,:),llT(lim(l),:));
                        [x,y]=bresenham(Pr(k),Pc(k),Tr(lim(l)),Tc(lim(l))); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                        pos=sub2ind(size(X),x,y);
                        f2=find(lsoM(pos)==2);
                        ll=length(x);
                        if length(f2)/length(x)>=0.2; % if 20% of grids are land
                            posdel(l)=1;
                        end
                        % Compute the line on the other direction
                        if Pc(k)<Tc(lim(l))
                            cP=Pc(k)+(size(Y,2)-Tc(lim(l)));
                            cT=1;
                            [x,y]=bresenham(Pr(k),cP,Tr(lim(l)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                            y=y-(size(Y,2)-Tc(lim(l)));
                            y(y<=0)=size(Y,2)+y(y<=0);
                            pos=sub2ind(size(X),x,y);
                            f2=find(lsoM(pos)==2);
                            if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                posdel2(l)=1;
                            end
                            
                            if length(x)<ll
                                L(l)=1; % transects shorter than the left to right direction
                            end
                        elseif Pc(k)>Tc(lim(l))
                            cT=Tc(lim(l))+(size(Y,2)-Pc(k));
                            cP=1;
                            [x,y]=bresenham(Pr(k),cP,Tr(lim(l)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                            y=y-(size(Y,2)-Pc(k));
                            y(y<=0)=size(Y,2)+y(y<=0);
                            pos=sub2ind(size(X),x,y);
                            f2=find(lsoM(pos)==2);
                            if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                posdel2(l)=1;
                            end
                            
                            if length(x)<ll
                                L(l)=1; % transects shorter than the left to right direction
                            end
                        end
                    end
                    
                    f1=find(L==0);
                    f=find(posdel(f1)==1);
                    f1=f1(f);
                    
                    f2=find(L==1);
                    f=find(posdel2(f2)==1);
                    f2=f2(f);
                    
                    L(unique([f1;f2]))=[];
                    lim(unique([f1;f2]))=[];
                    dist(unique([f1;f2]))=[];
                    
                    if isempty(dist)==0
                        % shortest distance and find closest with no big land mass in between
                        f=find(dist==min(dist));
                        if length(f)>1
                                 f=f(1);
                        end
                        if L(f)==1
                            if Pc(k)<Tc(lim(f))
                                cP=Pc(k)+(size(Y,2)-Tc(lim(f)));
                                cT=1;
                                [x,y]=bresenham(Pr(k),cP,Tr(lim(f)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                y=y-(size(Y,2)-Tc(lim(f)));
                                y(y<=0)=size(Y,2)+y(y<=0);
                            elseif Pc(k)>Tc(lim(f))
                                cT=Tc(lim(f))+(size(Y,2)-Pc(k));
                                cP=1;
                                [x,y]=bresenham(Pr(k),cP,Tr(lim(f)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                y=y-(size(Y,2)-Pc(k));
                                y(y<=0)=size(Y,2)+y(y<=0);
                            end
                        elseif L(f)==0
                            [x,y]=bresenham(Pr(k),Pc(k),Tr(lim(f)),Tc(lim(f))); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                        end
                        
                        Y(x,y);lat(1:length(x),k)=ans(:,1);
                        X(x,y);lon(1:length(x),k)=ans(1,:)';
                        d(x,y); div(1:length(x),k)=diag(ans)'; 
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    else
                        
                           lat(1,k)=llP(k,1);
                           lon(1,k)=llP(k,2);
                           div(1,k)=d(Pr(k),Pc(k));
                    end
                else
                    f=find(abs(abs(llP(k,1))-LonDeg(:,1))==min(abs(abs(llP(k,1))-LonDeg(:,1)))); % which lon are we in? ~ size of search window
                    lon2=30*LonDeg(f,2); %search in a rectangle of 60Lon*40Lat deg at the equator
                    if lon2>360 %the whole circunference (cases near the poles)
                        lon2=360;
                    end
                    lim=find(llT(:,2)<=llP(k,2)+lon2 & llT(:,2)>=llP(k,2)-lon2 & llT(:,1)<=llP(k,1)+20 & llT(:,1)>=llP(k,1)-20);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if isempty(lim)==0
                        dist=NaN(length(lim),1);
                        posdel=zeros(length(lim),1);
                        posdel2=zeros(length(lim),1);
                        L=posdel;
                        for l=1:length(lim)
                            dist(l)=lldistkm(llP(k,:),llT(lim(l),:));
                            [x,y]=bresenham(Pr(k),Pc(k),Tr(lim(l)),Tc(lim(l))); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                            pos=sub2ind(size(X),x,y);
                            f2=find(lsoM(pos)==2);
                            ll=length(x);
                            if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                posdel(l)=1;
                            end
                            % Compute the line on the other direction
                            if Pc(k)<Tc(lim(l))
                                cP=Pc(k)+(size(Y,2)-Tc(lim(l)));
                                cT=1;
                                [x,y]=bresenham(Pr(k),cP,Tr(lim(l)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                y=y-(size(Y,2)-Tc(lim(l)));
                                y(y<=0)=size(Y,2)+y(y<=0);
                                pos=sub2ind(size(X),x,y);
                                f2=find(lsoM(pos)==2);
                                if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                    posdel2(l)=1;
                                end
                                
                                if length(x)<ll
                                    L(l)=1; % transects shorter than the left to right direction
                                end
                            elseif Pc(k)>Tc(lim(l))
                                cT=Tc(lim(l))+(size(Y,2)-Pc(k));
                                cP=1;
                                [x,y]=bresenham(Pr(k),cP,Tr(lim(l)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                y=y-(size(Y,2)-Pc(k));
                                y(y<=0)=size(Y,2)+y(y<=0);
                                pos=sub2ind(size(X),x,y);
                                f2=find(lsoM(pos)==2);
                                if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                    posdel2(l)=1;
                                end
                                
                                if length(x)<ll
                                    L(l)=1; % transects shorter than the left to right direction
                                end
                            end
                        end
                        
                        f1=find(L==0);
                        f=find(posdel(f1)==1);
                        f1=f1(f);
                        
                        f2=find(L==1);
                        f=find(posdel2(f2)==1);
                        f2=f2(f);
                        
                        L(unique([f1;f2]))=[];
                        lim(unique([f1;f2]))=[];
                        dist(unique([f1;f2]))=[];
                        
                        if isempty(dist)==0
                            % shortest distance and find closest with no big land mass in between
                            f=find(dist==min(dist));
                            if length(f)>1
                                 f=f(1);
                            end
                            if L(f)==1
                                if Pc(k)<Tc(lim(f))
                                    cP=Pc(k)+(size(Y,2)-Tc(lim(f)));
                                    cT=1;
                                    [x,y]=bresenham(Pr(k),cP,Tr(lim(f)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                    y=y-(size(Y,2)-Tc(lim(f)));
                                    y(y<=0)=size(Y,2)+y(y<=0);
                                elseif Pc(k)>Tc(lim(f))
                                    cT=Tc(lim(f))+(size(Y,2)-Pc(k));
                                    cP=1;
                                    [x,y]=bresenham(Pr(k),cP,Tr(lim(f)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                    y=y-(size(Y,2)-Pc(k));
                                    y(y<=0)=size(Y,2)+y(y<=0);
                                end
                            elseif L(f)==0
                                [x,y]=bresenham(Pr(k),Pc(k),Tr(lim(f)),Tc(lim(f))); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                            end
                            
                            Y(x,y);lat(1:length(x),k)=ans(:,1);
                            X(x,y);lon(1:length(x),k)=ans(1,:)';
                            d(x,y); div(1:length(x),k)=diag(ans)';
                        else
                            
                           lat(1,k)=llP(k,1);
                           lon(1,k)=llP(k,2);
                           div(1,k)=d(Pr(k),Pc(k));                   
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    else
                        % break
                        f=find(abs(abs(llP(k,1))-LonDeg(:,1))==min(abs(abs(llP(k,1))-LonDeg(:,1)))); % which lon are we in? ~ size of search window
                        lon2=50*LonDeg(f,2); %search in a rectangle of 100Lon*40Lat deg at the equator
                        if lon2>360 %the whole circunference (cases near the poles)
                            lon2=360;
                        end
                        lim=find(llT(:,2)<=llP(k,2)+lon2 & llT(:,2)>=llP(k,2)-lon2 & llT(:,1)<=llP(k,1)+20 & llT(:,1)>=llP(k,1)-20);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%
                        if isempty(lim)==0
                            dist=NaN(length(lim),1);
                            posdel=zeros(length(lim),1);
                            posdel2=zeros(length(lim),1);
                            L=posdel;
                            for l=1:length(lim)
                                dist(l)=lldistkm(llP(k,:),llT(lim(l),:));
                                [x,y]=bresenham(Pr(k),Pc(k),Tr(lim(l)),Tc(lim(l))); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                pos=sub2ind(size(X),x,y);
                                f2=find(lsoM(pos)==2);
                                ll=length(x);
                                if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                    posdel(l)=1;
                                end
                                % Compute the line on the other direction
                                if Pc(k)<Tc(lim(l))
                                    cP=Pc(k)+(size(Y,2)-Tc(lim(l)));
                                    cT=1;
                                    [x,y]=bresenham(Pr(k),cP,Tr(lim(l)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                    y=y-(size(Y,2)-Tc(lim(l)));
                                    y(y<=0)=size(Y,2)+y(y<=0);
                                    pos=sub2ind(size(X),x,y);
                                    f2=find(lsoM(pos)==2);
                                    if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                        posdel2(l)=1;
                                    end
                                    
                                    if length(x)<ll
                                        L(l)=1; % transects shorter than the left to right direction
                                    end
                                elseif Pc(k)>Tc(lim(l))
                                    cT=Tc(lim(l))+(size(Y,2)-Pc(k));
                                    cP=1;
                                    [x,y]=bresenham(Pr(k),cP,Tr(lim(l)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                    y=y-(size(Y,2)-Pc(k));
                                    y(y<=0)=size(Y,2)+y(y<=0);
                                    pos=sub2ind(size(X),x,y);
                                    f2=find(lsoM(pos)==2);
                                    if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                        posdel2(l)=1;
                                    end
                                    
                                    if length(x)<ll
                                        L(l)=1; % transects shorter than the left to right direction
                                    end
                                end
                            end
                            
                            f1=find(L==0);
                            f=find(posdel(f1)==1);
                            f1=f1(f);
                            
                            f2=find(L==1);
                            f=find(posdel2(f2)==1);
                            f2=f2(f);
                            
                            L(unique([f1;f2]))=[];
                            lim(unique([f1;f2]))=[];
                            dist(unique([f1;f2]))=[];
                            
                            if isempty(dist)==0
                                % shortest distance and find closest with no big land mass in between
                                f=find(dist==min(dist));
                                if length(f)>1
                                 f=f(1);
                                end
                                if L(f)==1
                                    if Pc(k)<Tc(lim(f))
                                        cP=Pc(k)+(size(Y,2)-Tc(lim(f)));
                                        cT=1;
                                        [x,y]=bresenham(Pr(k),cP,Tr(lim(f)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                        y=y-(size(Y,2)-Tc(lim(f)));
                                        y(y<=0)=size(Y,2)+y(y<=0);
                                    elseif Pc(k)>Tc(lim(f))
                                        cT=Tc(lim(f))+(size(Y,2)-Pc(k));
                                        cP=1;
                                        [x,y]=bresenham(Pr(k),cP,Tr(lim(f)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                        y=y-(size(Y,2)-Pc(k));
                                        y(y<=0)=size(Y,2)+y(y<=0);
                                    end
                                elseif L(f)==0
                                    [x,y]=bresenham(Pr(k),Pc(k),Tr(lim(f)),Tc(lim(f))); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                end
                                
                                Y(x,y);lat(1:length(x),k)=ans(:,1);
                                X(x,y);lon(1:length(x),k)=ans(1,:)';
                                d(x,y); div(1:length(x),k)=diag(ans)';
                            else
                                
                                    f=find(abs(abs(llP(k,1))-LonDeg(:,1))==min(abs(abs(llP(k,1))-LonDeg(:,1)))); % which lon are we in? ~ size of search window
                                    lon2=100*LonDeg(f,2); %search in a rectangle of 200Lon*80Lat deg at the equator
                                    if lon2>360 %the whole circunference (cases near the poles)
                                        lon2=360;
                                    end  
                                    lim=find(llT(:,2)<=llP(k,2)+lon2 & llT(:,2)>=llP(k,2)-lon2 & llT(:,1)<=llP(k,1)+40 & llT(:,1)>=llP(k,1)-40);
                                    %%%%%%%%%%%%%%%%
                                    if isempty(lim)==0
                                        dist=NaN(length(lim),1);
                                        posdel=zeros(length(lim),1);
                                        posdel2=zeros(length(lim),1);
                                        L=posdel;
                                        for l=1:length(lim)
                                            dist(l)=lldistkm(llP(k,:),llT(lim(l),:));
                                            [x,y]=bresenham(Pr(k),Pc(k),Tr(lim(l)),Tc(lim(l))); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                            pos=sub2ind(size(X),x,y);
                                            f2=find(lsoM(pos)==2);
                                            ll=length(x);
                                            if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                                posdel(l)=1;
                                            end
                                            % Compute the line on the other direction
                                            if Pc(k)<Tc(lim(l))
                                                cP=Pc(k)+(size(Y,2)-Tc(lim(l)));
                                                cT=1;
                                                [x,y]=bresenham(Pr(k),cP,Tr(lim(l)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                                y=y-(size(Y,2)-Tc(lim(l)));
                                                y(y<=0)=size(Y,2)+y(y<=0);
                                                pos=sub2ind(size(X),x,y);
                                                f2=find(lsoM(pos)==2);
                                                if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                                    posdel2(l)=1;
                                                end
                                                
                                                if length(x)<ll
                                                    L(l)=1; % transects shorter than the left to right direction
                                                end
                                            elseif Pc(k)>Tc(lim(l))
                                                cT=Tc(lim(l))+(size(Y,2)-Pc(k));
                                                cP=1;
                                                [x,y]=bresenham(Pr(k),cP,Tr(lim(l)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                                y=y-(size(Y,2)-Pc(k));
                                                y(y<=0)=size(Y,2)+y(y<=0);
                                                pos=sub2ind(size(X),x,y);
                                                f2=find(lsoM(pos)==2);
                                                if length(f2)/length(x)>=0.2; % if 10% of grids are land
                                                    posdel2(l)=1;
                                                end
                                                
                                                if length(x)<ll
                                                    L(l)=1; % transects shorter than the left to right direction
                                                end
                                            end
                                        end
                                        
                                        f1=find(L==0);
                                        f=find(posdel(f1)==1);
                                        f1=f1(f);
                                        
                                        f2=find(L==1);
                                        f=find(posdel2(f2)==1);
                                        f2=f2(f);
                                        
                                        L(unique([f1;f2]))=[];
                                        lim(unique([f1;f2]))=[];
                                        dist(unique([f1;f2]))=[];
                                        
                                        if isempty(dist)==0
                                            % shortest distance and find closest with no big land mass in between
                                            f=find(dist==min(dist));
                                            if length(f)>1
                                              f=f(1);
                                            end
                                            if L(f)==1
                                                if Pc(k)<Tc(lim(f))
                                                    cP=Pc(k)+(size(Y,2)-Tc(lim(f)));
                                                    cT=1;
                                                    [x,y]=bresenham(Pr(k),cP,Tr(lim(f)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                                    y=y-(size(Y,2)-Tc(lim(f)));
                                                    y(y<=0)=size(Y,2)+y(y<=0);
                                                elseif Pc(k)>Tc(lim(f))
                                                    cT=Tc(lim(f))+(size(Y,2)-Pc(k));
                                                    cP=1;
                                                    [x,y]=bresenham(Pr(k),cP,Tr(lim(f)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                                    y=y-(size(Y,2)-Pc(k));
                                                    y(y<=0)=size(Y,2)+y(y<=0);
                                                end
                                            elseif L(f)==0
                                                [x,y]=bresenham(Pr(k),Pc(k),Tr(lim(f)),Tc(lim(f))); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                            end
                                            
                                            Y(x,y);lat(1:length(x),k)=ans(:,1);
                                            X(x,y);lon(1:length(x),k)=ans(1,:)';
                                            d(x,y); div(1:length(x),k)=diag(ans)';
                                        else
                                        end
                                    end    
                            end               
                        else
                            f=find(abs(abs(llP(k,1))-LonDeg(:,1))==min(abs(abs(llP(k,1))-LonDeg(:,1)))); % which lon are we in? ~ size of search window
                            lon2=50*LonDeg(f,2); %search in a rectangle of 100Lon*40Lat deg at the equator
                            if lon2>360 %the whole circunference (cases near the poles)
                                lon2=360;
                            end
                            lim=find(llT(:,2)<=llP(k,2)+lon2 & llT(:,2)>=llP(k,2)-lon2 & llT(:,1)<=llP(k,1)+20 & llT(:,1)>=llP(k,1)-20);
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%
                            if isempty(lim)==0
                                dist=NaN(length(lim),1);
                                posdel=zeros(length(lim),1);
                                posdel2=zeros(length(lim),1);
                                L=posdel;
                                for l=1:length(lim)
                                    dist(l)=lldistkm(llP(k,:),llT(lim(l),:));
                                    [x,y]=bresenham(Pr(k),Pc(k),Tr(lim(l)),Tc(lim(l))); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                    pos=sub2ind(size(X),x,y);
                                    f2=find(lsoM(pos)==2);
                                    ll=length(x);
                                    if length(f2)/length(x)>=0.2; % if 10% of grids are land
                                        posdel(l)=1;
                                    end
                                    % Compute the line on the other direction
                                    if Pc(k)<Tc(lim(l))
                                        cP=Pc(k)+(size(Y,2)-Tc(lim(l)));
                                        cT=1;
                                        [x,y]=bresenham(Pr(k),cP,Tr(lim(l)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                        y=y-(size(Y,2)-Tc(lim(l)));
                                        y(y<=0)=size(Y,2)+y(y<=0);
                                        pos=sub2ind(size(X),x,y);
                                        f2=find(lsoM(pos)==2);
                                        if length(f2)/length(x)>=0.2; % if 10% of grids are land
                                            posdel2(l)=1;
                                        end
                                        
                                        if length(x)<ll
                                            L(l)=1; % transects shorter than the left to right direction
                                        end
                                    elseif Pc(k)>Tc(lim(l))
                                        cT=Tc(lim(l))+(size(Y,2)-Pc(k));
                                        cP=1;
                                        [x,y]=bresenham(Pr(k),cP,Tr(lim(l)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                        y=y-(size(Y,2)-Pc(k));
                                        y(y<=0)=size(Y,2)+y(y<=0);
                                        pos=sub2ind(size(X),x,y);
                                        f2=find(lsoM(pos)==2);
                                        if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                            posdel2(l)=1;
                                        end
                                        
                                        if length(x)<ll
                                            L(l)=1; % transects shorter than the left to right direction
                                        end
                                    end
                                end
                                
                                f1=find(L==0);
                                f=find(posdel(f1)==1);
                                f1=f1(f);
                                
                                f2=find(L==1);
                                f=find(posdel2(f2)==1);
                                f2=f2(f);
                                
                                L(unique([f1;f2]))=[];
                                lim(unique([f1;f2]))=[];
                                dist(unique([f1;f2]))=[];
                                
                                if isempty(dist)==0
                                    % shortest distance and find closest with no big land mass in between
                                    f=find(dist==min(dist));
                                    if length(f)>1
                                       f=f(1);
                                    end
                                    if L(f)==1
                                        if Pc(k)<Tc(lim(f))
                                            cP=Pc(k)+(size(Y,2)-Tc(lim(f)));
                                            cT=1;
                                            [x,y]=bresenham(Pr(k),cP,Tr(lim(f)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                            y=y-(size(Y,2)-Tc(lim(f)));
                                            y(y<=0)=size(Y,2)+y(y<=0);
                                        elseif Pc(k)>Tc(lim(f))
                                            cT=Tc(lim(f))+(size(Y,2)-Pc(k));
                                            cP=1;
                                            [x,y]=bresenham(Pr(k),cP,Tr(lim(f)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                            y=y-(size(Y,2)-Pc(k));
                                            y(y<=0)=size(Y,2)+y(y<=0);
                                        end
                                    elseif L(f)==0
                                        [x,y]=bresenham(Pr(k),Pc(k),Tr(lim(f)),Tc(lim(f))); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                    end
                                    
                                    Y(x,y);lat(1:length(x),k)=ans(:,1);
                                    X(x,y);lon(1:length(x),k)=ans(1,:)';
                                    d(x,y); div(1:length(x),k)=diag(ans)';
                                else
                                          f=find(abs(abs(llP(k,1))-LonDeg(:,1))==min(abs(abs(llP(k,1))-LonDeg(:,1)))); % which lon are we in? ~ size of search window
                                          lon2=100*LonDeg(f,2);
                                          if lon2>360 %the whole circunference (cases near the poles)
                                            lon2=360;
                                          end
                                          lim=find(llT(:,2)<=llP(k,2)+lon2 & llT(:,2)>=llP(k,2)-lon2 & llT(:,1)<=llP(k,1)+40 & llT(:,1)>=llP(k,1)-40);
                                    %%%%%%%%%%%%%%%%
                                        if isempty(lim)==0
                                        dist=NaN(length(lim),1);
                                        posdel=zeros(length(lim),1);
                                        posdel2=zeros(length(lim),1);
                                        L=posdel;
                                            for l=1:length(lim)
                                            dist(l)=lldistkm(llP(k,:),llT(lim(l),:));
                                            [x,y]=bresenham(Pr(k),Pc(k),Tr(lim(l)),Tc(lim(l))); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                            pos=sub2ind(size(X),x,y);
                                            f2=find(lsoM(pos)==2);
                                            ll=length(x);
                                                if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                                posdel(l)=1;
                                                end
                                            % Compute the line on the other direction
                                                if Pc(k)<Tc(lim(l))
                                                cP=Pc(k)+(size(Y,2)-Tc(lim(l)));
                                                cT=1;
                                                [x,y]=bresenham(Pr(k),cP,Tr(lim(l)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                                y=y-(size(Y,2)-Tc(lim(l)));
                                                y(y<=0)=size(Y,2)+y(y<=0);
                                                pos=sub2ind(size(X),x,y);
                                                f2=find(lsoM(pos)==2);
                                                    if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                                    posdel2(l)=1;
                                                    end
                                                
                                                    if length(x)<ll
                                                    L(l)=1; % transects shorter than the left to right direction
                                                    end
                                                elseif Pc(k)>Tc(lim(l))
                                                cT=Tc(lim(l))+(size(Y,2)-Pc(k));
                                                cP=1;
                                                [x,y]=bresenham(Pr(k),cP,Tr(lim(l)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                                y=y-(size(Y,2)-Pc(k));
                                                y(y<=0)=size(Y,2)+y(y<=0);
                                                pos=sub2ind(size(X),x,y);
                                                f2=find(lsoM(pos)==2);
                                                    if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                                    posdel2(l)=1;
                                                    end
                                                
                                                    if length(x)<ll
                                                    L(l)=1; % transects shorter than the left to right direction
                                                    end
                                                end
                                            end
                                        
                                        f1=find(L==0);
                                        f=find(posdel(f1)==1);
                                        f1=f1(f);
                                        
                                        f2=find(L==1);
                                        f=find(posdel2(f2)==1);
                                        f2=f2(f);
                                        
                                        L(unique([f1;f2]))=[];
                                        lim(unique([f1;f2]))=[];
                                        dist(unique([f1;f2]))=[];
                                        
                                          if isempty(dist)==0
                                            % shortest distance and find closest with no big land mass in between
                                            f=find(dist==min(dist));
                                            if length(f)>1
                                               f=f(1);
                                            end
                                            if L(f)==1
                                                if Pc(k)<Tc(lim(f))
                                                    cP=Pc(k)+(size(Y,2)-Tc(lim(f)));
                                                    cT=1;
                                                    [x,y]=bresenham(Pr(k),cP,Tr(lim(f)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                                    y=y-(size(Y,2)-Tc(lim(f)));
                                                    y(y<=0)=size(Y,2)+y(y<=0);
                                                elseif Pc(k)>Tc(lim(f))
                                                    cT=Tc(lim(f))+(size(Y,2)-Pc(k));
                                                    cP=1;
                                                    [x,y]=bresenham(Pr(k),cP,Tr(lim(f)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                                    y=y-(size(Y,2)-Pc(k));
                                                    y(y<=0)=size(Y,2)+y(y<=0);
                                                end
                                            elseif L(f)==0
                                                [x,y]=bresenham(Pr(k),Pc(k),Tr(lim(f)),Tc(lim(f))); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                            end
                                            
                                            Y(x,y);lat(1:length(x),k)=ans(:,1);
                                            X(x,y);lon(1:length(x),k)=ans(1,:)';
                                            d(x,y); div(1:length(x),k)=diag(ans)';
                                          else
                                             
                                               lat(1,k)=llP(k,1);
                                               lon(1,k)=llP(k,2);
                                               div(1,k)=d(Pr(k),Pc(k));
                                         end
                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        end          
                                end
                            else
                                f=find(abs(abs(llP(k,1))-LonDeg(:,1))==min(abs(abs(llP(k,1))-LonDeg(:,1)))); % which lon are we in? ~ size of search window
                                lon2=60*LonDeg(f,2); %search in a rectangle of 120Lon*80Lat deg at the equator
                                if lon2>360 %the whole circunference (cases near the poles)
                                    lon2=360;
                                end
                                lim=find(llT(:,2)<=llP(k,2)+lon2 & llT(:,2)>=llP(k,2)-lon2 & llT(:,1)<=llP(k,1)+40 & llT(:,1)>=llP(k,1)-40);
                                
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                                if isempty(lim)==0
                                    dist=NaN(length(lim),1);
                                    posdel=zeros(length(lim),1);
                                    posdel2=zeros(length(lim),1);
                                    L=posdel;
                                    for l=1:length(lim)
                                        dist(l)=lldistkm(llP(k,:),llT(lim(l),:));
                                        [x,y]=bresenham(Pr(k),Pc(k),Tr(lim(l)),Tc(lim(l))); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                        pos=sub2ind(size(X),x,y);
                                        f2=find(lsoM(pos)==2);
                                        ll=length(x);
                                        if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                            posdel(l)=1;
                                        end
                                        % Compute the line on the other direction
                                        if Pc(k)<Tc(lim(l))
                                            cP=Pc(k)+(size(Y,2)-Tc(lim(l)));
                                            cT=1;
                                            [x,y]=bresenham(Pr(k),cP,Tr(lim(l)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                            y=y-(size(Y,2)-Tc(lim(l)));
                                            y(y<=0)=size(Y,2)+y(y<=0);
                                            pos=sub2ind(size(X),x,y);
                                            f2=find(lsoM(pos)==2);
                                            if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                                posdel2(l)=1;
                                            end
                                            
                                            if length(x)<ll
                                                L(l)=1; % transects shorter than the left to right direction
                                            end
                                        elseif Pc(k)>Tc(lim(l))
                                            cT=Tc(lim(l))+(size(Y,2)-Pc(k));
                                            cP=1;
                                            [x,y]=bresenham(Pr(k),cP,Tr(lim(l)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                            y=y-(size(Y,2)-Pc(k));
                                            y(y<=0)=size(Y,2)+y(y<=0);
                                            pos=sub2ind(size(X),x,y);
                                            f2=find(lsoM(pos)==2);
                                            if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                                posdel2(l)=1;
                                            end
                                            
                                            if length(x)<ll
                                                L(l)=1; % transects shorter than the left to right direction
                                            end
                                        end
                                    end
                                    
                                    f1=find(L==0);
                                    f=find(posdel(f1)==1);
                                    f1=f1(f);
                                    
                                    f2=find(L==1);
                                    f=find(posdel2(f2)==1);
                                    f2=f2(f);
                                    
                                    L(unique([f1;f2]))=[];
                                    lim(unique([f1;f2]))=[];
                                    dist(unique([f1;f2]))=[];
                                    
                                    if isempty(dist)==0
                                        % shortest distance and find closest with no big land mass in between
                                        f=find(dist==min(dist));
                                        if length(f)>1
                                           f=f(1);
                                        end
                                        if L(f)==1
                                            if Pc(k)<Tc(lim(f))
                                                cP=Pc(k)+(size(Y,2)-Tc(lim(f)));
                                                cT=1;
                                                [x,y]=bresenham(Pr(k),cP,Tr(lim(f)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                                y=y-(size(Y,2)-Tc(lim(f)));
                                                y(y<=0)=size(Y,2)+y(y<=0);
                                            elseif Pc(k)>Tc(lim(f))
                                                cT=Tc(lim(f))+(size(Y,2)-Pc(k));
                                                cP=1;
                                                [x,y]=bresenham(Pr(k),cP,Tr(lim(f)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                                y=y-(size(Y,2)-Pc(k));
                                                y(y<=0)=size(Y,2)+y(y<=0);
                                            end
                                        elseif L(f)==0
                                            [x,y]=bresenham(Pr(k),Pc(k),Tr(lim(f)),Tc(lim(f))); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                        end
                                        
                                        Y(x,y);lat(1:length(x),k)=ans(:,1);
                                        X(x,y);lon(1:length(x),k)=ans(1,:)';
                                        d(x,y); div(1:length(x),k)=diag(ans)';
                                    else                                     
                                          f=find(abs(abs(llP(k,1))-LonDeg(:,1))==min(abs(abs(llP(k,1))-LonDeg(:,1)))); % which lon are we in? ~ size of search window
                                          lon2=100*LonDeg(f,2);
                                          if lon2>360 %the whole circunference (cases near the poles)
                                            lon2=360;
                                          end
                                        lim=find(llT(:,2)<=llP(k,2)+lon2 & llT(:,2)>=llP(k,2)-lon2 & llT(:,1)<=llP(k,1)+40 & llT(:,1)>=llP(k,1)-40);
                                        %%%%%%%%%%%%%%%%
                                       if isempty(lim)==0
                                        dist=NaN(length(lim),1);
                                        posdel=zeros(length(lim),1);
                                        posdel2=zeros(length(lim),1);
                                        L=posdel;
                                         for l=1:length(lim)
                                            dist(l)=lldistkm(llP(k,:),llT(lim(l),:));
                                            [x,y]=bresenham(Pr(k),Pc(k),Tr(lim(l)),Tc(lim(l))); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                            pos=sub2ind(size(X),x,y);
                                            f2=find(lsoM(pos)==2);
                                            ll=length(x);
                                            if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                                posdel(l)=1;
                                            end
                                            % Compute the line on the other direction
                                            if Pc(k)<Tc(lim(l))
                                                cP=Pc(k)+(size(Y,2)-Tc(lim(l)));
                                                cT=1;
                                                [x,y]=bresenham(Pr(k),cP,Tr(lim(l)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                                y=y-(size(Y,2)-Tc(lim(l)));
                                                y(y<=0)=size(Y,2)+y(y<=0);
                                                pos=sub2ind(size(X),x,y);
                                                f2=find(lsoM(pos)==2);
                                                if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                                    posdel2(l)=1;
                                                end
                                                
                                                if length(x)<ll
                                                    L(l)=1; % transects shorter than the left to right direction
                                                end
                                            elseif Pc(k)>Tc(lim(l))
                                                cT=Tc(lim(l))+(size(Y,2)-Pc(k));
                                                cP=1;
                                                [x,y]=bresenham(Pr(k),cP,Tr(lim(l)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                                y=y-(size(Y,2)-Pc(k));
                                                y(y<=0)=size(Y,2)+y(y<=0);
                                                pos=sub2ind(size(X),x,y);
                                                f2=find(lsoM(pos)==2);
                                                if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                                    posdel2(l)=1;
                                                end
                                                
                                                if length(x)<ll
                                                    L(l)=1; % transects shorter than the left to right direction
                                                end
                                            end
                                          end
                                        
                                        f1=find(L==0);
                                        f=find(posdel(f1)==1);
                                        f1=f1(f);
                                        
                                        f2=find(L==1);
                                        f=find(posdel2(f2)==1);
                                        f2=f2(f);
                                        
                                        L(unique([f1;f2]))=[];
                                        lim(unique([f1;f2]))=[];
                                        dist(unique([f1;f2]))=[];
                                        
                                         if isempty(dist)==0
                                            % shortest distance and find closest with no big land mass in between
                                            f=find(dist==min(dist));
                                            if length(f)>1
                                               f=f(1);
                                            end
                                            if L(f)==1
                                                if Pc(k)<Tc(lim(f))
                                                    cP=Pc(k)+(size(Y,2)-Tc(lim(f)));
                                                    cT=1;
                                                    [x,y]=bresenham(Pr(k),cP,Tr(lim(f)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                                    y=y-(size(Y,2)-Tc(lim(f)));
                                                    y(y<=0)=size(Y,2)+y(y<=0);
                                                elseif Pc(k)>Tc(lim(f))
                                                    cT=Tc(lim(f))+(size(Y,2)-Pc(k));
                                                    cP=1;
                                                    [x,y]=bresenham(Pr(k),cP,Tr(lim(f)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                                    y=y-(size(Y,2)-Pc(k));
                                                    y(y<=0)=size(Y,2)+y(y<=0);
                                                end
                                            elseif L(f)==0
                                                [x,y]=bresenham(Pr(k),Pc(k),Tr(lim(f)),Tc(lim(f))); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                            end
                                            
                                            Y(x,y);lat(1:length(x),k)=ans(:,1);
                                            X(x,y);lon(1:length(x),k)=ans(1,:)';
                                            d(x,y); div(1:length(x),k)=diag(ans)';
                                         else
                                               lat(1,k)=llP(k,1);
                                               lon(1,k)=llP(k,2);
                                               div(1,k)=d(Pr(k),Pc(k)); 
                                          end
                                       else 
                                               lat(1,k)=llP(k,1);
                                               lon(1,k)=llP(k,2);
                                               div(1,k)=d(Pr(k),Pc(k));
                                       end
                                    end
                                else
                                    if lon>360 %the whole circunference (cases near the poles)
                                        lon=360;
                                    end
                                    lon2=100*LonDeg(f,2); %search in a rectangle of 200Lon*80Lat deg at the equator
                                    lim=find(llT(:,2)<=llP(k,2)+lon2 & llT(:,2)>=llP(k,2)-lon2 & llT(:,1)<=llP(k,1)+40 & llT(:,1)>=llP(k,1)-40);
                                    %%%%%%%%%%%%%%%%
                                    if isempty(lim)==0
                                        dist=NaN(length(lim),1);
                                        posdel=zeros(length(lim),1);
                                        posdel2=zeros(length(lim),1);
                                        L=posdel;
                                        for l=1:length(lim)
                                            dist(l)=lldistkm(llP(k,:),llT(lim(l),:));
                                            [x,y]=bresenham(Pr(k),Pc(k),Tr(lim(l)),Tc(lim(l))); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                            pos=sub2ind(size(X),x,y);
                                            f2=find(lsoM(pos)==2);
                                            ll=length(x);
                                            if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                                posdel(l)=1;
                                            end
                                            % Compute the line on the other direction
                                            if Pc(k)<Tc(lim(l))
                                                cP=Pc(k)+(size(Y,2)-Tc(lim(l)));
                                                cT=1;
                                                [x,y]=bresenham(Pr(k),cP,Tr(lim(l)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                                y=y-(size(Y,2)-Tc(lim(l)));
                                                y(y<=0)=size(Y,2)+y(y<=0);
                                                pos=sub2ind(size(X),x,y);
                                                f2=find(lsoM(pos)==2);
                                                if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                                    posdel2(l)=1;
                                                end
                                                
                                                if length(x)<ll
                                                    L(l)=1; % transects shorter than the left to right direction
                                                end
                                            elseif Pc(k)>Tc(lim(l))
                                                cT=Tc(lim(l))+(size(Y,2)-Pc(k));
                                                cP=1;
                                                [x,y]=bresenham(Pr(k),cP,Tr(lim(l)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                                y=y-(size(Y,2)-Pc(k));
                                                y(y<=0)=size(Y,2)+y(y<=0);
                                                pos=sub2ind(size(X),x,y);
                                                f2=find(lsoM(pos)==2);
                                                if length(f2)/length(x)>=0.2; % if 20% of grids are land
                                                    posdel2(l)=1;
                                                end
                                                
                                                if length(x)<ll
                                                    L(l)=1; % transects shorter than the left to right direction
                                                end
                                            end
                                        end
                                        
                                        f1=find(L==0);
                                        f=find(posdel(f1)==1);
                                        f1=f1(f);
                                        
                                        f2=find(L==1);
                                        f=find(posdel2(f2)==1);
                                        f2=f2(f);
                                        
                                        L(unique([f1;f2]))=[];
                                        lim(unique([f1;f2]))=[];
                                        dist(unique([f1;f2]))=[];
                                        
                                        if isempty(dist)==0
                                            % shortest distance and find closest with no big land mass in between
                                            f=find(dist==min(dist));
                                            if length(f)>1
                                               f=f(1);
                                            end
                                            if L(f)==1
                                                if Pc(k)<Tc(lim(f))
                                                    cP=Pc(k)+(size(Y,2)-Tc(lim(f)));
                                                    cT=1;
                                                    [x,y]=bresenham(Pr(k),cP,Tr(lim(f)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                                    y=y-(size(Y,2)-Tc(lim(f)));
                                                    y(y<=0)=size(Y,2)+y(y<=0);
                                                elseif Pc(k)>Tc(lim(f))
                                                    cT=Tc(lim(f))+(size(Y,2)-Pc(k));
                                                    cP=1;
                                                    [x,y]=bresenham(Pr(k),cP,Tr(lim(f)),cT); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                                    y=y-(size(Y,2)-Pc(k));
                                                    y(y<=0)=size(Y,2)+y(y<=0);
                                                end
                                            elseif L(f)==0
                                                [x,y]=bresenham(Pr(k),Pc(k),Tr(lim(f)),Tc(lim(f))); %Bresenhams line algorithm (I compute it on the grid coordinates to avoid rounding of LatLon coordinates (then I backtransform to obtain distances between grids along transects)
                                            end
                                            
                                            Y(x,y);lat(1:length(x),k)=ans(:,1);
                                            X(x,y);lon(1:length(x),k)=ans(1,:)';
                                            d(x,y); div(1:length(x),k)=diag(ans)';
                                        else
                                               lat(1,k)=llP(k,1);
                                               lon(1,k)=llP(k,2);
                                               div(1,k)=d(Pr(k),Pc(k)); 
                                        end
                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    else
                                               lat(1,k)=llP(k,1);
                                               lon(1,k)=llP(k,2);
                                               div(1,k)=d(Pr(k),Pc(k));
                                    end
                                end
                            end
                        end
                    end
                end 
            end
        end
          LAT{j}=lat;
          LON{j}=lon;
          DIV{j}=div; 
          clear dist Pr Pc Tr Tc llP llT lat lon div 
    end
    a=list(i).name;
     f=strfind(a,'_');
    a=a(1:f(1));
   
  eval(['save ' a 'transects D label_cases Point_timeslices X Y LAT LON DIV'])
   end
return

