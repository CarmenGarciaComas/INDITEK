function [] = inditek_gridding_alphadiv(label_cases,ext_pattern)

disp('** inditek_gridding_alphadiv.m **')

% To calculate diversity by grid from particle diversity:
% I consider shelf particles for the shelf grids and ocean particles for the ocean grids independently to better keep shelf-ocean gradient

load data/landShelfOceanMask % 0=ocean/1=shelf/2=land
landShelfOceanMask=flip(landShelfOceanMask,3);
[X,Y]=meshgrid(landShelfOcean_Lon,landShelfOcean_Lat);

Xocean=cell(size(landShelfOceanMask,3),1);
Yocean=Xocean;
Xshelf=Xocean;
Yshelf=Xocean;
posocean=Xocean;
posshelf=Xocean;

if ext_pattern==3
    ini=3;
else
    ini=1;
end

for i=1:size(landShelfOceanMask,3)
  z=landShelfOceanMask(:,:,i)';    
  f=find(z==0);
  posocean{i}=f;
  Xocean{i}=X(f);Yocean{i}=Y(f);
  f=find(z==1);
  posshelf{i}=f;
  Xshelf{i}=X(f);Yshelf{i}=Y(f);
end

for i=1:length(label_cases)
    eval(['load INDITEK' label_cases{i}  '_alpha'])
    D=NaN(size(X,1),size(X,2),size(D_ocean,2));
    i
    for j=1:size(D_ocean,2)
      d=NaN(size(X));
      if j<=ini
        D(:,:,j)=1;
      else
        z=D_ocean(:,j);
        x=ocean_lonlatAge(:,j,1);
        y=ocean_lonlatAge(:,j,2);
	
        f=find(isnan(z)==0);
        z=z(f);x=x(f);y=y(f);
        f=find(isnan(x)==0); %for points initiated at t-1 to add diversity forward
        z=z(f);x=x(f);y=y(f);
        F = scatteredInterpolant(x,y,z);
        F.Method = 'natural';
        F.ExtrapolationMethod = 'linear';
        vq = F(Xocean{j},Yocean{j});
        d(posocean{j})=vq;

        z=D_shelf(:,j);
        x=shelf_lonlatAge(:,j,1);
        y=shelf_lonlatAge(:,j,2);
        f=find(isnan(z)==0);
        z=z(f);x=x(f);y=y(f);
        f=find(isnan(x)==0);
        z=z(f);x=x(f);y=y(f);
        F = scatteredInterpolant(x,y,z);
        F.Method = 'natural';
        F.ExtrapolationMethod = 'linear';
        vq = F(Xshelf{j},Yshelf{j});
        d(posshelf{j})=vq;
        D(:,:,j)=d;
      end     
      % Avoid negative values due to extrapolation from very low Ds
      f=find(D<0);
      D(f)=0;
      
    end
    eval(['save INDITEK' label_cases{i} '_grid'  ' D X Y Point_timeslices']);
end

return

