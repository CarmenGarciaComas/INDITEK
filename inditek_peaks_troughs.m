function [] = inditek_peaks_troughs(label_cases,ext_pattern)

disp('** inditek_peaks_troughs.m **')

% output location of ocean ridges (troughs) and local maxima from extrema2 with D>=75% of D in the local maxima (peaks) 

load data/POSoceanRidge
PosMIN=POSoceanRidge; % Position of troughs is the location of new ocean crust
if ext_pattern==3
    ini=3;
else
    ini=1;
end
tf=82;
PosMAX=cell(tf,2);
for i=1:length(label_cases)
    i
    % Calculate and save location and values of peaks by map
    eval(['load INDITEK' label_cases{i} '_grid']);
    PosMAX=cell(82,1);
    for j=ini+1:size(D,3) %jump 1st time frames that have a single value Do=1 (depending on ext. pattern it will be different);
        d=D(:,:,j);
        [Dmax,Posmax,Dmin,Posmin] = extrema2(d); % local maxima and minima and their position in the grids
        f1=Dmax>=quantile(Dmax,0.75);
        PosMAX{j}=Posmax(f1);% Position of peaks is the location of local maxima with D>75% of maxima in the time frame i
    end
   
    eval(['save INDITEK' label_cases{i} '_PT PosM* D X Y label_cases Point_timeslices']);
end

return


