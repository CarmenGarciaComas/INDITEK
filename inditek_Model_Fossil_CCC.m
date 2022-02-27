
function [CCC,Div1,Div2] = inditek_Model_Fossil_CCC(Q10,kfood,rhomin,rhomax,Kmax,Kmin,ext_pattern)
disp('** inditek_Model_Fossil_CCC.m **')

list=dir('*global.mat');
eval(['load ' list(1).name ' Point_timeslices']);


if ext_pattern==3
    ini=3;
else
    ini=1;
end

for i=1:length(list)
    [d,T]=xlsread('data/FossilTimeSeries');
    
    if ext_pattern==1
        t=d(:,1);
        d=d(:,2);
    elseif ext_pattern==2
        t=d(:,3);
        d=d(:,4);
    elseif ext_pattern==3
        t=d(:,5);
        d=d(:,6);
    end
    
    d=d(isnan(d)==0);t=t(isnan(d)==0);
    d=interp1(t,d,-1*Point_timeslices,'linear');
    
    eval(['load ' list(i).name]);
    D=gammaD;
    
    if ext_pattern==1
        NaNpos=[16:18,26:29,41:43,48:49,71];
    elseif ext_pattern==2
        NaNpos=[13:16,26:28,31:32,39:45,49:52,67:68];
    elseif ext_pattern==3
        NaNpos=[15:17,26:28,31:33,41:44,49:50,70];
    end
    
    D(ini)=1;
    Div=D;
    d(NaNpos)=NaN;D(NaNpos)=NaN;
    
    d=d(ini:end)';
    D=D(ini:end);
    % standardise t.s. to 0-1
    D=(D-nanmin(D))/(nanmax(D)-nanmin(D));
    d=(d-nanmin(d))/(nanmax(d)-nanmin(d));
    
    nptos = length(d);
    meanx = nanmean(d);
    meany = nanmean(D);
    sigmax = nanstd(d);
    sigmay = nanstd(D);
    rho_xy = corrcoef(d,D,'rows','pairwise');%corr(d,D);
    
    r= (2.0 * rho_xy * sigmax * sigmay) / (sigmax^2 + sigmay^2 + (meanx - meany)^2);
    CCC(i)=r(2);
    eval(['Div' num2str(i) '=Div;'])
end
return

