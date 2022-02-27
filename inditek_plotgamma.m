function [] = inditek_plotgamma(ext_pattern)

disp('** inditek_plotgamma.m **')
  
% Relationship of model diversity time series with fossil diversity time series

list=dir('*global.mat');

eval(['load ' list(1).name ' Point_timeslices']);
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
if ext_pattern==1
    NaNpos=[16:18,26:29,41:43,48:49,71];
elseif ext_pattern==2
    NaNpos=[13:16,26:28,31:32,39:45,49:52,67:68];
elseif ext_pattern==3
    NaNpos=[15:17,26:28,31:33,41:44,49:50,70];
end

d=d(isnan(d)==0);t=t(isnan(d)==0);
d=interp1(t,d,-1*Point_timeslices,'linear');
d(NaNpos)=NaN;
dnorm=(d-nanmin(d))./(nanmax(d)-nanmin(d));

if ext_pattern==3
    ini=3;
else
    ini=1;
end
    
label1_model={'Exponential model','Logistic model'};
label2_model={'Exponential model normalised diversity','Logistic model normalised diversity'};
n=0;

figure(1)
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [12 7]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 12 7]);

for i=1:length(list)

   if i == 1 %exponential
    
     eval(['load ' list(i).name]); 
     D=gammaD;
     D(1:ini)=1; 
     D(NaNpos)=NaN;
     M=nanmax(D)+0.1*nanmax(D); 
     Dnorm=(D-nanmin(D))./(nanmax(D)-nanmin(D));
     
    nptos = length(dnorm(ini:end));
    meanx = nanmean(dnorm(ini:end));
    meany = nanmean(Dnorm(ini:end));
    sigmax = nanstd(dnorm(ini:end));
    sigmay = nanstd(Dnorm(ini:end));
    rho_xy = corrcoef(dnorm(ini:end)',Dnorm(ini:end),'rows','pairwise');
    
    ccc= (2.0 * rho_xy(2) * sigmax * sigmay) / (sigmax^2 + sigmay^2 + (meanx - meany)^2);
    
      n=n+1;  
  
     subplot(2,3,1)
     
     plot(-1.*Point_timeslices,D,'k.-','markerfacecolor','k')			
     hold on
     
     if ext_pattern==1
     x_points = [-451,-451,-441,-441];  
     elseif ext_pattern==2
     x_points = [-471,-471,-452,-452];
     elseif ext_pattern==3
     x_points = [-456,-456,-442,-442];
     end
     y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);

     hold on
     if ext_pattern==1
     x_points = [-391,-391,-367,-367];
     elseif ext_pattern==2
     x_points = [-399,-399,-378,-378];
     elseif ext_pattern==3
     x_points = [-396,-396,-374,-374];
     end
     y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);
     
     hold on
     if ext_pattern==1
     x_points = [-269,-269,-251,-251];
     elseif ext_pattern==2
     x_points = [-275,-275,-244,-244];
     elseif ext_pattern==3
     x_points = [-343,-343,-317,-317];
     end
     y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);
     
     hold on
     if ext_pattern==1
     x_points = [-225,-225,-202,-202];
     elseif ext_pattern==2
     x_points = [-207,-207,-185,-185];   
     elseif ext_pattern==3
     x_points = [-269,-269,-250,-250];  
     end
     y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);
     hold on
     if ext_pattern==3
        x_points = [-215,-215,-199,-199];
     end
      y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);
     hold on
     if ext_pattern==1
     x_points = [-64,-64,-62,-62]; 
     elseif ext_pattern==2
     x_points = [-86,-86,-75,-75];
      elseif ext_pattern==3
     x_points = [-69,-69,-65,-65];
     end
     y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);
   
     a=ylabel('# Genera');
     set(a,'FontName','times new roman','fontsize',12);
     a=xlabel('Time (MA)');
     set(a,'FontName','times new roman','fontsize',12);
     a=title(label1_model{1});
     set(a,'FontName','times new roman','fontsize',12);
     set(gca,'xlim',[-541,0])
     set(gca,'ylim',[0,M],'ytick',[0:1000:M],'xtick',[-500:100:0],'FontName','times new roman','fontsize',10)

     subplot(2,3,2)
     
     plot(-1.*Point_timeslices,Dnorm,'-','color',[0,0,0])
     hold on
     plot(-1.*Point_timeslices,dnorm,'-','color',[0.6,0.6,0.6])
     hold on
     if ext_pattern==1
     x_points = [-451,-451,-441,-441];  
     elseif ext_pattern==2
     x_points = [-471,-471,-452,-452];
     elseif ext_pattern==3
     x_points = [-456,-456,-442,-442]; 
     end
     y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);

     hold on
     if ext_pattern==1
     x_points = [-391,-391,-367,-367];
     elseif ext_pattern==2
     x_points = [-399,-399,-378,-378];
     elseif ext_pattern==3
     x_points = [-396,-396,-374,-374];
     end
     y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);
     
     hold on
     if ext_pattern==1
     x_points = [-269,-269,-251,-251];
     elseif ext_pattern==2
     x_points = [-275,-275,-244,-244];
     elseif ext_pattern==3
     x_points = [-343,-343,-317,-317];
     end
     y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);
     
     hold on
     if ext_pattern==1
     x_points = [-225,-225,-202,-202];
     elseif ext_pattern==2
     x_points = [-207,-207,-185,-185];   
     elseif ext_pattern==3
     x_points = [-269,-269,-250,-250]; 
     end
     y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);
     hold on
     if ext_pattern==3
        x_points = [-215,-215,-199,-199];
     end
      y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);
     hold on
     if ext_pattern==1
     x_points = [-64,-64,-62,-62]; 
     elseif ext_pattern==2
     x_points = [-86,-86,-75,-75];
      elseif ext_pattern==3
     x_points = [-69,-69,-65,-65];
     end
     y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);
     
     a=ylabel('Normalised diversity');
     set(a,'FontName','times new roman','fontsize',12);
     a=xlabel('Time(MA)');
     set(a,'FontName','times new roman','fontsize',12);
     a=title(label1_model{1});
     set(a,'FontName','times new roman','fontsize',12);
     a=legend({'model','fossil'},'Location','northwest');
     legend boxoff
     set(a,'FontName','times new roman','fontsize',10);
     set(gca,'ylim',[0,1.1],'xlim',[-541,0],'xtick',[-500:100:0],'FontName','times new roman','fontsize',10)
     
     subplot(2,3,3)
     
     plot(dnorm(ini+1:end),Dnorm(ini+1:end),'ko','markerfacecolor','k')
     hold on
     p2 = polyfit(dnorm(ini+1:end)',Dnorm(ini+1:end),1);
     y=polyval(p2,dnorm(ini+1:end)');
     plot(dnorm(ini+1:end),y,'k-')
     a=text(0.02,0.9,['r^2= ' num2str(rho_xy(2)^2,2)]);
     set(a,'FontName','times new roman','fontsize',12,'color','k');
     a=text(0.02,0.83,['ccc= ' num2str(ccc,2)]);
     set(a,'FontName','times new roman','fontsize',12);
     a=xlabel('Fossil normalised diversity');
     set(a,'FontName','times new roman','fontsize',12);
     a=ylabel(label2_model{n});
     set(a,'FontName','times new roman','fontsize',12);
     set(gca,'ylim',[-0.01,1.01],'xlim',[-0.01,1.01],'ytick',[0,0.5,1],'xtick',[0,0.5,1],'FontName','times new roman','fontsize',10)
     
  elseif i == 2 %logistic
    
     eval(['load ' list(i).name]);
     D=gammaD;
     D(1:ini)=1;
     
     Dnorm=(D-min(D))./(max(D)-min(D));
     
     
    nptos = length(dnorm(ini:end));
    meanx = nanmean(dnorm(ini:end));
    meany = nanmean(Dnorm(ini:end));
    sigmax = nanstd(dnorm(ini:end));
    sigmay = nanstd(Dnorm(ini:end));
    rho_xy = corrcoef(dnorm(ini:end)',Dnorm(ini:end),'rows','pairwise');
    
    ccc= (2.0 * rho_xy(2) * sigmax * sigmay) / (sigmax^2 + sigmay^2 + (meanx - meany)^2);
    
     n=n+1;  
     
     subplot(2,3,4)
     plot(-1.*Point_timeslices,D,'k.-','markerfacecolor','k')
     
     hold on
     
     if ext_pattern==1
     x_points = [-451,-451,-441,-441];  
     elseif ext_pattern==2
     x_points = [-471,-471,-452,-452];
     elseif ext_pattern==3
     x_points = [-456,-456,-442,-442];
     end
     y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);

     hold on
     if ext_pattern==1
     x_points = [-391,-391,-367,-367];
     elseif ext_pattern==2
     x_points = [-399,-399,-378,-378];
     elseif ext_pattern==3
     x_points = [-396,-396,-374,-374];
     end
     y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);
     
     hold on
     if ext_pattern==1
     x_points = [-269,-269,-251,-251];
     elseif ext_pattern==2
     x_points = [-275,-275,-244,-244];
     elseif ext_pattern==3
     x_points = [-343,-343,-317,-317];
     end
     y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);
     
     hold on
     if ext_pattern==1
     x_points = [-225,-225,-202,-202];
     elseif ext_pattern==2
     x_points = [-207,-207,-185,-185];   
     elseif ext_pattern==3
     x_points = [-269,-269,-250,-250];  
     end
     y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);
     hold on
     if ext_pattern==3
        x_points = [-215,-215,-199,-199];
     end
      y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);
     hold on
     if ext_pattern==1
     x_points = [-64,-64,-62,-62]; 
     elseif ext_pattern==2
     x_points = [-86,-86,-75,-75];
      elseif ext_pattern==3
     x_points = [-69,-69,-65,-65];
     end
     y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);
     
     a=ylabel('# Genera');
     set(a,'FontName','times new roman','fontsize',12);
     a=title(label1_model{2});
     set(a,'FontName','times new roman','fontsize',12);
     a=xlabel('Time (MA)');
     
     set(a,'FontName','times new roman','fontsize',12);
     set(gca,'ylim',[0,M],'xlim',[-541,0],'ytick',[0:1000:M],'xtick',[-500:100:0],'FontName','times new roman','fontsize',10)


     subplot(2,3,5)
     plot(-1.*Point_timeslices,Dnorm,'k-','color',[0,0,0])
     hold on
     plot(-1.*Point_timeslices,dnorm,'-','color',[0.6,0.6,0.6])
     
     hold on
     
     if ext_pattern==1
     x_points = [-451,-451,-441,-441];  
     elseif ext_pattern==2
     x_points = [-471,-471,-452,-452];
     elseif ext_pattern==3
     x_points = [-456,-456,-442,-442]; 
     end
     y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);

     hold on
     if ext_pattern==1
     x_points = [-391,-391,-367,-367];
     elseif ext_pattern==2
     x_points = [-399,-399,-378,-378];
     elseif ext_pattern==3
     x_points = [-396,-396,-374,-374];
     end
     y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);
     
     hold on
     if ext_pattern==1
     x_points = [-269,-269,-251,-251];
     elseif ext_pattern==2
     x_points = [-275,-275,-244,-244];
     elseif ext_pattern==3
     x_points = [-343,-343,-317,-317];
     end
     y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);
     
     hold on
     if ext_pattern==1
     x_points = [-225,-225,-202,-202];
     elseif ext_pattern==2
     x_points = [-207,-207,-185,-185];   
     elseif ext_pattern==3
     x_points = [-269,-269,-250,-250];  
     end
     y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);
     hold on
     if ext_pattern==3
        x_points = [-215,-215,-199,-199];
     end
      y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);
     hold on
     if ext_pattern==1
     x_points = [-64,-64,-62,-62]; 
     elseif ext_pattern==2
     x_points = [-86,-86,-75,-75];
      elseif ext_pattern==3
     x_points = [-69,-69,-65,-65];
     end
     y_points = [0, M, M, 0];
     color = [0.7,0.7,0.7];
     hold on;
     a = fill(x_points, y_points, color,'LineStyle','none');
     alpha(a,.1);
     
     a=ylabel('Normalised diversity');
     set(a,'FontName','times new roman','fontsize',12);
     a=xlabel('Time(MA)');
     set(a,'FontName','times new roman','fontsize',12);
     a=title(label1_model{n});
     set(a,'FontName','times new roman','fontsize',12);
     a=legend({'model','fossil'},'Location','northeast');
     legend boxoff
     set(a,'FontName','times new roman','fontsize',10);
     set(gca,'ylim',[0,1.1],'xlim',[-541,0],'xtick',[-500:100:0],'FontName','times new roman','fontsize',10)
     
     subplot(2,3,6)
     plot(dnorm(ini+1:end),Dnorm(ini+1:end),'ko','markerfacecolor','k')
     hold on
     p2 = polyfit(dnorm(ini+1:end)',Dnorm(ini+1:end),1);
     y=polyval(p2,dnorm(ini+1:end)');
     plot(dnorm(ini+1:end),y,'k-')
     a=text(0.02,0.9,['r^2= ' num2str(rho_xy(2)^2,2)]);
     set(a,'FontName','times new roman','fontsize',12,'color','k');
     a=text(0.02,0.83,['ccc= ' num2str(ccc,2)]);
     set(a,'FontName','times new roman','fontsize',12);
     a=xlabel('Fossil normalised diversity');
     set(a,'FontName','times new roman','fontsize',12);
     a=ylabel(label2_model{i});
     set(a,'FontName','times new roman','fontsize',12);
     set(gca,'ylim',[-0.01,1.01],'xlim',[-0.01,1.01],'ytick',[0,0.5,1],'xtick',[0,0.5,1],'FontName','times new roman','fontsize',10)
  end 
end 

if ext_pattern==1
    name='INDITEK_ModelvsFossil_Zaffos';
elseif ext_pattern==2
    name='INDITEK_ModelvsFossil_Alroy';
elseif ext_pattern==3
    name='INDITEK_ModelvsFossil_Sepkoski';
end

print ('-dpng', '-r200', name)
return

