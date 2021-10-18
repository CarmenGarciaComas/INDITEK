close all
clear all
% folder with public matlab tools used in our code:
addpath(genpath('~/toolbox/'))
%
disp('CHECK (1) that in you have downloaded m_maps to your toolbox folder, and (2) that you are in the main folder with files: INDITEKexponential_grid.mat & INDITEKlogistic_grid.mat for figure 1 / INDITEKexponential_global.mat & INDITEKlogistic_global.mat for figure 2 / inditek_LogisticModel_KminKmax.mat for figure 3. If it is OK press any key to continue')
pause
%################################ FIGURE 1: DIVERSITY MAPS ################################
m_proj('Mollweide','clongitude',0);
timeframes=[25,36,69,82];
load INDITEKexponential_grid.mat
TF=[400,300,70,0];
Dexp=D(:,:,timeframes);
M=log10(256);
load INDITEKlogistic_grid.mat
Dlog=D(:,:,timeframes);
close all
figure
set(gcf, 'PaperUnits', 'centimeters','papersize',[15 20])
set(gcf,'paperposition',[0 0 15 20])
POS=[0.01 0.74 0.45 0.18;0.52 0.74 0.45 0.18;0.01 0.53 0.45 0.18;0.52 0.53 0.45 0.18;0.01 0.32 0.45 0.18;0.52 0.32 0.45 0.18;0.008 0.11 0.45 0.18;0.52 0.11 0.45 0.18];
Let={'a','e','b','f','c','g','d','h'};
sp=0;
for i=1:length(timeframes)
    sp=sp+1;
    subplot(4,2,sp)
    set(gca,'position',POS(sp,:))
    m_pcolor(X,Y,log10(Dlog(:,:,i)));
    m_grid('xaxis','middle','xticklabel',[],'yticklabel',[],'backcolor',[0.9 0.9 0.9],'LineWidth',0.6,'linestyle',':');
    a=title([num2str(TF(i)),' MA']);
    set(a,'FontName','arial','fontsize',12)
    a=text(-2,1,Let{sp});
    set(a,'FontName','arial','fontsize',12)
    caxis([log10(1),M])
    colormap(jet(50));
    if sp==1
        a=text(-0.96,1.8,'Logistic model');
        set(a,'FontName','arial','fontsize',12,'fontweight','bold')
    end
    sp=sp+1;
    subplot(4,2,sp)
    set(gca,'position',POS(sp,:))
    m_pcolor(X,Y,log10(Dexp(:,:,i)));
    m_grid('xaxis','middle','xticklabel',[],'yticklabel',[],'backcolor',[0.9 0.9 0.9],'LineWidth',0.1,'linestyle',':');
    a=title([num2str(TF(i)),' MA']);
    set(a,'FontName','arial','fontsize',12)
    a=text(-2,1,Let{sp});
    set(a,'FontName','arial','fontsize',12)
    caxis([log10(1),M])
    colormap(jet(50));
    if sp==2
        a=text(-1.18,1.8,'Exponential model');
        set(a,'FontName','arial','fontsize',12,'fontweight','bold')
    end
    if i==length(timeframes)
        h = colorbar('location','SouthOutside','Position',[0.55,5,0.35,0.03]);
        set(h,'xtick',[log10(1),log10(2),log10(4),log10(8),log10(16),log10(32),log10(64),log10(128),log10(256)],'xticklabel',{'1','2','4','8','16','32','64','128','256'})
        set(h,'position',[0.32,0.08,0.3,0.016],'FontName','arial','fontsize',9)
        text(-4,-2.1,'Diversity (# genera area^-^1)','FontName','arial','fontsize',12);
    end
    
end
print ('Inditek_Figure1','-djpeg', '-r200')
close
% ################################### FIGURE 2: TIME SERIES ##################

clear all;close all;clc
figure
set(gcf, 'PaperUnits', 'centimeters','papersize',[20 10])
set(gcf,'paperposition',[0 0 18 8])
load INDITEKexponential_global
Dexp=gammaD;
load INDITEKlogistic_global
Dlog=gammaD;
[d,T]=xlsread('data/FossilTimeSeries');
if ext_pattern==3
    ini=4;
else
    ini=2;
end
if ext_pattern==1
    t=d(:,1);
    d1=d(:,2);
elseif ext_pattern==2
    t=d(:,3);
    d1=d(:,4);
elseif ext_pattern==3
    t=d(:,5);
    d1=d(:,6);
end
d1=d1(isnan(d1)==0);t=t(isnan(d1)==0);
d1=interp1(t,d1,-1*Point_timeslices,'linear');
[AX H1 H2]=plotyy(-1*Point_timeslices, [Dexp,Dlog], -1*Point_timeslices, d1, 'plot');
set(AX(1),'linewidth',0.75,'FontName','arial','fontsize',12);
set(AX(2),'linewidth',0.75,'FontName','arial','fontsize',12);
set(H1(1),'LineStyle','-','Color','b','LineWidth',1,'Marker','o','MarkerSize',3,'MarkerEdgeColor','b','Markerfacecolor','b')
set(H1(2),'LineStyle','-','Color','r','LineWidth',1,'Marker','o','MarkerSize',3,'MarkerEdgeColor','r','Markerfacecolor','r')
set(H2(1),'LineStyle',':','Color',[0.45,0.45,0.45],'LineWidth',1.1)
set(AX,{'ycolor'},{'k';'k'})
set(AX,'xlim',[-541,0],'xtick',[])
set(AX(1),'xlim',[-541,0],'xtick',[-500:100:0])
M=5000;
if ext_pattern==1
    M2=2500;m2=500;
elseif ext_pattern==2
    M2=1000;m2=200;
elseif ext_pattern==3
    M2=5000;m2=1000;
end
set(get(AX(1),'Ylabel'),'String','Model diversity (# genera)','FontName','arial','fontsize',12);
set(AX(1),'ylim',[0,M],'ytick',[0:1000:M])
set(AX(2),'ylim',[0,M2],'ytick',[0:m2:M2])
set(get(AX(2),'Ylabel'),'String','Fossil diversity (# genera)','rotation',-90, 'VerticalAlignment','middle','FontName','arial','fontsize',12);
xh = get(AX(2),'Ylabel');
p = get(xh,'position');
p(1) = 1.15*p(1);
set(xh,'position',p);
a=xlabel('Time (MA)');
set(a,'FontName','arial','fontsize',12);
hold on
a=legend('Exponential diversification','Logistic diversification','Fossil record');
legend boxoff
set(a,'Location','Northwest','FontName','arial','fontsize',12);
if ext_pattern==1
    x_points = [-451,-451,-441,-441];
elseif ext_pattern==2
    x_points = [-471,-471,-452,-452];
elseif ext_pattern==3
    x_points = [-448,-448,-437,-437];
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
    x_points = [-393,-393,-353,-353];
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
    x_points = [-253,-253,-245,-245];
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
    x_points = [-220,-220,-200,-200];
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
    x_points = [-64,-64,-59,-59];
end
y_points = [0, M, M, 0];
color = [0.7,0.7,0.7];
hold on;
a = fill(x_points, y_points, color,'LineStyle','none');
alpha(a,.1);
if ext_pattern==1
    NaNpos=[16:18,26:29,41:43,48:49,71];
elseif ext_pattern==2
    NaNpos=[13:16,26:28,31:32,39:45,49:52,67:68];
elseif ext_pattern==3
    NaNpos=[17:19,26:30,44:45,48:50,60,71];
end
print ('Inditek_Figure2','-djpeg', '-r200')


%######################### FIGURE 3: KminKmax CCC ##################################

load 'inditek_LogisticModel_KminKmax'
figure
set(gcf, 'PaperUnits', 'centimeters','papersize',[10 8])
set(gcf,'paperposition',[0 0 9 7])
CCC(CCC==0)=NaN;
imAlpha=ones(size(CCC));
imAlpha(isnan(CCC))=0;
imagesc(CCC,'AlphaData',imAlpha);
imagesc(1:size(CCC,1),1:size(CCC,2),CCC,'AlphaData',imAlpha)
colormap(jet(50));
set(gca,'Fontname','arial','fontsize',9)
a=xlabel('Kmin (# genera area^-^1)');
set(a,'Fontname','arial','fontsize',9);
a=ylabel('Kmax (# genera area^-^1)');
set(a,'Fontname','arial','fontsize',9);
set(gca, 'ytick',[1:size(CCC,1)],'yticklabel',{'256','128','64','32','16','8','4','2'},'xtick',[1:size(R,1)],'xticklabel',{'2','4','8','16','32','64','128','256'},'linewidth',1.2);
caxis([0,1])
colorbar
hold on
print ('Inditek_Figure3','-djpeg', '-r200')

clear all;close all;clc
