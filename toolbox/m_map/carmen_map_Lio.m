clear all
close all

figure(1)
color=colormaplus('jet',100,'log+');
colormap(flipud(color))
m_proj('albers equal-area','lat',[24 34],'long',[117 130]);
m_gshhs_f('patch',[.6 .6 .6]);
m_elev('contourf',[-10:-100:-8000],'EdgeColor','none') 
set(gca,'color','none')
m_grid('box','fancy')
xlabel('Longitude')
ylabel('Latitude')
colorbar
