
%% Load data: 

filename = 'AERODEM_Korsgaard_mosaic.nc'; 

x = ncread(filename,'x'); 
y = ncread(filename,'y'); 
sfz = permute(ncread(filename,'surface'),[2 1]);
year = permute(ncread(filename,'year_mask'),[2 1]);
reliab = permute(ncread(filename,'reliability_mask'),[2 1]);

%% Plot data: 
% subsubplot, cmocean, and imagescn are from Climate Data Tools for Matlab. 

figure
subsubplot(1,3,1)
imagescn(x,y,sfz)
axis image off
bedmachine('gl','greenland','color',.5*[1 1 1],'linewidth',0.3)
title('surface','fontname','courier','interpreter','none')

subsubplot(1,3,2)
hr=imagescn(x,y,reliab);
hr.AlphaData = reliab>0; 
axis image off
cmocean algae
bedmachine('gl','greenland','color',.5*[1 1 1],'linewidth',0.3)
title('reliability_mask','fontname','courier','interpreter','none')

subsubplot(1,3,3)
h=imagescn(x,y,year); 
h.AlphaData = year>0; 
axis image off
bedmachine('gl','greenland','color',.5*[1 1 1],'linewidth',0.3)
caxis([79 87])
title('year_mask','fontname','courier','interpreter','none')

%exportgraphics(gcf,'aerodem_mosaic.jpg','resolution',600)
