%% Plot Fig. 02 - the study area with calibration grid cells
% by Tina Trautmann, Aug. 2022

%% set save pth
spth = [pwd '/figures/'];
if ~exist(spth, 'dir'), mkdir(spth),end


%% load data
pthIn = [pwd '/data/input/studyArea/globalBaseline_Constraints_1deg.mat']

load(pthIn, 'lat', 'lon', 'time')


%  the lat & lon of opti
opti = load([pwd '/data/input/opti/lat_lon_904.mat'], 'lat', 'lon');

%% Prep Plotting Stuff
% space stuff
pix_a = AreaGridLatLon(lat,lon,[1 1]);
pix_a = pix_a(:,1);

nPix   = size(pix_a,1);


% get the zones
load('data/input/ancillary/clusterRegions.mat', 'CLregions');
[pix_x,pix_y] = LatLon2PixelCoord(lon,lat,90,-180,1);
idx = sub2ind([180, 360],pix_y, pix_x);

KG_map   = CLregions;

KG_v        = KG_map(idx); %valids in map
KG_v(KG_v==0) = 4;
KG_v(KG_v==4) = 10;
KG_v(KG_v==2) = 20;
KG_v(KG_v==3) = 30;
KG_v(KG_v==1) = 40;
KG_v(KG_v==5) = 50;
KG_v = KG_v ./ 10;

zoneNames   = {'R1- Cold', 'R2- Temperate', 'R3- Humid',  'R4- Sub-humid',   'R5- Semi-arid'};

%% map of study area
figure,
set(gcf, 'Position', [5 1 24 10]);

% idx in map
[x,y]  = LatLon2PixelCoord(lon,lat,90,-180,1);
idxMap = sub2ind([180,360], y, x);

colKG       = [rgb('DarkCyan');rgb('YellowGreen');rgb('DarkGreen');rgb('Olive');rgb('Gold')];
colLimKG    = [0.5 5.5];
colLabelKG.cticks  = [1:1:5];

data        = zoneIdx;
Z           = NaN(180, 360); % prepare mapgrid
Z           = imbedm(lat, lon, data, Z, geoRaRef); % insert data in mapgrid
Z           = [Z ones(180,1)];  % add column for correct plotting of last column
Z           = [Z; ones(1,361)]; % add row for correct plotting of last row

axes(ha(1));
worldmap('World');
ha(1).Layer = 'top';
setm(gca,'ParallelLabel', 'off','MeridianLabel','off', 'FLineWidth', 0.5);
geoshow(ha(1), land, 'FaceColor', [0.85 0.85 0.85])
hold on
sm = surfm([-90 90], [-180 180], Z);
geoshow(ha(1), rivers, 'Color', [.25 .25 .25]);
colormap(ha(1),colKG);
caxis(ha(1),colLimKG);
% t   =   title('Regions of the Study Area','Fontsize', 10);

ha(1).Position= [0.001 0.5 0.49 0.6]

% add points for opti grids
sm2 = scatterm(opti.lat, opti.lon, 2*ones(size(opti.lat)),rgb('DarkRed'), 'filled');

cb  =   colorbar;
set(cb,'FontSize',8, 'Position', [0.05   0.95   0.4    0.007], 'Orientation', 'horizontal', 'Ticks', colLabelKG.cticks, 'TickLabels', zoneNames, 'AxisLocation', 'in');

keepLabels = ~cellfun(@isempty,cb.TickLabels);
cb.TickLabels = [];
ax = axes('Position', cb.Position,...
    'Color', 'none',...
    'YTick', [],...
    'XLim', cb.Limits,...
    'XTick', cb.Ticks(keepLabels),...
    'XTickLabel', {'Cold', 'Temperate', 'Humid', 'Sub-humid', 'Semi-arid'},...
    'FontSize', cb.FontSize);
xtickangle(ax, 35);

% boxplots
axes(ha(3))
h = histogram(KG_v)%,
b = bar(1:5, h.Values);
b.CData = colKG;
b.FaceColor = 'flat'
set(ha(3), 'Box', 'off', 'Color', 'none', 'XLim', [0.5 5.5], 'Fontsize', 8);
ha(3).YAxis.Visible = 'off';
ha(3).Position = [0.04 0.75 0.1 0.04];

an = annotation('textbox', 'String', '\bullet calibration grid cells', 'fontsize', 9,  'edgecolor', 'none', 'Position', [0.15 0.66 0.3 0.05], 'HorizontalAlignment', 'left', 'Color', rgb('DarkRed') );

print(gcf,[spth 'Fig2_StudyArea.png'],'-dpng','-r300');
