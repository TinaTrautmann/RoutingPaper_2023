%% S1 regional FlowAccumulation
% by Tina Trautmann, Feb. 2022

%% set save pth
spth = [pwd '/figures/supplement/'];
if ~exist(spth, 'dir'), mkdir(spth),end

%% studyArea
% load lat lon
pthIn = [pwd '/data/input/studyArea/globalBaseline_Constraints_1deg.mat']

load(pthIn, 'lat', 'lon')

% idx in map
[x,y]  = LatLon2PixelCoord(lon,lat,90,-180,1);
idxMap = sub2ind([180,360], y, x);


%% read the flow accumulation
data = load([pwd '/data/input/ancillary/RoutingData_1deg.mat'], 'flwacc');
% extract the data
flw_data   = data.flwacc(idxMap);

figure
PlotMapGlobal_noInfo('Flow Accumulation', [], ' ', lat, lon, flw_data, [], othercolor('YlGnBu9',250), [], 1, 1);
print(gcf, [spth 'S1a_FlowAccumulation_studyArea.png'], '-dpng', '-r300')


%% read the cluster zones
%% get the zones
load([pwd '/data/input/ancillary/clusterRegions.mat'], 'CLregions');

KG_map   = CLregions;

KG_v        = KG_map(idxMap); %valids in map
KG_v(KG_v==0) = 4; 
KG_v(KG_v==4) = 10; 
KG_v(KG_v==2) = 20; 
KG_v(KG_v==3) = 30; 
KG_v(KG_v==1) = 40; 
KG_v(KG_v==5) = 50; 
KG_v = KG_v ./ 10;

zonesID     = unique(KG_v);
zoneNames   = {'R1- Cold', 'R2- Temperate', 'R3- Humid',  'R4- Sub-humid',   'R5- Semi-arid'};
colKG       = [rgb('DarkCyan');rgb('YellowGreen');rgb('DarkGreen');rgb('Olive');rgb('Gold')];

%% Plot KG_v vs flw_data
figure, 
boxplot(flw_data,KG_v)
set(gca, 'ylim', [0 25],'XTickLabels', zoneNames)
title('FlowAccumulation per Region')


%5 heatmap
% classify flw_acc
flw_class = flw_data;
flw_class(flw_class>5 & flw_class<10) =6;
flw_class(flw_class>=10) = 7;

tmpT = table(KG_v,flw_class);

figure, 
h = heatmap(tmpT,'KG_v', 'flw_class')
h.YDisplayData = flipud(h.YDisplayData);

%make it proportional to number of grids in this zone
flw_prop = h.ColorData./sum(h.ColorData,1);
figure,
h2=heatmap(round(flw_prop,2))
h2.YDisplayData = flipud(h2.YDisplayData)
h2.XLabel = 'Cluster Regions';
h2.XDisplayLabels = zoneNames;
h2.YLabel = 'Flow Accumulation';
h2.YDisplayLabels = {'>10','>5 -9',  '5', '4', '3', '2', '1'}
colormap(othercolor('YlGnBu9',250))
print(gcf, [spth 'S0_FlowAccumulation_per_region_percentGridCount.png'], '-dpng', '-r300')

