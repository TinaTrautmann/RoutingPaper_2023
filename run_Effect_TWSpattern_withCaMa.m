%% Compare org GRACE and TWS_noRiver & TWS_noRiver_CaMa
% by Tina Trautmann, Aug. 2022

%% set save pth
spth = [pwd '/figures/'];
if ~exist(spth, 'dir'), mkdir(spth),end


%% load data
pthIn = [pwd '/data/input/studyArea/globalBaseline_Constraints_1deg.mat']

load(pthIn, 'lat', 'lon', 'time')

obs     = load(pthIn, 'TWS_GRACE', 'TWS_GRACE_unc');
noRiver = load(pthIn, 'TWS_noRiver_01', 'TWS_noRiver_05', 'TWS_noRiver_25', 'TWS_noRiver_CaMa');

%  the lat & lon of opti
opti = load([pwd '/data/input/opti/lat_lon_904.mat'], 'lat', 'lon');

%% Prep Plotting Stuff
TWScat = 1; %only consider TWSobs >= -500mm and <= 500mm

% space stuff
pix_a = AreaGridLatLon(lat,lon,[1 1]);
pix_a = pix_a(:,1);

% time stuff
[xMonth, ~, ~, ~, ~, M]     = createDateVector(time{1},time{2}, 'm');

tick_locations  = xMonth(M==1);

nPix   = size(pix_a,1);
nTix   = size(xMonth,2);

% colors
colLabel = [];

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



%% prepare the obs
oData = obs.TWS_GRACE;

% anomaly for TWS
if TWScat == 1
    oData(oData<=-500) = NaN;
    oData(oData>=500)  = NaN;
end
m     = nanmean(oData,2);
oData = oData - repmat(m,1,length(xMonth));

% MSC
oData_MSC   = calcMSC(oData, M);

% put into structure
dataObs.GRACE_TWS      = oData;
dataObsMSC.GRACE_TWS   = oData_MSC;
    
% uData = obs.TWS.unc;
uData = obs.TWS_GRACE_unc;
if TWScat == 1
    uData(oData<=-500) = NaN;
    uData(oData>=500)  = NaN;
end
uData_MSC   = calcMSC(uData, M);

% put into structure
dataUnc.wTWS     = uData;
dataUncMSC.wTWS  = uData_MSC;


%% prepare the TWS_noRiver
names       = {'TWS_noRiver_01', 'TWS_noRiver_05', 'TWS_noRiver_25', 'TWS_noRiver_CaMa'};
oDatas      = {'noRiver.TWS_noRiver_01', 'noRiver.TWS_noRiver_05', 'noRiver.TWS_noRiver_25', 'noRiver.TWS_noRiver_CaMa'};
for ii=1:numel(names)
    name  = names{ii};
    oData = eval([oDatas{ii}]);
    
    % anomaly for TWS
    if TWScat == 1
        oData(oData<=-500) = NaN;
        oData(oData>=500)  = NaN;
    end
    m     = nanmean(oData,2);
    oData = oData - repmat(m,1,length(xMonth));
    
    % MSC
    oData_MSC   = calcMSC(oData, M);
    
    % put into structure
    dataObs.(name)      = oData;
    dataObsMSC.(name)   = oData_MSC;
end

%% plot the maps
tmpExpNames = fieldnames(dataObsMSC);
tmpObs      = dataObsMSC.GRACE_TWS;
for fN=1:numel(tmpExpNames)
    expName = tmpExpNames{fN};
    mefMSC.(expName) = NaN(nPix,1);
    for nP=1:nPix
        mefMSC.(expName)(nP,1) = calcMEF(tmpObs(nP,:),dataObsMSC.(expName)(nP,:), ones(size(tmpObs(nP,:))));
    end
end


dataNames   = {'GRACE_TWS', 'TWS_noRiver_25', 'TWS_noRiver_05', 'TWS_noRiver_01', 'TWS_noRiver_CaMa'};
names       = {'TWS_noRiver_25', 'TWS_noRiver_05', 'TWS_noRiver_01', 'TWS_noRiver_CaMa'};
dataObs_tmp.wTWS_obs     = dataObsMSC.GRACE_TWS;

dataXX = {};
for op=1:numel(names)
    tmp.wTWS_obs = dataObsMSC.(names{op});
    dataXX = [dataXX tmp];
    clear tmp
end

mefObs      = mefMSC.GRACE_TWS;
mefMod      = mefMSC;

colEXP      = [rgb('Black'); rgb('Chocolate');  rgb('DarkMagenta'); rgb('Blue'); rgb('Crimson')];
unitNames   = {'mm'};
xSeason     = 1:1:12;


varNames  = fieldnames(dataXX{1});
expNames  = strrep(dataNames,'_','-');

zoneNames  = strrep(zoneNames,'_','-');
zoneNames2 = {'R1', 'R2', 'R3', 'R4', 'R5'};
zNames     = ['Global', zoneNames];
zNames2    = strrep(zNames,'-','');
zNames2    = strrep(zNames2,' ','');
zoneIdx    = KG_v;



% preps maps
colLim      = [0 1];
colLimDiff  = [0 0.1];

col         = othercolor('RdBu4',100);
colDiff     = othercolor('YlOrRd9',201); %othercolor('PuOr11',201);

land        = shaperead('landareas', 'UseGeoCoords', true);
rivers      = shaperead('worldrivers', 'UseGeoCoords', true);
geoRaRef    = georasterref('RasterSize', [180 360], 'RasterInterpretation', 'cells',  ...
    'LatitudeLimits', [-90 90], 'LongitudeLimits', [-180 180]);
Z1           = NaN(180, 360); % prepare mapgrid

% for boxplot dataTmp = [zones, exp, values]
dataTmp = NaN(numel(zoneNames)+1,numel(expNames),length(pix_a));

% metrics -> 3 tables, 1 per variable
tmpArray    = NaN(numel(dataNames)-1,numel(zNames));
varTypes    = cell(numel(zNames),1);
varTypes(:) = {'double'};
for vn=1:numel(varNames)
    T_corr.(varNames{vn}) = array2table(tmpArray,'VariableNames', zNames2, 'RowNames', dataNames(2:end));
    T_MEF.(varNames{vn})  = array2table(tmpArray,'VariableNames', zNames2, 'RowNames', dataNames(2:end));
    T_RMSE.(varNames{vn}) = array2table(tmpArray,'VariableNames', zNames2, 'RowNames', dataNames(2:end));
    T_corrMean.(varNames{vn}) = array2table(tmpArray,'VariableNames', zNames2, 'RowNames', dataNames(2:end));
    T_MEFMean.(varNames{vn})  = array2table(tmpArray,'VariableNames', zNames2, 'RowNames', dataNames(2:end));
end

%% FIGURE
nrows   = 4;
ncols   = 3;
eN      = numel(dataXX);
cnt = 1;
sp  = figure;
set(gcf, 'Position', [5 0 16 22]);
ha  = tight_subplot(nrows,ncols,[.05 .02],[.08 .08],[.05 .02]);

% A) TWS-noRiver 2
dataTmp(1,1,:)=mefObs;
for zN=1:numel(zoneNames)
    dataTmp(1+zN,1,1:length(find(zoneIdx==zN)))=mefObs(zoneIdx==zN);
end

% 
cc = 1 
rr = 1

data        = mefObs-mefMod.(dataNames{2});
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

colormap(ha(1),colDiff);
caxis(ha(1),colLimDiff);

title([expNames{2} ],'Fontsize', 10)

tmpExp1 = squeeze(dataTmp(1,1,:));
dataTmp(1,2,:)= tmpExp1 - mefMod.(dataNames{2}); % MEF of 1 - MEF of this
for zN=1:numel(zoneNames)
    tmpExp1 = squeeze(dataTmp(zN+1,1,:));
    dataTmp(1+zN,2,1:length(find(zoneIdx==zN)))=tmpExp1(1:length(find(zoneIdx==zN)),:) - mefMod.(dataNames{2})(zoneIdx==zN);
end

ha(1).Position=  [0.001 0.5 0.49 0.6]

% colorbar for difference
cb  =   colorbar;
set(cb,'FontSize',8, 'Position', [0.2   0.95   0.6    0.007], 'Orientation', 'horizontal','Ticks', [0:0.025:0.1], 'TickLabels', {'0', '', '0.05', '', '0.1'}, 'AxisLocation', 'in');
cb1 =   ylabel(cb, 'MEF difference', 'FontSize', 9);
cb1.Position(1) = 4;

delete(ha(3)) %!!!

% TWS-noRiver 3
cc = 2

data        = mefObs-mefMod.(dataNames{3});
Z           = NaN(180, 360); % prepare mapgrid
Z           = imbedm(lat, lon, data, Z, geoRaRef); % insert data in mapgrid
Z           = [Z ones(180,1)];  % add column for correct plotting of last column
Z           = [Z; ones(1,361)]; % add row for correct plotting of last row

axes(ha(2));
worldmap('World');
ha(2).Layer = 'top';
setm(gca,'ParallelLabel', 'off','MeridianLabel','off', 'FLineWidth', 0.5);
geoshow(ha(2), land, 'FaceColor', [0.85 0.85 0.85])
hold on
sm = surfm([-90 90], [-180 180], Z);
geoshow(ha(2), rivers, 'Color', [.25 .25 .25]);

colormap(ha(2),colDiff);
caxis(ha(2),colLimDiff);

title([expNames{3} ],'Fontsize', 10)

tmpExp1 = squeeze(dataTmp(1,1,:));
dataTmp(1,2,:)= tmpExp1 - mefMod.(dataNames{3}); % MEF of 1 - MEF of this
for zN=1:numel(zoneNames)
    tmpExp1 = squeeze(dataTmp(zN+1,1,:));
    dataTmp(1+zN,2,1:length(find(zoneIdx==zN)))=tmpExp1(1:length(find(zoneIdx==zN)),:) - mefMod.(dataNames{3})(zoneIdx==zN);
end

ha(2).Position=  [0.501 0.5 0.49 0.6]

% % colorbar for difference
% cb  =   colorbar;
% set(cb,'FontSize',8, 'Position', [0.55   0.95   0.4    0.007], 'Orientation', 'horizontal','Ticks', [0:0.025:0.1], 'TickLabels', {'0', '', '', '', '0.1'}, 'AxisLocation', 'in');
% cb1 =   ylabel(cb, 'MEF difference', 'FontSize', 9);
% cb1.Position(2) = 4;
% 
% TWS-noRiver 4 
rr = 2
a_tmp = (rr-1)*ncols;
delete(ha(6)) %!!!

cc = 1      
data        = mefObs-mefMod.(dataNames{4});
Z           = NaN(180, 360); % prepare mapgrid
Z           = imbedm(lat, lon, data, Z, geoRaRef); % insert data in mapgrid
Z           = [Z ones(180,1)];  % add column for correct plotting of last column
Z           = [Z; ones(1,361)]; % add row for correct plotting of last row

axes(ha(4));
worldmap('World');
ha(4).Layer = 'top';
setm(gca,'ParallelLabel', 'off','MeridianLabel','off', 'FLineWidth', 0.5);
geoshow(ha(4), land, 'FaceColor', [0.85 0.85 0.85])
hold on
sm = surfm([-90 90], [-180 180], Z);
geoshow(ha(4), rivers, 'Color', [.25 .25 .25]);

colormap(ha(4),colDiff);
caxis(ha(4),colLimDiff);

title([expNames{4} ],'Fontsize', 10)

tmpExp1 = squeeze(dataTmp(1,1,:));
dataTmp(1,3,:)= tmpExp1 - mefMod.(dataNames{4}); % MEF of 1 - MEF of this
for zN=1:numel(zoneNames)
    tmpExp1 = squeeze(dataTmp(zN+1,1,:));
    dataTmp(1+zN,3,1:length(find(zoneIdx==zN)))=tmpExp1(1:length(find(zoneIdx==zN)),:) - mefMod.(dataNames{4})(zoneIdx==zN);
end

ha(4).Position= [0.001 0.26 0.49 0.6]

% TWS-noRiver 5
cc = 2
data        = mefObs-mefMod.(dataNames{5});
Z           = NaN(180, 360); % prepare mapgrid
Z           = imbedm(lat, lon, data, Z, geoRaRef); % insert data in mapgrid
Z           = [Z ones(180,1)];  % add column for correct plotting of last column
Z           = [Z; ones(1,361)]; % add row for correct plotting of last row

axes(ha(5));
worldmap('World');
ha(5).Layer = 'top';
setm(gca,'ParallelLabel', 'off','MeridianLabel','off', 'FLineWidth', 0.5);
geoshow(ha(5), land, 'FaceColor', [0.85 0.85 0.85])
hold on
sm = surfm([-90 90], [-180 180], Z);
geoshow(ha(5), rivers, 'Color', [.25 .25 .25]);

colormap(ha(5),colDiff);
caxis(ha(5),colLimDiff);

title([expNames{5} ],'Fontsize', 10)

tmpExp1 = squeeze(dataTmp(1,1,:));
dataTmp(1,4,:)= tmpExp1 - mefMod.(dataNames{5}); % MEF of 1 - MEF of this
for zN=1:numel(zoneNames)
    tmpExp1 = squeeze(dataTmp(zN+1,1,:));
    dataTmp(1+zN,4,1:length(find(zoneIdx==zN)))=tmpExp1(1:length(find(zoneIdx==zN)),:) - mefMod.(dataNames{5})(zoneIdx==zN);
end
ha(5).Position= [0.501 0.26 0.49 0.6]



% B) MSC for global + zones
rr=4;
a_tmp = 6;

varN  = varNames{1};
unitN = unitNames{1};
dObs = dataObs_tmp.(varN);
for dN = 1:eN
    tmp = dataXX{dN};
    eval(char(['d' num2str(dN) '= tmp.(varN);']));
end

for cc=1:3
    if cc==1
        idxZ = 1:1:size(dObs,1);
    else
        idxZ = find(zoneIdx==cc-1);
    end
    axes(ha(a_tmp+cc));
    
    plot(xSeason,nanmeanArea(dObs(idxZ,:),pix_a(idxZ)), 'Linewidth', 1.75, 'LineStyle', ':', 'Color', colEXP(1,:)), hold on
    
    for dN = 1:eN
        eval(char(['dXX = d' num2str(dN) ';']))
        plot(xSeason,nanmeanArea(dXX(idxZ,:),pix_a(idxZ)), 'Linewidth', 1.25, 'LineStyle', '-', 'Color', colEXP(dN+1,:)), hold on
    end
    plot(xSeason,nanmeanArea(dObs(idxZ,:),pix_a(idxZ)), 'Linewidth', 1.75, 'LineStyle', ':', 'Color', colEXP(1,:)), hold on
    plot([0.5 12.5],[0 0], '-', 'color', rgb('Black'))
    set(gca, 'Fontsize', 9)
    title(zNames{cc}, 'FontSize',10)
    if cc==1; yl= ylabel(['TWS [' unitN ']'], 'Fontweight', 'b', 'Fontweight', 'b', 'FontSize', 10); yl.Position(1) = -1.5; end
    ha(cnt).XAxis.MinorTickValues =  [1:1:12];
    ha(a_tmp+cc).XLim = [0.5 12.5];
    
    ha(a_tmp+cc).Position(1) = ha(a_tmp+cc).Position(1)  + 0.03;
    ha(a_tmp+cc).Position(2) = 0.26;
    ha(a_tmp+cc).Position(3) = 0.25;
    ha(a_tmp+cc).Position(4) = 0.15;
    
    %calculate metrics
    tmpO = dObs(idxZ,:);
    for dN = 1:eN
        eval(char(['dXX = d' num2str(dN) ';']))
        tmp1 = dXX(idxZ,:);
        T_corr.(varN){dN,cc} = round(corr(tmpO(:),tmp1(:), 'rows', 'complete'),2);
        T_MEF.(varN){dN,cc}  = round(calcMEF(tmpO(:),tmp1(:),ones(size(tmp1(:)))),2);
        T_RMSE.(varN){dN,cc} = round(calcRMSE(tmpO(:),tmp1(:)),0);
        T_corrMean.(varN){dN,cc} = round(corr(nanmeanArea(dObs(idxZ,:),pix_a(idxZ))', nanmeanArea(dXX(idxZ,:),pix_a(idxZ))', 'rows', 'complete'),2);
        T_MEFMean.(varN){dN,cc} = round(calcMEF(nanmeanArea(dObs(idxZ,:),pix_a(idxZ)), nanmeanArea(dXX(idxZ,:),pix_a(idxZ)), ones(size(pix_a(idxZ),12))),2);
    end
    
    %add the correlation in the plot
    for dN = 1:eN
        %             t1 = text(0.97,-0.025+0.1*dN, ['r^2 = ' num2str(T_corr.(varN){dN,cc})],'Units','Normalized', 'Fontsize', 6, 'Color', col(dN+1,:), 'HorizontalAlignment', 'right');
        if cc == 5
            t1 = text(0.97,-0.05+0.1*dN, ['MEF = ' num2str(T_MEFMean.(varN){dN,cc})],'Units','Normalized', 'Fontsize', 8, 'Color', colEXP(dN+1,:), 'HorizontalAlignment', 'right');
        elseif cc == 2 || cc == 3
            t1 = text(0.03,-0.05+0.1*dN, ['MEF = ' num2str(T_MEFMean.(varN){dN,cc})],'Units','Normalized', 'Fontsize', 8, 'Color', colEXP(dN+1,:), 'HorizontalAlignment', 'left');            
        else
            t1 = text(0.97,0.55+0.1*dN, ['MEF = ' num2str(T_MEFMean.(varN){dN,cc})],'Units','Normalized', 'Fontsize', 8, 'Color', colEXP(dN+1,:), 'HorizontalAlignment', 'right');
        end
    end
end


rr=5;
a_tmp = 9;

for cc=4:6
    if cc==1
        idxZ = 1:1:size(dObs,1);
    else
        idxZ = find(zoneIdx==cc-1);
    end
    axes(ha(a_tmp+cc-3));
    
    plot(xSeason,nanmeanArea(dObs(idxZ,:),pix_a(idxZ)), 'Linewidth', 1.75, 'LineStyle', ':', 'Color', colEXP(1,:)), hold on
    
    for dN = 1:eN
        eval(char(['dXX = d' num2str(dN) ';']))
        plot(xSeason,nanmeanArea(dXX(idxZ,:),pix_a(idxZ)), 'Linewidth', 1.25, 'LineStyle', '-', 'Color', colEXP(dN+1,:)), hold on
    end
    plot(xSeason,nanmeanArea(dObs(idxZ,:),pix_a(idxZ)), 'Linewidth', 1.75, 'LineStyle', ':', 'Color', colEXP(1,:)), hold on
    plot([0.5 12.5],[0 0], '-', 'color', rgb('Black'))
    set(gca, 'Fontsize', 9)
    title(zNames{cc}, 'FontSize',10)
    if cc==4; yl= ylabel(['TWS [' unitN ']'], 'Fontweight', 'b', 'Fontweight', 'b', 'FontSize', 10); yl.Position(1) = -1.5; end
    ha(a_tmp+cc-3).XAxis.MinorTickValues =  [1:1:12];
    ha(a_tmp+cc-3).XLim = [0.5 12.5];
    
    ha(a_tmp+cc-3).Position(1) = ha(a_tmp+cc-3).Position(1)  + 0.03;
    ha(a_tmp+cc-3).Position(2) = 0.06;
    ha(a_tmp+cc-3).Position(3) = 0.25;
    ha(a_tmp+cc-3).Position(4) = 0.15;
    
    %calculate metrics
    tmpO = dObs(idxZ,:);
    for dN = 1:eN
        eval(char(['dXX = d' num2str(dN) ';']))
        tmp1 = dXX(idxZ,:);
        T_corr.(varN){dN,cc} = round(corr(tmpO(:),tmp1(:), 'rows', 'complete'),2);
        T_MEF.(varN){dN,cc}  = round(calcMEF(tmpO(:),tmp1(:),ones(size(tmp1(:)))),2);
        T_RMSE.(varN){dN,cc} = round(calcRMSE(tmpO(:),tmp1(:)),0);
        T_corrMean.(varN){dN,cc} = round(corr(nanmeanArea(dObs(idxZ,:),pix_a(idxZ))', nanmeanArea(dXX(idxZ,:),pix_a(idxZ))', 'rows', 'complete'),2);
        T_MEFMean.(varN){dN,cc} = round(calcMEF(nanmeanArea(dObs(idxZ,:),pix_a(idxZ)), nanmeanArea(dXX(idxZ,:),pix_a(idxZ)), ones(size(pix_a(idxZ),12))),2);
    end
    
    %add the correlation in the plot
    for dN = 1:eN
        %             t1 = text(0.97,-0.025+0.1*dN, ['r^2 = ' num2str(T_corr.(varN){dN,cc})],'Units','Normalized', 'Fontsize', 6, 'Color', col(dN+1,:), 'HorizontalAlignment', 'right');
        if cc == 5
            t1 = text(0.97,-0.05+0.1*dN, ['MEF = ' num2str(T_MEFMean.(varN){dN,cc})],'Units','Normalized', 'Fontsize', 8, 'Color', colEXP(dN+1,:), 'HorizontalAlignment', 'right');
        elseif cc == 2 || cc == 3 || cc == 6
            t1 = text(0.03,-0.05+0.1*dN, ['MEF = ' num2str(T_MEFMean.(varN){dN,cc})],'Units','Normalized', 'Fontsize', 8, 'Color', colEXP(dN+1,:), 'HorizontalAlignment', 'left');            
        else
            t1 = text(0.97,0.55+0.1*dN, ['MEF = ' num2str(T_MEFMean.(varN){dN,cc})],'Units','Normalized', 'Fontsize', 8, 'Color', colEXP(dN+1,:), 'HorizontalAlignment', 'right');
        end
    end
end

% set axis
set(ha(7:12), 'Xlim', [0.5 12.5], 'XTick', [3:3:12], 'GridAlpha', 0.25, 'XMinorGrid', 'on', 'MinorGridAlpha', 0.05, 'MinorGridLineStyle', '-',...
    'YGrid', 'on')

for tmp=7:12,ha(tmp).YLabel.FontSize = 10; ha(tmp).YAxis.FontSize = 8;ha(tmp).XAxis.FontSize = 8; end

l  = legend(expNames, 'Position', [.0 0 .99 .03], 'box', 'off', 'Orientation', 'Horizontal','Fontsize', 9, 'NumColumns',3);

% an = annotation('textbox',[0.01 .905 1 .1],'String', 'a)' ,'LineStyle','none', 'Fontsize', 11, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
% an = annotation('textbox',[0.05 .905 1 .1],'String', 'Regions of the Study Area' ,'LineStyle','none', 'Fontsize', 11, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

% an = annotation('textbox',[0.51 .905 1 .1],'String', 'b)' ,'LineStyle','none', 'Fontsize', 11, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
% an = annotation('textbox',[0.55 .905 1 .1],'String', ['Difference in MEF with GRACE TWS'] ,'LineStyle','none', 'Fontsize', 11, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

an = annotation('textbox',[0.01 .905 1 .1],'String', 'a)' ,'LineStyle','none', 'Fontsize', 11, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
an = annotation('textbox',[0.05 .905 1 .1],'String', ['Difference in MEF with GRACE TWS'] ,'LineStyle','none', 'Fontsize', 11, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');


an = annotation('textbox',[0.01 .36 1 .1],'String', 'b)' ,'LineStyle','none', 'Fontsize', 11, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
an = annotation('textbox',[0.05 .36 1 .1],'String', 'Mean Seasonal Dynamics' ,'LineStyle','none', 'Fontsize', 11, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');


print(gcf,[spth 'Fig3_TWSpattern_CaMa.png'],'-dpng','-r300');

