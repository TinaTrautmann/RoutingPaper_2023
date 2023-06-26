%% S4 Comparison with EartH2Observe models (Schellekens et al. 2017)
% by Tina Trautmann, Feb. 2022

%% set save pth
spth = [pwd '/figures/supplement/'];
if ~exist(spth, 'dir'), mkdir(spth),end

%% download and save the EartH2Observe data as .nc in one single folder, e.g.
try
    lp_data         = 'M:/people/uweber/_data/4TT/EarthH2O/Data/'
catch
    error('please  provide the folder path of the EartH2Observe model results. if needed download the data from  the Water Cycle Integrator portal (WCI, wci.earth2observe.eu)')
end


%% load wTWS & info of VEG_org
space_run = 'validCali';
expName   = 'E_B_bL_RD4';
rDate     = '20200629';

pth_in  = [pwd '/data/model_output/VEG/' ];
fName   = [pth_in space_run '_' expName '_' space_run '_wTotal.nc'];
mod.MOD.wTotal  =    ncread(fName,'wTotal');
infoIn.MOD      = load([pth_in space_run '_' expName '_' rDate '_info.mat']);
info            = infoIn.MOD.info;

% space stuff
pthIn   = [pwd '/data/input/studyArea/globalBaseline_Constraints_1deg.mat']
load(pthIn, 'lat', 'lon', 'time')

pix_a = AreaGridLatLon(lat,lon,[1 1]);
pix_a = pix_a(:,1);

[pix_x,pix_y] = LatLon2PixelCoord(lon,lat,90,-180,1);
idx = sub2ind([180, 360],pix_y, pix_x);

% time stuff
[xDay,  ~, ~, ~, ~, Md,D]   = createDateVector(info.tem.model.time.sDate,info.tem.model.time.eDate, 'd');
[xMonth, ~, ~, ~, ~, M]     = createDateVector(info.tem.model.time.sDate,info.tem.model.time.eDate, 'm');

tick_locations  = xMonth(M==1);
tick_locationsD = xDay(D==1);

nPix   = size(pix_a,1);
nTix   = size(xDay,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read the all the Earth2Observe models
% Directory
var_name_all    = {'CanopInt','GroundMoist','RootMoist','SurfMoist','SurfStor','SWE','TotMoist'}; % [mm/s] as average over past 3-hours %name in netcdf ; , ,'SWEVeg'
var_save_all    = {'CW','GW','RW','SM','SW','SWE','TM'}; % name that shall be saved ,'SWEVeg'

time_e20    = {'2002-01-01', '2012-12-31'};
[xMonth_e20, ~, ~, ~, ~, M]   = createDateVector(time_e20{1}, time_e20{2}, 'm');

% Time of data
time_in = {'1979-01-01', '2012-12-31'}

% Overlapping Time
[idx1, idx_in, time3, xData3] = calcConsisTime(time_e20,time_in, 'monthly');
[Y,M,D]     = datevec(xData3);

files_data      = dir(lp_data);
% Loop over files (each including data for 1 model + 1 variable)
for j=3:length(files_data)
    % read nc
    single_data = files_data(j).name;
    info_tmp    = ncinfo(strcat(lp_data,single_data));
    tmp         = strsplit(single_data,'_');
    %get institute
    inst        = char(tmp(2))
    %get variable
    var_name    = char(info_tmp.Variables(1).Name)
    data_tmp    = ncread(strcat(lp_data,single_data), var_name);
    % get timesteps of file
    ts              = size(data_tmp,3);
    data_allts      = NaN(nPix,ts);
    % loop over timesteps & extract data for valid pixel
    cnt = 1;
    for d=1:ts
        data_map           = data_tmp(:,:,d)';
        data_allts(:,cnt)  = data_map(idx);
        cnt  = cnt+1;
    end
    % extract relevant months
    data = data_allts(:,idx_in);
    
    for vv=1:numel(var_name_all)
        if strcmp(var_name,var_name_all{vv})==1
            var_save = var_save_all{vv}
        end
    end
    % put variable in structure of institute
    eval(char(['h2o.' inst '.' var_save ' = data;']))
    
    % Plot Map 
    PlotMapGlobal_noInfo([var_save ' - ' inst],'mean','',lat, lon, mean(data,2),[],jet,[],1,1)
end

% loop over all institutes and sum up to TWS
insts   = fieldnames(h2o);
TWS_all = zeros(nPix,length(xData3));

for ii=1:numel(insts)
    stors   = fieldnames(h2o.(insts{ii}))
    TWS     = zeros(size(h2o.(insts{ii}).(stors{1})));
    tws_int = [];
    h2o.(insts{ii}).tws_info = ' ';
    %loop over storages
    for ss=1:numel(stors)
        %sum up: TWS = SWE + CW + TM + SW + GW (add SWEVeg as well?)
        if strcmp(stors{ss},'SWE')==1 || strcmp(stors{ss},'CW')==1 || strcmp(stors{ss},'TM')==1 || strcmp(stors{ss},'SW')==1 || strcmp(stors{ss},'GW')==1
            TWS     = TWS + h2o.(insts{ii}).(stors{ss});
            tws_int = [tws_int,'+', stors{ss}]; %string to track which storages are included
            stors{ss}
        end
    end
    h2o.(insts{ii}).wTotal = TWS;
    h2o.(insts{ii}).tws_info = tws_int;
        
    
    % average of all institutes
    TWS_all = cat(3, TWS_all, TWS);
        
    clear TWS tws_int TWSa m 
end
h2o.meanAll.wTotal = mean(TWS_all,3);

clear data

% Plot Map
figure
PlotMapGlobal_noInfo('avg TWS of all EartH2O models','2002-2012','',lat, lon, mean(mean(TWS_all,3),2),[],jet,[],1,1)

% get consistent time
[idx_e2o, idx_mo, time3, xMonth3] = calcConsisTime(time_e20, {info.tem.model.time.sDate info.tem.model.time.eDate}, 'monthly');

%% get the zones
load('data/input/ancillary/clusterRegions.mat', 'CLregions');

KG_map   = CLregions;

KG_v        = KG_map(idx); %valids in map
KG_v(KG_v==0) = 4; 
KG_v(KG_v==4) = 10; 
KG_v(KG_v==2) = 20; 
KG_v(KG_v==3) = 30; 
KG_v(KG_v==1) = 40; 
KG_v(KG_v==5) = 50; 
KG_v = KG_v ./ 10;

zonesID     = unique(KG_v);
zoneNames   = {'R1- Cold', 'R2- Temperate', 'R3- Humid',  'R4- Sub-humid',   'R5- Semi-arid'};

TWScat = 1;

%% load obs
obs     = load(pthIn, 'TWS_GRACE', 'TWS_GRACE_unc');
oData   = obs.TWS_GRACE(:,idx_mo);
if TWScat == 1
    oData(oData<=-500) = NaN;
    oData(oData>=500)  = NaN;
end
m       = nanmean(oData,2);
oData   = oData - repmat(m,1,length(xMonth3));
oData_MSC   = calcMSC(oData, M);

uData = obs.TWS_GRACE_unc(:,idx_mo);
if TWScat == 1
    uData(oData<=-500) = NaN;
    uData(oData>=500)  = NaN;
end
uData_MSC   = calcMSC(uData, M);

% put into structure
dataObs.wTWS     = oData;
dataObsMSC.wTWS  = oData_MSC;

dataUnc.wTWS     = uData;
dataUncMSC.wTWS  = uData_MSC;

%% prepare mod
expNames    = ['MOD', fieldnames(h2o)'];

NaNidx  = find(isnan(dataObs.wTWS));

for eN=1:numel(expNames)
    expN    = expNames{eN};
    if strcmp(expN,'MOD')
        sData_d = mod.(expN).wTotal;   
        sData   = aggDay2Mon(sData_d, info.tem.model.time.sDate,info.tem.model.time.eDate);
        sData   = sData(:,idx_mo);
    else
        sData   = h2o.(expN).wTotal;
    end
    
    % monthly aggregation
    % only consistent data points
    sData(NaNidx) = NaN;
    
    % anomaly if TWS
    m = nanmean(sData,2);
    sData = sData - repmat(m,1,length(xMonth3));
    
    % calculate MSC
    sData_MSC   = calcMSC(sData, M);
    
    % put into structure
    data.(expN).wTWS      = sData;
    dataMSC.(expN).wTWS   = sData_MSC;
end

dataNames = ['obs' expNames];

% GRACE vs wTWS
dataXX = {};
for op=1:numel(expNames)
    dataXX = [dataXX  dataMSC.(expNames{op})];
end
obsTmp.wTWS      = dataObsMSC.wTWS;


tmpExpNames = fieldnames(dataMSC);
tmpObs      = dataObsMSC.wTWS;
for fN=1:numel(tmpExpNames)
    expName = tmpExpNames{fN};
    mefMSC.(expName) = NaN(nPix,1);
    for nP=1:nPix
        mefMSC.(expName)(nP,1) = calcMEF(tmpObs(nP,:),dataMSC.(expName).wTWS(nP,:), ones(size(tmpObs(nP,:))));
    end
end


% colXX  = [rgb('Red'); rgb('green'); rgb('Blue');  rgb('DarkMagenta'); rgb('GoldenRod'); rgb('DimGray')];
% unitNames_tmp     = {'mm'};
% 
% tName       = 'Gridwise MEF with GRACE-TWS';
% Lim.colLim  = [-1 1];
% Lim.colLimDiff  = [-0.8 0.8];
% Lim.col       = othercolor('RdBu4',100);
% Lim.colDiff   = othercolor('PuOr11',201);
% 
% tmpExpNames1 = tmpExpNames(1:5)
% [sp] = plot_CombinedMSC_MEFmaps3_e2o(tmpExpNames1, obsTmp, dataXX, mefMSC.MOD, mefMSC, lat, lon, zoneNames, KG_v, pix_a, unitNames_tmp, colXX, Lim)
% print(gcf,[spth 'combinedPlot3b_MEF_1.png'],'-dpng','-r300');
% 
% tmpExpNames2 = [tmpExpNames(1); tmpExpNames(6:9)]
% [sp] = plot_CombinedMSC_MEFmaps3_e2o(tmpExpNames2, obsTmp, dataXX, mefMSC.MOD, mefMSC, lat, lon, zoneNames, KG_v, pix_a, unitNames_tmp, colXX, Lim)
% print(gcf,[spth 'combinedPlot3b_MEF_2.png'],'-dpng','-r300');
% 
% tmpExpNames3 = [tmpExpNames(1); tmpExpNames(10:11)]
% [sp] = plot_CombinedMSC_MEFmaps3_e2o(tmpExpNames3, obsTmp, dataXX, mefMSC.MOD, mefMSC, lat, lon, zoneNames, KG_v, pix_a, unitNames_tmp, colXX, Lim)
% print(gcf,[spth 'combinedPlot3b_MEF_3.png'],'-dpng','-r300');
% 
%% Boxplots of MEF for regions
colEXP = othercolor('Paired9', numel(expNames));

% for boxplot dataTmp = [zones, exp, values]
dataTmp = NaN(numel(zoneNames)+1,numel(expNames),length(pix_a));
for cc=1:numel(expNames)
    data = mefMSC.(expNames{cc});
    dataTmp(1,cc,:)=data;
    for zN=1:numel(zoneNames)
        dataTmp(1+zN,cc,1:length(find(KG_v==zN)))=data(KG_v==zN);
    end
end

figure
set(gcf, 'Position', [5 1 24 10]);

zoneNames2 = {'R1', 'R2', 'R3', 'R4', 'R5'}             
xTmp = 1:length(zoneNames)+1;
p = boxplot2(dataTmp,xTmp);
for ii = 1:numel(expNames)
    structfun(@(x) set(x(ii,:), 'color', colEXP(ii,:), ...
        'markeredgecolor', colEXP(ii,:)), p);
end
ha=gca
ha.XTickLabel = ['Global', zoneNames2];
ha.XLim = [0.5 length(zoneNames)+1.5];
ha.Box     = 'on';
ha.YTickLabelMode = 'auto';
ha.YLim = [-1 1]

set([p.box],'linewidth', 1.5);
set([p.med],'linewidth', 1.5);
set([p.lwhis p.uwhis], 'linestyle', '-', 'linewidth', 1.5);
set(p.out, 'marker', '.');
hold on,
rf=refline(0,0);
rf.Color =  rgb('Black');
title('Gridwise MEF of MOD & EartH2Observe Models')
legend(expNames, 'Location', 'EastOutside')
print(gcf,[spth 'S4a_Comparison_TWS-MEF_EartH2ObserveEnsemble_gridwiseBoxplots.png'],'-dpng','-r300');

close all

% MEF of average MSC for H2O models
avgMEF = NaN(6,numel(expNames));

tmpExpNames = fieldnames(dataMSC);
tmpObs      = dataObsMSC.wTWS;
for fN=1:numel(tmpExpNames)
    expName = tmpExpNames{fN};
    avgMEF(1,fN) = calcMEF(nanmean(tmpObs,1),nanmean(dataMSC.(expName).wTWS,1), ones(size(1,12)));
    for zn=2:6
        idxZ = find(KG_v==zn-1);
        avgMEF(zn,fN) = calcMEF(nanmean(tmpObs(idxZ,:),1),nanmean(dataMSC.(expName).wTWS(idxZ,:),1), ones(size(1,12)));
    end
end

figure
bar(-2,-2, 'FaceColor', rgb('GoldenRod')) %for legend
hold on, plot(-1,-1, '-', 'color', [.1 .1 .1]) %for legend
hold on, plot(-1,-1, '*', 'color', [.1 .1 .1]) %for legend
hold on

boxplot(avgMEF(:,2:end-1)', 'colors', rgb('GoldenRod'))

h = findobj('Tag', 'Box')
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),rgb('GoldenRod'),'FaceAlpha',.5, 'EdgeColor',rgb('GoldenRod'));
end
for j=1:length(h)
    jj = 7-j;
    tmp = get(h(jj),'XData');
    hold on, line([tmp(1) tmp(3)],[nanmean(avgMEF(j,2:end-1)) nanmean(avgMEF(j,2:end-1))], 'LineStyle', '-', 'Color', 'k','Linewidth', .75) %model mean
    hold on, plot(j, avgMEF(j,1), '*', 'color', [.1 .1 .1], 'Linewidth', 1) % our model
end
set(gca, 'Position', [0.15 0.25 0.74 0.68])
set(gca, 'XLim', [0.5 6.5], 'XTick', [1:1:6], 'XTickLabel', ['global' zoneNames2])
set(gca, 'YLim', [0 1])
ylabel('MEF with GRACE TWS', 'Fontsize', 8)

l=legend({'eartH2Observe model ensemble','model mean','MOD of this study'}, 'Orientation', 'Horizontal')
legend('boxoff')
set(l, 'fontsize', 8, 'Position', [0.24 0 0.56 0.14])

title('MEF of MOD & EartH2Observe Models')
print(gcf,[spth 'S5_Comparison_TWS-MEF_EartH2ObserveEnsemble.png'],'-dpng','-r300');
