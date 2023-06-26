%% Effect of wRiver on model validation with CaMa
% by Tina Trautmann, Aug. 2022

%% set save pth
spth = [pwd '/figures/'];
if ~exist(spth, 'dir'), mkdir(spth),end


%% load wTWS & info of VEG_org
space_run = 'validCali';
expName   = 'E_B_bL_RD4';
rDate     = '20200629';

pth_in  = [pwd '/data/model_output/VEG/' ];
fName   = [pth_in space_run '_' expName '_' space_run '_wTotal.nc'];
mod.MOD.wTotal    =    ncread(fName,'wTotal');
infoIn.MOD  = load([pth_in space_run '_' expName '_' rDate '_info.mat']);


%% load routed & necessary data
expNames             = {'VEG_25_8', 'VEG_05_5', 'VEG_01_7', 'VEG_CaMa_7'};
expNames_short       = {'VEG_25', 'VEG_05', 'VEG_01', 'VEG_CaMa' };
expNames_paper       = {'MOD_R25', 'MOD_R05', 'MOD_R01', 'MOD_CaMa'};
rDates               = {'20210531','20210531','20210531','20220817'};
space_run            = 'validCali';

varNames = {'wTotal', 'wRiver'}

% load info, wTWS & wRiver from 
for eN=1:numel(expNames)
    expName  = expNames{eN};
    expName2 = expNames_short{eN};
    expName3 = expNames_paper{eN};
    rDate    = rDates{eN};
    % path of model output
    pth_in  = [pwd '/data/model_output/' expName2 '/'];
    infoIn.(expName) = load([pth_in expName '_' rDate '_info.mat']);
    
    % load variables
    for vn=1:numel(varNames)
        vName = varNames{vn};
        fName   =  [pth_in  expName '_' space_run '_' vName '.nc'];
        fName2  =  [pth_in   space_run '_' expName2 '_' space_run '_' vName '.nc'];
        if isfile(fName)
            mod.(expName3).(vName) =    ncread(fName,vName);
        elseif isfile(fName2)
            mod.(expName3).(vName) =    ncread(fName2,vName);
        end
    end
end


%% Prep Plotting
info = infoIn.MOD.info;

TWScat = 1; %only consider TWSobs >= -500mm and <= 500mm

pthIn   = [pwd '/data/input/studyArea/globalBaseline_Constraints_1deg.mat']
obs     = load(pthIn, 'TWS_GRACE', 'TWS_GRACE_unc');

load(pthIn, 'lat', 'lon', 'time')


% space stuff
load([pwd '/data/input/studyArea/globalBaseline_Constraints_1deg.mat'], 'lat', 'lon')

pix_a = AreaGridLatLon(lat,lon,[1 1]);
pix_a = pix_a(:,1);

% time stuff
[xDay,  ~, ~, ~, ~, Md,D]   = createDateVector(time{1},time{2}, 'd');
[xMonth, ~, ~, ~, ~, M]     = createDateVector(time{1},time{2}, 'm');

tick_locations  = xMonth(M==1);
tick_locationsD = xDay(D==1);

nPix   = size(pix_a,1);
nTix   = size(xDay,2);


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

zonesID     = unique(KG_v);
zoneNames   = {'R1- Cold', 'R2- Temperate', 'R3- Humid',  'R4- Sub-humid',   'R5- Semi-arid'};

%% Comparison of wTWS + wTWS+wRiver with GRACE TWS
expNames = fieldnames(mod);

% squeeze wTotal of VEG
mod.MOD.wTotal = squeeze(mod.MOD.wTotal);

for eN=2:numel(expNames)
    eName = expNames{eN};
    
    % squeeze wTotal
    mod.(eName).wTotal = squeeze(mod.(eName).wTotal);
    
    %wRiver units
    mod.(eName).wRiver  = mod.(eName).wRiver ./ pix_a;
    %add wRiver to wTWS
    mod.(eName).wRiverTWS = mod.(eName).wRiver +  mod.(eName).wTotal;
end

% prepare obs


oData   = obs.TWS_GRACE;
if TWScat == 1
    oData(oData<=-500) = NaN;
    oData(oData>=500)  = NaN;
end

% only look at years, where we also have wrIver from CaMa flood 
% xMonth(35)
oData(:,1:34) = NaN;

m       = nanmean(oData,2);
oData   = oData - repmat(m,1,length(xMonth));
% % detrend
% oData = DetrendMatrix(oData);

oData_MSC   = calcMSC(oData, M);

uData = obs.TWS_GRACE_unc;
if TWScat == 1
    uData(oData<=-500) = NaN;
    uData(oData>=500)  = NaN;
end
uData(:,1:34) = NaN;
uData_MSC   = calcMSC(uData, M);

% put into structure
dataObs.wTWS     = oData;
dataObsMSC.wTWS  = oData_MSC;

dataUnc.wTWS     = uData;
dataUncMSC.wTWS  = uData_MSC;

%prepare mod
NaNidx  = find(isnan(dataObs.wTWS));

for eN=1:numel(expNames)
    expN    = expNames{eN};
    if strcmp(expN,'MOD')
        sData_d = mod.(expN).wTotal;   
    else
        sData_d = mod.(expN).wRiverTWS;
    end
    
%     %detrend data
%     sData_d = DetrendMatrix(sData_d);
% 
    % monthly aggregation
    sData = aggDay2Mon(sData_d, info.tem.model.time.sDate, info.tem.model.time.eDate);
    % only consistent data points
    sData(NaNidx) = NaN;
    
    % anomaly if TWS
    m = nanmean(sData,2);
    sData = sData - repmat(m,1,length(xMonth));
    
    % calculate MSC
    sData_MSC   = calcMSC(sData, M);
    
    % put into structure
    data.(expN).wTWS      = sData;
    dataMSC.(expN).wTWS   = sData_MSC;

end


%% combined msc & mef maps
dataNames = ['obs' expNames'];

% GRACE vs wTWS+wRiver
dataXX = {};
for op=1:numel(expNames)
    dataXX = [dataXX  dataMSC.(expNames{op})];
end
obsTmp.wTWS      = dataObsMSC.wTWS;


colXX  = [rgb('Black'); rgb('green'); rgb('Chocolate'); rgb('DarkMagenta');  rgb('Blue'); rgb('Crimson')];%; rgb('DimGray')];
unitNames_tmp     = {'mm'};

tmpExpNames = fieldnames(dataMSC);
tmpObs      = dataObsMSC.wTWS;
for fN=1:numel(tmpExpNames)
    expName = tmpExpNames{fN};
    mefMSC.(expName) = NaN(nPix,1);
    for nP=1:nPix
        mefMSC.(expName)(nP,1) = calcMEF(tmpObs(nP,:),dataMSC.(expName).wTWS(nP,:), ones(size(tmpObs(nP,:))));
    end
end

tName       = 'Gridwise MEF with GRACE-TWS';
Lim.colLim  = [-1 1];
Lim.colLimDiff  = [-0.2 0.2];
Lim.col       = othercolor('RdBu11',100);
Lim.colDiff   = othercolor('PuOr11',201);

[sp] = plot_CombinedMSC_MEFmaps4B_withCaMa(tmpExpNames, obsTmp, dataXX, mefMSC.MOD, mefMSC, lat, lon, zoneNames, KG_v, pix_a, unitNames_tmp, colXX, Lim)
print(gcf,[spth 'Fig5_Effect_Validation_withCaMa.png'],'-dpng','-r300');

    
