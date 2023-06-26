%% S3 Analysis of Effect of Selection of Calibration Grid Cells 
% by Tina Trautmann, Feb. 2022


%% ANALYSIS FOR REL RIVER
%% ----- OPTI Parameters
%% set save pth
spth = [pwd '/figures/supplement/'];
if ~exist(spth, 'dir'), mkdir(spth),end

%which experiments?
expNames    = {'VEG', 'VEG_25', 'VEG_05', 'VEG_01'};
rDates      = {'20220110','20220110','20220110','20220110','20220110'};

space_run   = 'validCali';

expNamesNEW = {};

%loop over experiments & load the info
for eN =1:numel(expNames)
    expName = expNames{eN};
    rDate   = rDates{eN};
    
    expNameNEW  = [expName '_relRiv'];
    expNamesNEW = [expNamesNEW expNameNEW]
    
    % path of model output
    pth_in  = [pwd '/data/model_output/relRiv_study/'];
    
    %load info
    infoIn.(expNameNEW) = load([pth_in  expNameNEW '_' rDate '_info.mat']);
end

%% load the original VEG model
% add the orig VEG model from Trautmann et al. 2022, HESS
expName2    = 'E_B_bL_RD4';
%load info of MOD
infoIn.(expName2) = load([pwd '/data/model_output/VEG/' space_run '_' expName2 '_20200629_info.mat']);


%% compare the values (=adding to a table)
%load one optimization info
infoOpti = load([pwd '/data/model_output/relRiv_opti/VEG_relRiv_20220108_info.mat']);

% make the comparison table
p_name = infoOpti.info.opti.params.names
ps     =  infoOpti.info.opti.params;
p_Defvalue  = ps.defaults';
p_lowBound  = ps.lBounds';
p_upBound   = ps.uBounds';
p_name2     = cell(numel(p_name),1);
for pn = 1:numel(p_name)
    paraName    = p_name{pn};
    tmp         = strsplit(paraName,'.');
    p_name2{pn} = tmp{3};
end
T_param    = table(p_Defvalue, p_lowBound, p_upBound, 'RowNames', p_name2);


% add the calibrated values from the validCali forward run
expNamesNEW = [expName2 expNamesNEW ];

for eN =2:numel(expNamesNEW)
    expName = expNamesNEW{eN};
        
    % add the optimized parameter values in the table
    nrow = size(T_param,1);
    T_param.(expName)= NaN(nrow, 1);
    
    for ii=1:numel(p_name)
        paraName    = p_name{ii};
        tmp         = strsplit(paraName,'.');
        % p structure
        T_param.(expName)(  tmp{3} ) = infoIn.(expName).info.tem.params.(tmp{2}).(tmp{3})(1);
    end
end

% write the param comparison table
T_param
writetable(T_param,[spth 'S3_Compare_optimizedParams_relRiv.xls'],'WriteRowNames',true)


%% --------------------------------------------------------------
%% EFFECT ON MODEL CALIBRATION
%% load all valid Cali forward runs
expNames    
expNames3   = {'MOD','MOD_R25','MOD_R05','MOD_R01'}
rDate       = '20220110'
space_run   = 'validCali';

varNames    = {'wTotal', 'evapTotal', 'roTotal'};

    
for ee   = 1:length(expNames)
    expName         = expNames{ee};
    expName3        = expNames3{ee};
    
    expName2 = [expName '_relRiv'];
    expName4 = [expName3 '_relRiv'];
    % load variables
    for vn=1:numel(varNames)
        vName = varNames{vn};
        fName = [pth_in expName2 '_' space_run '_' vName '.nc'];
        try
            mod.(expName4).(vName)          =    ncread(fName,vName);
        catch
            disp(['not a file! : ' fName]);
        end
    end
    
end

% add the orig VEG model
expName2    = 'E_B_bL_RD4';

% load variables
for vn=1:numel(varNames)
    vName = varNames{vn};
    fName = [pwd '/data/model_output/VEG/'  space_run '_' expName2 '_' space_run '_' vName '.nc'];
    mod.MOD.(vName)   =    ncread(fName,vName);
end


%load info
info = infoIn.E_B_bL_RD4.info;

% load T_cost
T_cost = readtable(['data/model_output/' 'Compare_finalCosts_relRiv.xls'],'ReadRowNames',true)

% load obs
info.opti.constraints.oneDataPath    = [pwd '/data/input/studyArea/globalBaseline_Constraints_1deg.mat'];
obs = readObservation(info);


TWScat = 1; %only consider TWSobs >= -500mm and <= 500mm

% space stuff
load('data/input/studyArea/globalBaseline_Constraints_1deg.mat', 'lat', 'lon', 'time')
pix_a = AreaGridLatLon(lat,lon,[1 1]);
pix_a = pix_a(:,1);

% time stuff
[xDay,  ~, ~, ~, ~, Md,D]   = createDateVector(info.tem.model.time.sDate,info.tem.model.time.eDate, 'd');
[xMonth, ~, ~, ~, ~, M]     = createDateVector(info.tem.model.time.sDate,info.tem.model.time.eDate, 'm');

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

% what variables to compare
dataNames    = ['obs', expNames3];
names       = {'wTWS', 'ET', 'Q'} ;
sDatas      = {'reshape(mod.XX.wTotal, info.tem.helpers.sizes.nPix, info.tem.helpers.sizes.nTix)',...
    'mod.XX.evapTotal', 'mod.XX.roTotal'};
oDatas      = {'obs.TWS.data', 'obs.Evap.data', 'obs.Q.data'};
oDatas_unc  = {'obs.TWS.unc',  'obs.Evap.unc', 'obs.Q.unc'};
unitNames   = {'mm', 'mm d^-^1', 'mm d^-^1'};

% Loop over variables & put everything into structures for obs, exp1,
% exp2
for ii=1:numel(names)
    name  = names{ii};
    oData = eval([oDatas{ii}]);
    oData_unc = eval([oDatas_unc{ii}]);
    
    % anomaly if TWS
    if strcmp(names{ii},'wTWS')
        if TWScat == 1
            oData(oData<=-500) = NaN;
            oData(oData>=500)  = NaN;
        end
        m = nanmean(oData,2);
        oData = oData - repmat(m,1,length(xMonth));
    end
    NaNidx  = find(isnan(oData));
    
    % MSC
    oData_MSC   = calcMSC(oData, M);
    oData_MSC_unc = calcMSC(oData_unc,M);
    
    % put into structure
    dataObs.(name)      = oData;
    dataObsMSC.(name)   = oData_MSC;
    dataUnc.(name)      = oData_unc;
    dataUncMSC.(name)   = oData_MSC_unc;
    
    % loop over experiments
    for ee=1:numel(expNames3)
        optN = [expNames3{ee} '_relRiv'];
        tmp  = strrep(sDatas{ii}, 'XX', [optN]);
        dataMSC.(optN).(name) = NaN(nPix,12);
        
        sData_d = eval(tmp);
        % monthly aggregation
        sData = aggDay2Mon(sData_d, info.tem.model.time.sDate, info.tem.model.time.eDate);
        % only consistent data points
        sData(NaNidx) = NaN;
        
        % anomaly if TWS
        if strcmp(names{ii},'wTWS')
            m = nanmean(sData,2);
            sData = sData - repmat(m,1,length(xMonth));
        end
        
        % calculate MSC
        sData_MSC   = calcMSC(sData, M);
        
        % put into structure
        data.(optN).(name)      = sData;
        dataMSC.(optN).(name)   = sData_MSC;
        
    end
end

% add the org VEG model
names       = {'wTWS', 'ET', 'Q'} ;
sDatas      = {'reshape(mod.XX.wTotal, info.tem.helpers.sizes.nPix, info.tem.helpers.sizes.nTix)',...
    'mod.XX.evapTotal', 'mod.XX.roTotal'};

optN     = 'MOD';
for ii=1:numel(names)
    name  = names{ii};    
    tmp = strrep(sDatas{ii}, 'XX', [optN]);
    
    sData_d = eval(tmp);
    % monthly aggregation
    sData = aggDay2Mon(sData_d, info.tem.model.time.sDate, info.tem.model.time.eDate);
    % only consistent data points
    sData(NaNidx) = NaN;
    
    % anomaly if TWS
    if strcmp(names{ii},'wTWS')
        m = nanmean(sData,2);
        sData = sData - repmat(m,1,length(xMonth));
    end
    
    % calculate MSC
    sData_MSC   = calcMSC(sData, M);
    
    % put into structure
    data.(optN).(name)      = sData;
    dataMSC.(optN).(name)   = sData_MSC;
     
end


dataNames    = {'obs', 'MOD', 'MOD_relRiv', 'MOD_R25_relRiv', 'MOD_R05_relRiv','MOD_R01_relRiv'};

%% MSC with 'uncertainty bands' -> plot for each experiment the run with lowest costs (instead of median)
expNames3 = dataNames(2:end)
dataXX = {};
for op=1:numel(expNames3)
    dataXX = [dataXX dataMSC.(expNames3{op})];
end
colXX  = [rgb('Black'); rgb('Red'); rgb('green'); rgb('Chocolate'); rgb('DarkMagenta'); rgb('Blue')];

dataObsMSC_tmp.wTWS = dataObsMSC.wTWS;
dataObsMSC_tmp.ET   = dataObsMSC.ET;
dataObsMSC_tmp.Q    = dataObsMSC.Q;
unitNames_tmp       = unitNames;

% dataUncMSC_tmp = [];
% [sname] = plotMSCvarsZones_xExp(dataNames, dataObsMSC_tmp, dataXX, dataUncMSC_tmp, zoneNames, KG_v, pix_a, unitNames_tmp, colXX);
% print(gcf,[spth sname '.png'],'-dpng','-r300');

%metric with MOD instead of obs
dataUncMSC_tmp = [];
[sname] = plotMSCvarsZones_xExp_bounded_MEF1Exp_transposed_relRiv(dataNames, dataObsMSC_tmp, dataXX, dataUncMSC_tmp, zoneNames, KG_v, pix_a, unitNames_tmp, colXX);
print(gcf,[spth 'S3_RelRiv_Calibration_MEFwithMOD.png'],'-dpng','-r300');

%% ---------- VALIDATION
gone

spth = [pwd '/figures/supplement/'];
if ~exist(spth, 'dir'), mkdir(spth),end

%which experiments?
expNames    = {'VEG', 'VEG_25', 'VEG_05', 'VEG_01'};
expNames2   = {'MOD', 'MOD_R25', 'MOD_R05', 'MOD_R01'};

rDates      = {'20220110','20220110','20220110','20220110','20220110'};

space_run   = 'validCali';
varNames    = {'wTotal', 'wRiver'};

expNamesNEW = {};
%loop over experiments & load the info
for eN =1:numel(expNames)
    expName = expNames{eN};
    expName2 = expNames2{eN};
    rDate   = rDates{eN};
    
    expNameNEW  = [expName '_relRiv'];
    expNamesNEW = [expNamesNEW expNameNEW]
    
    expName4    = [expName2 '_relRiv']
    % path of model output
    pth_in  = [pwd '/data/model_output/relRiv_study/'];
    
    % load variables
    for vn=1:numel(varNames)
        vName = varNames{vn};
        fName = [pth_in expNameNEW '_' space_run '_' vName '.nc'];
        try
            mod.(expName4).(vName)          =    ncread(fName,vName);
        catch
            disp(['not a file! : ' fName]);
        end
    end
    
    %load info
    infoIn.(expNameNEW) = load([pth_in  expNameNEW '_' rDate '_info.mat']);
end

%% load the original VEG model
% add the orig VEG model from Trautmann et al. 2022, HESS
expName2    = 'E_B_bL_RD4';
%load info of MOD
infoIn.MOD = load([pwd '/data/model_output/VEG/' space_run '_' expName2 '_20200629_info.mat']);

% load variables
for vn=1:numel(varNames)
    vName = varNames{vn};
    fName = [pwd '/data/model_output/' space_run '_' expName2 '_' space_run '_20200629/modelOutput/'  space_run '_' expName2 '_' space_run '_' vName '.nc'];
    mod.MOD.(vName)   =    ncread(fName,vName);
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
expNames = [expNames(end) expNames(1:end-1)']

% squeeze wTotal of VEG
mod.MOD.wTotal        = squeeze(mod.MOD.wTotal);
mod.MOD_relRiv.wTotal = squeeze(mod.MOD_relRiv.wTotal);

for eN=3:numel(expNames)
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
m       = nanmean(oData,2);
oData   = oData - repmat(m,1,length(xMonth));

oData_MSC   = calcMSC(oData, M);

uData = obs.TWS_GRACE_unc;
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

%prepare mod
NaNidx  = find(isnan(dataObs.wTWS));

for eN=1:numel(expNames)
    expN    = expNames{eN};
    if strcmp(expN,'MOD')
        sData_d = mod.(expN).wTotal;   
    elseif strcmp(expN,'MOD_relRiv')
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


colXX  = [rgb('Black');  rgb('DarkRed'); rgb('green'); rgb('Chocolate'); rgb('DarkMagenta');  rgb('Blue')];%; rgb('DimGray')];
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

[sp] = plot_CombinedMSC_MEFmaps4C(tmpExpNames, obsTmp, dataXX, mefMSC.MOD, mefMSC, lat, lon, zoneNames, KG_v, pix_a, unitNames_tmp, colXX, Lim)
print(gcf,[spth 'S4_relRiv_Effect_Validation.png'],'-dpng','-r300');

