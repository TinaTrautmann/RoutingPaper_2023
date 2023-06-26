%% Effect of wRiver on model calibration - including CaMa Flood
% by Tina Trautmann, Aug. 2022

%% set save pth
spth = [pwd '/figures/'];

% spth = [pwd '/figures/manuscript/'];
if ~exist(spth, 'dir'), mkdir(spth),end

space_run   = 'validCali';
nloops      = 10;
expNames    = {'VEG','VEG_25','VEG_05','VEG_01','VEG_CaMa'};
expNames3   = {'MOD','MOD_R25','MOD_R05','MOD_R01','MOD_CaMa'};
rDate       = '20210531'

pth_in      = [pwd '/data/model_output/']

varNames    = {'wTotal', 'evapTotal', 'roTotal'};

for ee   = 1:length(expNames)
    expName         = expNames{ee};
    expName3        = expNames3{ee};
    %for n-Loops
    for nL=1:nloops
        expName2 = [expName '_' num2str(nL)];
        expName4 = [expName3 '_' num2str(nL)];
        
        % load variables
        for vn=1:numel(varNames)
            vName = varNames{vn};
            fName = [pth_in expName '/' expName2 '_' space_run '_' vName '.nc'];
            try
                mod.(expName3).(expName4).(vName)          =    ncread(fName,vName);
            catch
                disp(['not a file! : ' fName]);
            end
        end
        
    end
end

% add the orig VEG model from Trautmann et al. 2022, HESS
expName2    = 'E_B_bL_RD4';
% load variables
for vn=1:numel(varNames)
    vName = varNames{vn};
    fName = [pth_in 'VEG/'  space_run '_' expName2 '_' space_run '_' vName '.nc'];
    mod.MOD.(expName2).(vName)   =    ncread(fName,vName);
end

%load info of MOD
load([pth_in 'VEG/' space_run '_' expName2 '_20200629_info.mat']);

% load T_cost
T_cost = readtable([pth_in 'Compare_finalCosts_withCaMa.xls'],'ReadRowNames',true)

% load obs
pth_obs  = [pwd '/data/input/studyArea/globalBaseline_Constraints_1deg.mat'];
info.opti.constraints.oneDataPath  = pth_obs;
obs = readObservation(info);


% space stuff
load(pth_obs, 'lat', 'lon', 'time')
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
sDatas      = {'reshape(mod.eX.XX.wTotal, info.tem.helpers.sizes.nPix, info.tem.helpers.sizes.nTix)',...
    'mod.eX.XX.evapTotal', 'mod.eX.XX.roTotal'};
oDatas      = {'obs.TWS.data', 'obs.Evap.data', 'obs.Q.data'};
oDatas_unc  = {'obs.TWS.unc',  'obs.Evap.unc', 'obs.Q.unc'};
unitNames   = {'mm', 'mm d^-^1', 'mm d^-^1'};

% Loop over variables & put everything into structures 
TWScat = 1; %only consider TWSobs >= -500mm and <= 500mm
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
        optN = expNames3{ee};
        tmp_eX = strrep(sDatas{ii}, 'eX', [optN]);
        dataMSC.(optN).(name).all = NaN(nPix,12,nloops);
        for nL=1:nloops
            expName2 = [optN '_' num2str(nL)];
            tmp = strrep(tmp_eX, 'XX', [expName2]);
            
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
            data.(optN).(expName2).(name)      = sData;
            dataMSC.(optN).(expName2).(name)   = sData_MSC;
            
            % put each n-loop into 3rd dimension
            dataMSC.(optN).(name).all(:,:,nL) = sData_MSC;
        end
    end
end

% add the org VEG model
names       = {'wTWS', 'ET', 'Q'} ;
sDatas      = {'reshape(mod.eX.XX.wTotal, info.tem.helpers.sizes.nPix, info.tem.helpers.sizes.nTix)',...
    'mod.eX.XX.evapTotal', 'mod.eX.XX.roTotal'};

expName2 = 'E_B_bL_RD4';
optN     = 'MOD';
for ii=1:numel(names)
    name  = names{ii};

    tmp_eX = strrep(sDatas{ii}, 'eX', [optN]);
    
    tmp = strrep(tmp_eX, 'XX', [expName2]);
    
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
    data.MOD.(expName2).(name)      = sData;
    dataMSC.MOD.(expName2).(name)   = sData_MSC;
    
    % put each n-loop into 3rd dimension
    dataMSC.MOD.(name).all = cat(3,dataMSC.MOD.(name).all,sData_MSC);
 
end

% get he min, max + mean for each timestep + experiment
%%% Select for each experiment the run with the lowes total costs
expNamesNEW = T_cost.Properties.RowNames;
minCost_Exp = [];
maxCost_Exp = [];
for eN = 1:numel(expNames)
    if eN==1
        [~, idxMin ] = min(T_cost.Total(1:11));
        [~, idxMax ] = max(T_cost.Total(1:11));
        tmp = expNamesNEW(1:11);
    else
        [~, idxMin ] = min(T_cost.Total(2+(eN-1)*10:1+eN*10));
        [~, idxMax ] = max(T_cost.Total(2+(eN-1)*10:1+eN*10));
        tmp = expNamesNEW(2+(eN-1)*10:1+eN*10);
    end
    minCost_Exp = [minCost_Exp tmp(idxMin)]
    maxCost_Exp = [maxCost_Exp tmp(idxMax)]
end

minCost_Exp=strrep(minCost_Exp,'VEG','MOD')
minCost_Exp=[strrep(minCost_Exp(1:4),'MOD_','MOD_R') minCost_Exp(5)]

minCost_Exp = [minCost_Exp(1), minCost_Exp(4), minCost_Exp(3), minCost_Exp(2), minCost_Exp(5)]

% loop over experiments
for ii=1:numel(names)
    name  = names{ii};
    for ee=1:numel(expNames3)
        optN = expNames3{ee};
        dataMSC.(optN).(name).minTotal = dataMSC.(optN).(minCost_Exp{ee}).(name);
    end    
end


%% MSC with 'uncertainty bands' -> plot for each experiment the run with lowest costs 
dataXX = {};
for op=1:numel(expNames3)
    dataXX = [dataXX dataMSC.(expNames3{op})];
end
colXX  = [rgb('Black'); rgb('green'); rgb('Chocolate'); rgb('DarkMagenta'); rgb('Blue'); rgb('Crimson')];


dataObsMSC_tmp.wTWS = dataObsMSC.wTWS;
dataObsMSC_tmp.ET   = dataObsMSC.ET;
dataObsMSC_tmp.Q    = dataObsMSC.Q;
unitNames_tmp       = unitNames;

% dataUncMSC_tmp = [];
% [sname] = plotMSCvarsZones_xExp_bounded(dataNames, dataObsMSC_tmp, dataXX, dataUncMSC_tmp, zoneNames, KG_v, pix_a, unitNames_tmp, colXX);
% print(gcf,[spth 'Fig3a_Effect_Calibration_nLoops_MEFwithObs.png'],'-dpng','-r300');

%metric with MOD instead of obs
dataUncMSC_tmp = [];
[sname] = plotMSCvarsZones_xExp_bounded_MEF1Exp_transposed_CaMa(dataNames, dataObsMSC_tmp, dataXX, dataUncMSC_tmp, zoneNames, KG_v, pix_a, unitNames_tmp, colXX);
print(gcf,[spth 'Fig4_Effect_Calibration_nLoops_MEFwithMOD_CaMa.png'],'-dpng','-r300');
