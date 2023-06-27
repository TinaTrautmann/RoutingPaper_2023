%% S5 Evaluation of simulated discharge at GRDC stations - with CaMa
% by Tina Trautmann, Aug. 2022 & Oct. 2022


%% set save pth
spth = [pwd '/figures/supplement/'];
if ~exist(spth, 'dir'), mkdir(spth),end

%% download and save the GRDC data as .nc e.g.
try
    pth_lp = 'M:/data/DataStructureMDI/DATA/Incoming/GRDC/2016/Data/'
catch
    error('please provide the folder path of the GRDC discharge measurements.')
end

%% load 0.25 deg spatial data (CaMa)
glob_25     = load([pwd '/data/input/ancillary/CaMa_Flood_up_area_validCali.mat'], 'lat_25', 'lon_25', 'study_mask_25');
glob_25.idx = find(~isnan(glob_25.study_mask_25));

load('data/input/ancillary/CaMa_Flood_up_area_validCali.mat', 'uparea_cama')
glob_25.uparea = uparea_cama;
clear uparea_cama

%% load 1 deg spatial data (TRIP)
glob_pth = [pwd '/data/input/studyArea/globalBaseline_Constraints_1deg.mat']
glob = load(glob_pth, 'lat', 'lon');

rrData     = load('data/input/ancillary/RoutingData_1deg.mat');

[pix_x,pix_y] = LatLon2PixelCoord(glob.lon,glob.lat,90,-180,1);
glob.idx = sub2ind([180, 360],pix_y, pix_x);

glob.flwacc       = rrData.flwacc(glob.idx);

%% which stations?
station_names   = {'Lena', 'Yenisey', 'Ob',   'Danube',  'Mississippi',  'Amazonas', 'Congo', 'Zambesi'}
station_numbers = {2903430, 2909150, 2912600,   6742900,  4127800,  3629001, 1147010, 1291100}


T_all  = readtable('data/input/ancillary/GRDC_stations.xlsx', 'ReadRowNames', true);

for nS=1:numel(station_numbers)
    S_name = station_names{nS};
    S_num  = station_numbers{nS};
    f_num  = num2str(S_num);
    f_in   = [pth_lp  f_num(1) '/' num2str(S_num) '.nc'];
    
    
    Q_d = ncread(f_in, 'CALC_day');
    dtime  = ncread(f_in, 'dtime');

    
    
    Q_lat = T_all{f_num,'lat'}
    Q_lon = T_all{f_num,'long'}
    
    ST_area   = T_all{f_num,'area'}; %catchment size in km²
    ST_name = T_all{f_num,'station'};
    if strcmp(S_name, 'Zambesi')
        ST_name = strsplit(ST_name{1},'(');
        ST_name = ST_name{1};
    end
    

    % the closest lat/lon for 0.25deg
    lat_diff = glob_25.lat_25-Q_lat;
    lon_diff = glob_25.lon_25-Q_lon;
    dist = sqrt((lon_diff).^2 + (lat_diff).^2);
    dist(dist == 0) = Inf;
    [~,closest_id]  = min(dist);
    Q_lat_25 =  glob_25.lat_25(closest_id);
    Q_lon_25 =  glob_25.lon_25(closest_id);

    
    % the closest lat/lon for 1deg
    lat_diff = glob.lat-Q_lat;
    lon_diff = glob.lon-Q_lon;
    dist = sqrt((lon_diff).^2 + (lat_diff).^2);
    dist(dist == 0) = Inf;
    [~,closest_id]  = min(dist);
    Q_lat_1 =  glob.lat(closest_id);
    Q_lon_1 =  glob.lon(closest_id);
    
    % the coordinates
    tmp_lat = fix(Q_lat);
    if sign(tmp_lat)==1
        Q_lat = tmp_lat+0.5;
    else
        Q_lat = tmp_lat-0.5;
    end
    
    tmp_lon = fix(Q_lon);
    if sign(tmp_lon)==1
        Q_lon = tmp_lon+0.5;
    else
        Q_lon = tmp_lon-0.5;
    end
    
    
    Q_datetime        = datetime(1582,10,15) + days(dtime);
    Q_datetime.Format = 'yyyy-MM-dd';
    Q_time            = {datestr(Q_datetime(1),'yyyy-mm-dd') datestr(Q_datetime(end),'yyyy-mm-dd')}
    
    Q_d(Q_d==-9999) = NaN;
    
    T_GRDC.(S_name).Q_d  = Q_d;
    T_GRDC.(S_name).time = Q_time;
    T_GRDC.(S_name).lat  = Q_lat;
    T_GRDC.(S_name).lon  = Q_lon;
    T_GRDC.(S_name).lat_1   = Q_lat_1;
    T_GRDC.(S_name).lon_1   = Q_lon_1;
    T_GRDC.(S_name).lat_25  = Q_lat_25;
    T_GRDC.(S_name).lon_25  = Q_lon_25;
    T_GRDC.(S_name).datetime = Q_datetime;
    T_GRDC.(S_name).stationName = ST_name;
    T_GRDC.(S_name).area    = ST_area;
    T_GRDC.(S_name).area_m  = ST_area .* 1000000;
end

% manually change some coordinates to fit the river sequence in TRIP
T_GRDC.Mississippi.lon  = -91.5;
T_GRDC.Amazonas.lat     = -2.5;
    
% manually change some coordinates to fit the upstream area in CaMa
T_GRDC.Danube.lat_25    = 45.375; % T_GRDC.Danube.lat_25 + 0.25
T_GRDC.Lena.lon_25      = 127.1250; % T_GRDC.Lena.lon_25 +  0.25

% test_lat  = T_GRDC.Lena.lat_25;
% test_lon0 = T_GRDC.Lena.lon_25;
% test_lon1 = T_GRDC.Lena.lon_25 +  0.25;
% test_lon2 = T_GRDC.Lena.lon_25 -  0.25;
% 
% test_idx0 = find(glob_25.lat_25 == test_lat & glob_25.lon_25 == test_lon0)
% test_idx1 = find(glob_25.lat_25 == test_lat & glob_25.lon_25 == test_lon1)
% test_idx2 = find(glob_25.lat_25 == test_lat & glob_25.lon_25 == test_lon2)
% 
% test = [glob_25.uparea(test_idx0), glob_25.uparea(test_idx1), glob_25.uparea(test_idx2)]
% [~,test_max] = max(test) % -> test_idx1 = test_lon1
% % (not test for the max but for the most similar upstream area)

% check for highest upstream area


%% Plot CaMa with Stations 
data    = glob_24.uparea;
land    = shaperead('landareas', 'UseGeoCoords', true);
rivers  = shaperead('worldrivers', 'UseGeoCoords', true);

geoRaRef    = georasterref('RasterSize', [720 1440], 'RasterInterpretation', 'cells',  ...
    'LatitudeLimits', [-90 90], 'LongitudeLimits', [-180 180]);
Z           = NaN(720, 1440); % prepare mapgrid
Z           = imbedm(glob_25.lat_25, glob_25.lon_25, data, Z, geoRaRef); % insert data in mapgrid
Z           = [Z ones(720,1)];  % add column for correct plotting of last column
Z           = [Z; ones(1,1441)]; % add row for correct plotting of last row

figure
ax=worldmap('World');
set(ax, 'layer', 'top','xcolor', 'w', 'ycolor', 'w'); %'Position', [0.1 0.12 0.8 0.75]
axis off; grid on; axis tight;
setm(gca,'ParallelLabel', 'off','MeridianLabel','off', 'FLineWidth', 0.5);
geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
hold on
sm = surfm([-90 90], [-180 180], Z);
geoshow(ax, rivers, 'Color', [.25 .25 .25]);
caxis([0 prctile(data,80)]);
colormap(othercolor('PuBu9'));
cb  =   colorbar;
set(cb,'FontSize',8, 'Position', [0.1 0.09 0.8 0.02],...
    'Orientation', 'horizontal');

for nS=1:numel(station_numbers)
    S_name = station_names{nS};
    plotm(T_GRDC.(S_name).lat_25,  T_GRDC.(S_name).lon_25, 'r*', 'MarkerSize', 10); hold on
end
t   =   title('Upstream Area & GRDC Station Location');

print(gcf,[spth 'S5a_Upstream area CaMa & GRDC Station Location_025.png'],'-dpng','-r300');



%% Plot TRIP flow accumulation with stations
data    = glob.flwacc;
land    = shaperead('landareas', 'UseGeoCoords', true);
rivers  = shaperead('worldrivers', 'UseGeoCoords', true);

geoRaRef    = georasterref('RasterSize', [180 360], 'RasterInterpretation', 'cells',  ...
    'LatitudeLimits', [-90 90], 'LongitudeLimits', [-180 180]);
Z           = NaN(180, 360); % prepare mapgrid
Z           = imbedm(glob.lat, glob.lon, data, Z, geoRaRef); % insert data in mapgrid
Z           = [Z ones(180,1)];  % add column for correct plotting of last column
Z           = [Z; ones(1,361)]; % add row for correct plotting of last row

figure
ax=worldmap('World');
set(ax, 'layer', 'top','xcolor', 'w', 'ycolor', 'w'); %'Position', [0.1 0.12 0.8 0.75]
axis off; grid on; axis tight;
setm(gca,'ParallelLabel', 'off','MeridianLabel','off', 'FLineWidth', 0.5);
geoshow(ax, land, 'FaceColor', [0.85 0.85 0.85]);
hold on
sm = surfm([-90 90], [-180 180], Z);
geoshow(ax, rivers, 'Color', [.25 .25 .25]);
caxis([0 100]);
colormap(othercolor('PuBu9'));
cb  =   colorbar;
set(cb,'FontSize',8, 'Position', [0.1 0.09 0.8 0.02],...
    'Orientation', 'horizontal');

for nS=1:numel(station_numbers)
    S_name = station_names{nS};
    plotm(T_GRDC.(S_name).lat,  T_GRDC.(S_name).lon, 'r*', 'MarkerSize', 10); hold on
end
t   =   title('Flow Accumulation & GRDC Station Location');

print(gcf,[spth 'S5a_Flow Accumulation & GRDC Station Location_2.png'],'-dpng','-r300');

%% load the modelled discharge
expNames             = {'VEG_25_8', 'VEG_05_5', 'VEG_01_7', 'VEG_CaMa_7'};
expNames_short       = {'VEG_25', 'VEG_05', 'VEG_01', 'VEG_CaMa' };
expNames_paper       = {'MOD_R25', 'MOD_R05', 'MOD_R01', 'MOD_CaMa'};
rDates               = {'20210531','20210531','20210531', '20220817'};
space_run            = 'validCali';

varNames = {'Qriver'}

pix_a = AreaGridLatLon(glob.lat,glob.lon,[1 1]);
pix_a = pix_a(:,1);

% load info, Qriver  
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
        if strcmp(expName2,'VEG_CaMa') %load the 0.25 version
            fName   =  [pth_in  expName '_' space_run '_' vName '_025deg.nc'];
            fName2  =  [pth_in   space_run '_' expName2 '_' space_run '_' vName '_025deg.nc'];
        else
            fName   =  [pth_in  expName '_' space_run '_' vName '.nc'];
            fName2  =  [pth_in   space_run '_' expName2 '_' space_run '_' vName '.nc'];
        end
        
        if isfile(fName)
            mod.(expName3).(vName) =    ncread(fName,vName);
        elseif isfile(fName2)
            mod.(expName3).(vName) =    ncread(fName2,vName);
        end
        

    end
end


load(glob_pth, 'time' )

% Plot each time series seperate
for sn=1:numel(station_names)
    stationName = station_names{sn};
    
    time_grdc = T_GRDC.(stationName).time
    [consTime_grdc, consTime_run, time_c, xMonth_c] = calcConsisTime(time_grdc,time, 'daily');
    
    % get the station index in the model data
    statIdx         =   find(glob.lat == T_GRDC.(stationName).lat & glob.lon ==T_GRDC.(stationName).lon);   %idx in nPix with the lat & lon of the GRDC station
    statIdx_025     =   find(glob_25.lat_25 == T_GRDC.(stationName).lat_25 & glob_25.lon_25 ==T_GRDC.(stationName).lon_25);     %idx in nPix with the lat & lon of the GRDC station
    
    t_steps = 'd';
    obsname = ['GRDC-' stationName];
    dataObs = T_GRDC.(stationName).Q_d(consTime_grdc)';
    % unit conversion m3/s into m/s
    dataObs = dataObs ./ T_GRDC.(stationName).area_m;
    % unit conversion m/s into mm/d
    dataObs = dataObs .* 1000 .* 86400;
    
    unit    = 'mm/d';
    col     = othercolor('GnBu9',250);
    colLim  = [];
    
    % comparing discharge of the experiments
    
    string_tmp = [];
    for op = 1:numel(expNames_paper)
        optN    = expNames_paper{op};
        if strcmp(optN,'MOD_CaMa') % data & unit conversion for CaMa
            % unit conversion m3/s into m/s
            dataD   = mod.(optN).Qriver(statIdx_025,consTime_run)./ glob_25.uparea(statIdx_025);% 
            % unit conversion m/s into mm/d
            dataD   = dataD .* 1000 .* 86400;
        else
            dataD   = mod.(optN).Qriver(statIdx,consTime_run)./ T_GRDC.(stationName).area_m;
        end
        eval(char(['data'  num2str(op) ' = dataD;']))
        
        tmp = ['data' num2str(op) ];
        string_tmp = [string_tmp  ',' tmp];
        
    end
    
    leg = ['obs', expNames_paper];
    % comparing observed and modelled discharge of all experiments
    name = ['Discharge - GRDC station ' stationName ' (lat ' num2str(T_GRDC.(stationName).lat) ' / lon ' num2str(T_GRDC.(stationName).lon) ')']
    [sname] = eval(char(['plotTimeSeries2(name, ' char(39) t_steps char(39) ', time_c, leg, dataObs ' string_tmp ');']))
    ylabel(unit)
%     print(gcf,[spth 'S5b_TimeSeries_' sname '.png'],'-dpng','-r300');
    
    
end


%% Plot MSC for all stations in subplots
%load the Qriver of MOD-Rs

colXX  = [rgb('Chocolate'); rgb('DarkMagenta');  rgb('Blue'); rgb('Crimson')];

nrows = 2;
ncols = numel(station_names) / nrows;

figure
set(gcf, 'Position', [5 1 ncols*6 nrows*6]);
ha  = tight_subplot(nrows,ncols,[.1 .05],[.1 .05],[.05 .05]);

for sn=1:numel(station_names)
    axes(ha(sn));
    stationName = station_names{sn};
    
    time_grdc = T_GRDC.(stationName).time
    [consTime_grdc, consTime_run, time_c, xMonth_c] = calcConsisTime(time_grdc,time, 'daily');
    [~,M] = datevec(xMonth_c);
    
    statIdx         =   find(glob.lat == T_GRDC.(stationName).lat & glob.lon ==T_GRDC.(stationName).lon);%idx in nPix with the lat & lon of the GRDC station
    statIdx_025     =   find(glob_25.lat_25 == T_GRDC.(stationName).lat_25 & glob_25.lon_25 ==T_GRDC.(stationName).lon_25);     %idx in nPix with the lat & lon of the GRDC station
    
    t_steps = 'd';
    obsname = ['GRDC-' stationName];
    dataObs = T_GRDC.(stationName).Q_d(consTime_grdc)';
    % unit conversion m3/s into m/s
    dataObs = dataObs ./ T_GRDC.(stationName).area_m;
    % unit conversion m/s into mm/d
    dataObs = dataObs .* 1000 .* 86400;

    unit    = 'Q_D_i_s [mm d^-1]';
    
    dataObs_MSC = calcMSC(dataObs,M);
    plot(1:12,dataObs_MSC, ':', 'Color', rgb('Black'), 'Linewidth', 1.75), hold on
    
    % comparing discharge of the experiments
    for op = 1:numel(expNames_paper)
        optN    = expNames_paper{op};
        
        if strcmp(optN,'MOD_CaMa') % data & unit conversion for CaMa
            % unit conversion m3/s into m/s
            dataD   = mod.(optN).Qriver(statIdx_025,consTime_run)./ glob_25.uparea(statIdx_025);%
            % unit conversion m/s into mm/d
            dataD   = dataD .* 1000 .* 86400;
        else
            dataD   = mod.(optN).Qriver(statIdx,consTime_run);%./ pix_a(statIdx); % in the Qrout_TRIP.m Qin is multiplied by the grid cells area;
            dataD   = dataD ./ T_GRDC.(stationName).area_m; % in the Qrout_TRIP.m Qin is multiplied by the grid cells area;
        end       
%         dataD   = mod.(optN).Qriver(statIdx,consTime_run)./ 86400;% .* info.tem.helpers.areaPix(statIdx);
        dataD_MSC = calcMSC(dataD,M);
        
        plot(1:12,dataD_MSC, '-', 'Color', colXX(op,:), 'Linewidth', 1.25); hold on
        tmpCorr = corr(dataObs_MSC(:),dataD_MSC(:));
        tmpRMSE = calcRMSE(dataObs_MSC(:),dataD_MSC(:));
        tmpPBIAS = mean(dataD_MSC(:)-dataObs_MSC(:))*100./mean(dataObs_MSC(:));

        t1 = text(0.03,1-0.06*op, [ 'corr = ' num2str(round(tmpCorr,2)) ' | PBIAS = ' num2str(round(tmpPBIAS,0)) '%'],'Units','Normalized', 'Fontsize', 8, 'Color', colXX(op,:), 'HorizontalAlignment', 'left');
    end
    if sn==1 || sn==ncols+1
        ylabel(unit, 'FontWeight', 'bold')
    end
    ha(sn).YLim(2) = ha(sn).YLim(2)+ha(sn).YLim(2)/4;
    ha(sn).XLim(1) = 1;
    set(gca, 'FontSize', 8)
    title([stationName  ' (' char(T_GRDC.(stationName).stationName) ')'], 'Fontsize', 10)
    
end

leg = ['obs', cellfun(@(c)[c '_b_e_s_t'],strrep(expNames_paper,'_','-'), 'uni', false)];
l = legend(leg, 'Position', [0.1 0.01 0.8 0.05], 'Orientation', 'Horizontal', 'FontSize', 10, 'Box', 'off')
print(gcf,[spth 'S6_Comparison_MSC_Qdis_GRDCStations_withCaMa.png'],'-dpng','-r300');


