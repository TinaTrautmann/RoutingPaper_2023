%% SELECTING CALIBRATION GRIDS REL RIVER
% by Tina Trautmann, Feb. 2022

%% set save pth
spth = [pwd '/figures/supplement/'];
if ~exist(spth, 'dir'), mkdir(spth),end

% use the wRiver from routed GRUN with eff_vel = 0.5
globData = load([pwd '/data/input/ancillary/GRUNrouted_wRiver_05.mat']) %wRiver, Qriver, lat, lon, xDays

% calculate average wRiver
wRiver_avg = mean(globData.wRiver,2);

figure, histogram(wRiver_avg)
set(gca, 'XLim', [0 200])
prctile(wRiver_avg, [50 70 75 80 85 90 95 99])
mean(wRiver_avg, 'omitnan')


test_lim = prctile(wRiver_avg, 50);
data     = wRiver_avg;
data(data<test_lim) = NaN;
data(data>=test_lim)= 1;

figure,
PlotMapGlobal_noInfo(['Avg wRiver (eff-vel=0.5) | > ' num2str(test_lim) ' | count = ' num2str(sum(data(:), 'omitnan'))],[ ],'', globData.lat, globData.lon, data,[],rgb('DarkBlue'),[],[],[])
print(gcf,[spth 'globalMap_avgRiver_gtMedian.png'],'-dpng','-r300');

dataRiver = data;

%% only consider the validCali grid cells
%load lat lon for validCali

load('input/opti_relRiv/globalBaseline_Constraints_1deg.mat', 'lat', 'lon')

%  get idx in dataRiver (=all) that are part of validCali (=sub)
[all_x,all_y]  = LatLon2PixelCoord(globData.lon,globData.lat,90,-180,1);
idx_all        = sub2ind([180, 360],all_y, all_x);

[sub_x,sub_y]   = LatLon2PixelCoord(lon,lat,90,-180,1);
idx_sub         = sub2ind([180, 360],sub_y, sub_x);

vec_idx_all  = find(ismember(idx_all, idx_sub));
vec_idx_sub  = NaN(length(idx_sub),1);
for ii=1:length(idx_sub)
    vec_idx_sub(ii)  = find(idx_sub == idx_all(vec_idx_all(ii)));
end
%
dataRiver_valid = dataRiver(vec_idx_all,:);

figure,
PlotMapGlobal_noInfo(['Avg wRiver (eff-vel=0.5) | > ' num2str(test_lim) ' | count = ' num2str(sum(dataRiver_valid(:), 'omitnan')) ' / ' num2str(length(lat))],[ ],'', lat, lon, dataRiver_valid,[],rgb('DarkBlue'),[],[],[])
print(gcf,[spth 'validMap_avgRiver_gtMedian.png'],'-dpng','-r300');

%% calculate area of all valid grids
area_v      = AreaGridLatLon(lat,lon,[1 1]);
area_v      = area_v(:,1);
areaTotal_v = sum(area_v(:));
pix         = length(lat);

% -> only grids with dataRiver_valid ==1 are valid now
idx_valids = find(dataRiver_valid == 1);

%% KG zones
load('data\input\ancillary\KoeppenGeiger_1deg.mat','KG');
KG_map   = KG.Agg.KG_new';
KG_v     = KG_map(idx_sub); %valids in map 

% some statistics on global KG
zones = unique(KG_v);

% Plot Percentage # of each class
countPix      = NaN(length(zones),1); 
percentCount  = NaN(length(zones),1); 

for i =1:length(zones)
    id = zones(i);
    countPix(i,1)       = length(find(KG_v==id));
    percentCount(i,1)   = (countPix(i,1).*100)/pix;
end

figure
pie(countPix(:,1))

% Plot Percentage Area of each class
areaKG      = NaN(length(zones),1); %col1=zoneID, col2= area, col3=percentage of area
percentArea = NaN(length(zones),1);

for i =1:length(zones)
    id = zones(i);
    areaKG(i,1)      = sum(area_v(KG_v==id));
    percentArea(i,1) = (areaKG(i,1).*100)./areaTotal_v;
end

figure
pie(areaKG(:,1))

KGnames         = KG.Agg.CodesAgg(ismember(KG.Agg.CodesID,zones))
percentCount    = round(percentCount,2);
percentArea     = round(percentArea,2);


% some statistics on global KG for idx_valid
KG_valid     = KG_v(idx_valids); %valids in vector 
area_valid       = AreaGridLatLon(lat(idx_valids),lon(idx_valids),[1 1]);
area_valid       = area_valid (:,1);
areaTotal_valid  = sum(area_valid (:));
pix_valid        = length(KG_valid);
lat_valid        = lat(idx_valids);
lon_valid        = lon(idx_valids);
zones = unique(KG_valid);

% Plot Percentage # of each class
countPix_valid       = NaN(length(zones),1); 
percentCount_valid   = NaN(length(zones),1); 

for i =1:length(zones)
    id = zones(i);
    countPix_valid (i,1)       = length(find(KG_valid==id));
    percentCount_valid (i,1)   = (countPix_valid (i,1).*100)/pix_valid ;
end

figure
pie(countPix_valid (:,1))

% Plot Percentage Area of each class
areaKG_valid      = NaN(length(zones),1); %col1=zoneID, col2= area, col3=percentage of area
percentArea_valid = NaN(length(zones),1);

for i =1:length(zones)
    id = zones(i);
    areaKG_valid(i,1)      = sum(area_v(KG_valid==id));
    percentArea_valid(i,1) = (areaKG_valid(i,1).*100)./areaTotal_valid;
end

figure
pie(areaKG_valid(:,1))

KGnames         = KG.Agg.CodesAgg(ismember(KG.Agg.CodesID,zones))
percentCount_valid    = round(percentCount_valid,2);
percentArea_valid     = round(percentArea_valid,2);


T = table(zones,countPix,percentCount,areaKG,percentArea,percentCount_valid,percentArea_valid,'RowNames',strrep(KGnames,'_','-'))
writetable(T,[spth '/percentageKG_validCalibration_NEWKG.xls'], 'WriteRowNames', true)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% potential calibration grids % 1st col idx, 2nd col area, 3rd col class; first max of grids
potPix(:,1) = idx_valids;
potPix(:,2) = area_valid;
potPix(:,3) = KG_valid;

%% choose the calibration grids
% 10% of total valid grids 

class_percent = percentArea; %the percentArea of all validCali
class_ID      = KG_valid; 
% area of nCal that should be covered by each KG zone
perc_covered   = 10; %8 percent of area should be covered
areaTotal_cal  = areaTotal_valid/100*perc_covered
areaClass_cal  = round(areaTotal_cal/100 * class_percent)


% calPix = NaN(nCal,2); %1st col: idx in lat_v/lon_v; 2nd col: class ID
allPixArea = area_valid;
allPixIdx  = [1:1:pix_valid]';
calPix     = NaN(pix_valid,3); % 1st col idx, 2nd col area, 3rd col class; first max of grids
cnt    = 1;
% loop over classes and randomly extract the number of grids until their area sums up to the wanted KG area 
class_n     = unique(class_ID);
for cl = 1:length(class_n)
   class_cl             = class_n(cl);
   % extract pixel areas of the current class == the remaining ones from
   % which to extract
   leftArea     = allPixArea(class_ID == class_cl);
   leftIdx      = allPixIdx(class_ID == class_cl);
   % needed area
   neededArea   = areaClass_cal(cl);
   % current area
   isArea       = 0;

   % loop over the remaining subset, randomly extract a grid and sum up the
   % area until it reaches the neededArea
   while isArea < neededArea
       % randomly get a grid
       idx = randi(length(leftArea),1);
       % sum up the area
       isArea = isArea + leftArea(idx);
       
       % write the idx, the KG and the area of that grid within the valid pix vector
       calPix(cnt,1) = leftIdx(idx);
       calPix(cnt,2) = leftArea(idx);
       calPix(cnt,3) = class_cl;
       cnt = cnt+1;
       
       % remove the current grid from the remaining grids
       leftArea(idx) = [];
       leftIdx(idx)  = [];
   end   

end

% remove the NaNs
row_rm = find(isnan(calPix(:,1)));
calPix(row_rm,:) = [];
nCal = size(calPix,1)

sum(calPix(:,2))
areaTotal_cal

% select the calibration pix
lat_c = lat_valid(calPix(:,1));
lon_c = lon_valid(calPix(:,1));

data  = KG_valid;
data  = data(calPix(:,1));

colLabel.size = 5;
colLabel.cticks  = [1:1:8];
colLabel.clabels = strrep(KGnames,'_','-');
figure, PlotMapGlobal2(['Aggregated KG of ' num2str(nCal) ' calibration grids'],['water<50%, permsnow<10%, artificial<10%, no direct human TWS trend, relevant river storage; based on % of agg KG zones'],[],lat_c, lon_c, data, [0.5 8.5],othercolor('Dark28',8),colLabel)
print(gcf,[spth 'aggKG_grids for calibrationAREA' num2str(nCal) '.png'],'-dpng','-r300');

neededArea = areaClass_cal;
isArea     = NaN(numel(zones),1);
for kg = 1:numel(zones)
    tmp_idx     = find(calPix(:,3)==kg);
    tmp         = calPix(tmp_idx,:);
    isArea(kg,1)= sum(tmp(:,2));
end

diffArea = neededArea-isArea;

% Plot Percentage # of each class
isPix      = NaN(length(zones),1); 
perPix     = NaN(length(zones),1); 

for i =1:length(zones)
    id = zones(i);
    isPix(i,1)       = length(find(calPix(:,3)==id));
    perPix(i,1)      = (isPix(i,1).*100)/length(calPix(:,3));
end

figure
pie(isPix(:,1))

T = table(zones,countPix,percentCount,areaKG,percentArea,percentCount_valid,percentArea_valid,neededArea,isArea, diffArea,isPix, perPix,'RowNames',strrep(KGnames,'_','-'))
writetable(T,[spth 'percentageKG_validCalibrationArea' num2str(nCal) '.xls'], 'WriteRowNames', true)
