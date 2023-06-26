%% S3 Analysis of Effect of Selection of Calibration Grid Cells - percentage wRiver with CaMa
% by Tina Trautmann, Aug. 2022

%% set save pth
spth = [pwd '/figures/supplement/'];
if ~exist(spth, 'dir'), mkdir(spth),end

%% 
%% Distribution of Calibration grids vs. all potentiall calibration
% constraints of rel river
calR    = load([pwd '/data/input/opti_relRiv/globalBaseline_Constraints_1deg.mat'], 'lat', 'lon');

% constraints of opti grids
cal     = load([pwd '/data/input/opti/globalBaseline_Constraints_1deg.mat'], 'lat', 'lon');

% Global study area
glob    = load([pwd '/data/input/studyArea/globalBaseline_Constraints_1deg.mat'], 'lat', 'lon')



%% Compare the prctiles of wRiver
data_names = {'wRiver_01', 'wRiver_05', 'wRiver_25', 'wRiver_CaMa'}
pth_in   = [pwd '/data/input/ancillary/']
prt      = [75, 90, 95, 98];

% get idx in vector
[all_x,all_y]  = LatLon2PixelCoord(glob.lon,glob.lat,90,-180,1);
idx_glob       = sub2ind([180, 360],all_y, all_x);

[cal_x,cal_y]  = LatLon2PixelCoord(cal.lon,cal.lat,90,-180,1);
idx_cal        = sub2ind([180, 360],cal_y, cal_x);

[calR_x,calR_y]  = LatLon2PixelCoord(calR.lon,calR.lat,90,-180,1);
idx_calR        = sub2ind([180, 360],calR_y, calR_x);


figure, set(gcf, 'Position', [5 5 10 15])
for vn=1:numel(data_names)
    data_name = data_names{vn} %'TWS_GRACE';
    if strcmp(data_name, 'wRiver_CaMa')
        allData         = load([pth_in 'CaMa_routed_storge.mat'], 'wRiver_cama_sum', 'lat', 'lon');       
        allData.wRiver  = allData.wRiver_cama_sum;
        allData.wRiver_cama_sum = [];
    else
        allData         = load([pth_in 'GRUNrouted_' data_name '.mat'], 'wRiver', 'lat', 'lon');
    end
    [allData_x,allData_y]  = LatLon2PixelCoord(allData.lon,allData.lat,90,-180,1);
    idx_all             = sub2ind([180, 360],allData_y, allData_x);

    vec_idx_glob  = find(ismember(idx_all, idx_glob)); %idx in allData vector that is member of glob
    vec_idx_cal   = find(ismember(idx_all, idx_cal)); %idx in allData vector that is member of cal
    vec_idx_calR  = find(ismember(idx_all, idx_calR)); %idx in allData vector that is member of calR


    glob_data = allData.wRiver(vec_idx_glob,:); 
    cal_data  = allData.wRiver(vec_idx_cal,:);    
    calR_data = allData.wRiver(vec_idx_calR,:); 

    
    name     = strrep(data_name,'_','-');

    dataTmp  = NaN(3,length(prt)); %dataTmp = (m,n); m=groups=glob,cal,calR ; n=bars=prctiles
    
    for pn=1:length(prt)
        dataTmp(1,pn) = prctile(glob_data(:), prt(pn));
        dataTmp(2,pn) = prctile(cal_data(:), prt(pn));
        dataTmp(3,pn) = prctile(calR_data(:), prt(pn));
    end
    
    % create subplot
    subplot(numel(data_names),1,vn)
    b=bar(dataTmp', 'FaceColor', 'flat')
    set(gca, 'XTickLabel', prt)
    title(strrep(data_name,'_','-'))
    b(1).CData = rgb('DarkBlue');
    b(2).CData = rgb('Salmon');
    b(3).CData = rgb('GoldenRod');
    ylabel('daily value [mm]', 'FontSize',8);
end
xlabel('percentiles of all  data points in space and time', 'FontSize', 8)
legend({'global', 'cal', 'cal-relRiv'},  'Orientation', 'Horizontal', 'Box', 'off', 'Fontsize', 8, 'Position', [0.15 0.01 0.85 0.03]);
an = annotation('textbox',[0.05 .96 1 .04],'String', 'Percentiles of wRiver for different spatial subsets' ,'LineStyle','none', 'Fontsize', 11, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top');
print(gcf, [spth 'S3_wRiver_prctiles_global_vs_opti_CaMa.png'], '-dpng', '-r300')
    
%%

