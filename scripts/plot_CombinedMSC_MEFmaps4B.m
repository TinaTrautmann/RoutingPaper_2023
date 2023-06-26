function [sp] = plot_CombinedMSC_MEFmaps4B(dataNames, dataObs, dataXX, mefObs, mefMod, lat, lon, zoneNames, zoneIdx, pix_a, unitNames, colEXP, Lim)

% big figure with
% A) difference maps of gridwise MEF for different regions 
% +  boxplots of gridwise MEF difference of different experiments globbally + for
% different zones
% B) MSC for global + zones

if isempty(colEXP)
    colEXP         = [rgb('DarkRed'); rgb('Blue'); rgb('Green')];
end

% should only be one varNames
try
    varNames  = fieldnames(dataObs);
catch
    varNames = fieldnames(dataXX{1});
end


expNames = strrep(dataNames,'_','-');
zoneNames = strrep(zoneNames,'_','-');

eN = numel(dataXX);


zNames    = ['Global', zoneNames];
zNames2   = strrep(zNames,'-','');
zNames2   = strrep(zNames2,' ','');

nrows   = 6;
ncols   = 3;
xSeason = 1:1:12;

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

% preps maps
% preps
colLim      = Lim.colLim;
colLimDiff  = Lim.colLimDiff;

col         = Lim.col;      %othercolor('RdBu4',100);
colDiff     = Lim.colDiff;  %othercolor('PuOr11',201);

land        = shaperead('landareas', 'UseGeoCoords', true);
rivers      = shaperead('worldrivers', 'UseGeoCoords', true);
geoRaRef    = georasterref('RasterSize', [180 360], 'RasterInterpretation', 'cells',  ...
    'LatitudeLimits', [-90 90], 'LongitudeLimits', [-180 180]);
Z1           = NaN(180, 360); % prepare mapgrid

% for boxplot dataTmp = [zones, exp, values]
dataTmp = NaN(numel(zoneNames)+1,numel(expNames),length(pix_a));


%% FIGURE
sp  = figure;
% set(gcf, 'Position', [5 1 ncols*7 nrows*6]);
set(gcf, 'Position', [5 0 16 27]);

ha  = tight_subplot(nrows,ncols,[.05 .02],[.08 .08],[.05 .02]);

cnt = 1;

% A) MAPS
rr=1;

a_tmp = (rr-1)*ncols;
delete(ha(a_tmp+length(dataNames)./2+1:a_tmp+ncols))

%map of obs
cc=1

data        = mefObs;
Z           = NaN(180, 360); % prepare mapgrid
Z           = imbedm(lat, lon, data, Z, geoRaRef); % insert data in mapgrid
Z           = [Z ones(180,1)];  % add column for correct plotting of last column
Z           = [Z; ones(1,361)]; % add row for correct plotting of last row

axes(ha(a_tmp+cc));
worldmap('World');
ha(a_tmp+cc).Layer = 'top';
setm(gca,'ParallelLabel', 'off','MeridianLabel','off', 'FLineWidth', 0.5);
geoshow(ha(a_tmp+cc), land, 'FaceColor', [0.85 0.85 0.85])
hold on
sm = surfm([-90 90], [-180 180], Z);
geoshow(ha(a_tmp+cc), rivers, 'Color', [.25 .25 .25]);

colormap(ha(a_tmp+cc),col);
caxis(ha(a_tmp+cc),colLim);

t   =   title([ expNames{cc} '_b_e_s_t'],'Fontsize', 10);

dataTmp(1,1,:)=data;
for zN=1:numel(zoneNames)
    dataTmp(1+zN,1,1:length(find(zoneIdx==zN)))=data(zoneIdx==zN);
end

ha(a_tmp+cc).Position= [ha(a_tmp+cc).Position(1) 0.5 0.45 0.7]

cb  =   colorbar;
set(cb,'FontSize',8, 'Position', [0.07   0.96   0.4    0.005], 'Orientation', 'horizontal', 'TickLabels', {'-1', '', '', '', '1'}, 'AxisLocation', 'in');
cb1 =   ylabel(cb, 'MEF', 'FontSize', 9);
cb1.Position(2) = 4;

% MOD 2
cc = 2

data        = mefObs-mefMod.(dataNames{2});
Z           = NaN(180, 360); % prepare mapgrid
Z           = imbedm(lat, lon, data, Z, geoRaRef); % insert data in mapgrid
Z           = [Z ones(180,1)];  % add column for correct plotting of last column
Z           = [Z; ones(1,361)]; % add row for correct plotting of last row

axes(ha(a_tmp+cc));
worldmap('World');
ha(a_tmp+cc).Layer = 'top';
setm(gca,'ParallelLabel', 'off','MeridianLabel','off', 'FLineWidth', 0.5);
geoshow(ha(a_tmp+cc), land, 'FaceColor', [0.85 0.85 0.85])
hold on
sm = surfm([-90 90], [-180 180], Z);
geoshow(ha(a_tmp+cc), rivers, 'Color', [.25 .25 .25]);

colormap(ha(a_tmp+cc),colDiff);
caxis(ha(a_tmp+cc),colLimDiff);

title([expNames{2} '_b_e_s_t'],'Fontsize', 10)

tmpExp1 = squeeze(dataTmp(1,1,:));
dataTmp(1,2,:)= tmpExp1 - mefMod.(dataNames{2}); % MEF of 1 - MEF of this
for zN=1:numel(zoneNames)
    tmpExp1 = squeeze(dataTmp(zN+1,1,:));
    dataTmp(1+zN,2,1:length(find(zoneIdx==zN)))=tmpExp1(1:length(find(zoneIdx==zN)),:) - mefMod.(dataNames{2})(zoneIdx==zN);
end

ha(a_tmp+cc).Position= [ha(a_tmp+cc).Position(1)+0.15 0.5 0.45 0.7]

%colorbar for difference
cb  =   colorbar;
set(cb,'FontSize',8, 'Position', [0.55   0.96   0.4    0.005], 'Orientation', 'horizontal', 'TickLabels', {'-0.2', '', '', '', '0.2'}, 'AxisLocation', 'in');
cb1 =   ylabel(cb, 'MEF difference', 'FontSize', 9);
cb1.Position(2) = 4;

% mod 2+3 
rr = 2
a_tmp = (rr-1)*ncols;
delete(ha(6)) %!!!

cc = 1      
data        = mefObs-mefMod.(dataNames{3});
Z           = NaN(180, 360); % prepare mapgrid
Z           = imbedm(lat, lon, data, Z, geoRaRef); % insert data in mapgrid
Z           = [Z ones(180,1)];  % add column for correct plotting of last column
Z           = [Z; ones(1,361)]; % add row for correct plotting of last row

axes(ha(a_tmp+cc));
worldmap('World');
ha(a_tmp+cc).Layer = 'top';
setm(gca,'ParallelLabel', 'off','MeridianLabel','off', 'FLineWidth', 0.5);
geoshow(ha(a_tmp+cc), land, 'FaceColor', [0.85 0.85 0.85])
hold on
sm = surfm([-90 90], [-180 180], Z);
geoshow(ha(a_tmp+cc), rivers, 'Color', [.25 .25 .25]);

colormap(ha(a_tmp+cc),colDiff);
caxis(ha(a_tmp+cc),colLimDiff);

title([expNames{3} '_b_e_s_t'],'Fontsize', 10)

tmpExp1 = squeeze(dataTmp(1,1,:));
dataTmp(1,3,:)= tmpExp1 - mefMod.(dataNames{3}); % MEF of 1 - MEF of this
for zN=1:numel(zoneNames)
    tmpExp1 = squeeze(dataTmp(zN+1,1,:));
    dataTmp(1+zN,3,1:length(find(zoneIdx==zN)))=tmpExp1(1:length(find(zoneIdx==zN)),:) - mefMod.(dataNames{3})(zoneIdx==zN);
end
        
ha(a_tmp+cc).Position= [ha(a_tmp+cc).Position(1) 0.27 0.45 0.7]

cc = 2
data        = mefObs-mefMod.(dataNames{4});
Z           = NaN(180, 360); % prepare mapgrid
Z           = imbedm(lat, lon, data, Z, geoRaRef); % insert data in mapgrid
Z           = [Z ones(180,1)];  % add column for correct plotting of last column
Z           = [Z; ones(1,361)]; % add row for correct plotting of last row

axes(ha(a_tmp+cc));
worldmap('World');
ha(a_tmp+cc).Layer = 'top';
setm(gca,'ParallelLabel', 'off','MeridianLabel','off', 'FLineWidth', 0.5);
geoshow(ha(a_tmp+cc), land, 'FaceColor', [0.85 0.85 0.85])
hold on
sm = surfm([-90 90], [-180 180], Z);
geoshow(ha(a_tmp+cc), rivers, 'Color', [.25 .25 .25]);

colormap(ha(a_tmp+cc),colDiff);
caxis(ha(a_tmp+cc),colLimDiff);

title([expNames{4} '_b_e_s_t'],'Fontsize', 10)

tmpExp1 = squeeze(dataTmp(1,1,:));
dataTmp(1,4,:)= tmpExp1 - mefMod.(dataNames{4}); % MEF of 1 - MEF of this
for zN=1:numel(zoneNames)
    tmpExp1 = squeeze(dataTmp(zN+1,1,:));
    dataTmp(1+zN,4,1:length(find(zoneIdx==zN)))=tmpExp1(1:length(find(zoneIdx==zN)),:) - mefMod.(dataNames{4})(zoneIdx==zN);
end
ha(a_tmp+cc).Position= [ha(a_tmp+cc).Position(1)+0.15 0.27 0.45 0.7]



% B) boxplots of gridwise MEF of different experiments globaly + for
% different zones
rr=4;
% plot the boxplots in row below
a_tmp = 6;%(rr-1)*ncols;
delete(ha(9))

zoneNames2 = {'R1', 'R2', 'R3', 'R4', 'R5'}

cc=1
axes(ha(a_tmp+cc));
xTmp = 1:length(zoneNames)+1;
p = boxplot2(dataTmp(:,1,:),xTmp, 'barwidth', 0.5);
ha(a_tmp+cc).XTickLabel = ['Global', zoneNames2];
ha(a_tmp+cc).XTick     = [1:1:length(zoneNames)+1];
ha(a_tmp+cc).XLim = [0.5 length(zoneNames)+1.5];
ha(a_tmp+cc).FontSize = 8;
ha(a_tmp+cc).Box     = 'on';
ha(a_tmp+cc).YTickLabelMode = 'auto';
structfun(@(x) set(x(1,:), 'color', colEXP(2,:), ...
    'markeredgecolor', colEXP(2,:)), p);
set([p.box],'linewidth', 1.5);
set([p.med],'linewidth', 1.5);
set([p.lwhis p.uwhis], 'linestyle', '-', 'linewidth', 1.5);
set(p.out, 'marker', 'none');
hold on,
rf=refline(0,0);
rf.Color =  rgb('Black');
ha(a_tmp+cc).YLim = colLim;

ha(a_tmp+cc).Position = [0.12+(cc-1)*0.23 0.75 .3 0.05]; %+cc*
% ax2=axes('Position', get(ha(a_tmp+cc),'outerpos'), 'XColor', rgb('White'), 'XTick', [], 'YTick', [], 'YColor', rgb('White'))
% uistack(ax2,'down',1)

cc=2
axes(ha(a_tmp+cc));

xTmp = 1:length(zoneNames)+1;
p = boxplot2(dataTmp(:,2,:),xTmp, 'barwidth', 0.5);
ha(a_tmp+cc).XTickLabel = ['Global', zoneNames2];
ha(a_tmp+cc).XTick     = [1:1:length(zoneNames)+1];
ha(a_tmp+cc).XLim = [0.5 length(zoneNames)+1.5];
ha(a_tmp+cc).FontSize = 8;
ha(a_tmp+cc).Box     = 'on';
ha(a_tmp+cc).YTickLabelMode = 'auto';
structfun(@(x) set(x(1,:), 'color', colEXP(3,:), ...
    'markeredgecolor', colEXP(3,:)), p);
set([p.box],'linewidth', 1.5);
set([p.med],'linewidth', 1.5);
set([p.lwhis p.uwhis], 'linestyle', '-', 'linewidth', 1.5);
set(p.out, 'marker', 'none');
hold on,
rf=refline(0,0);
rf.Color =  rgb('Black');
ha(a_tmp+cc).YLim = colLim;
tmpTest = dataTmp(end,2,:);
ha(a_tmp+cc).YLim = [prctile(tmpTest(:),15) prctile(tmpTest(:),85)]%colLim;
% ha(a_tmp+cc).YLim = [min(prctile(tmpTest(:),10),-0.05) max(prctile(tmpTest(:),90),0.05)];
ha(a_tmp+cc).Position = [0.12+(cc-1)*0.47 0.75 .3 0.05];
        
%         ax2=axes('Position', get(ha(a_tmp+cc),'outerpos'), 'XColor', rgb('White'), 'XTick', [], 'YTick', [], 'YColor', rgb('White'))
%         uistack(ax2,'down',1)

rr=5;
% plot the boxplots in row below
a_tmp = 9;%(rr-1)*ncols;
delete(ha(12))

cc=1
axes(ha(a_tmp+cc));
xTmp = 1:length(zoneNames)+1;
p = boxplot2(dataTmp(:,3,:),xTmp, 'barwidth', 0.5);
ha(a_tmp+cc).XTickLabel = ['Global', zoneNames2];
ha(a_tmp+cc).XTick     = [1:1:length(zoneNames)+1];
ha(a_tmp+cc).XLim = [0.5 length(zoneNames)+1.5];
ha(a_tmp+cc).FontSize = 8;
ha(a_tmp+cc).Box     = 'on';
ha(a_tmp+cc).YTickLabelMode = 'auto';
structfun(@(x) set(x(1,:), 'color', colEXP(4,:), ...
    'markeredgecolor', colEXP(4,:)), p);
set([p.box],'linewidth', 1.5);
set([p.med],'linewidth', 1.5);
set([p.lwhis p.uwhis], 'linestyle', '-', 'linewidth', 1.5);
set(p.out, 'marker', 'none');
hold on,
rf=refline(0,0);
rf.Color =  rgb('Black');
tmpTest = dataTmp(end,3,:);
ha(a_tmp+cc).YLim = [min(prctile(tmpTest(:),15),-0.05) prctile(tmpTest(:),85)]%colLim;
% ha(a_tmp+cc).YLim = [min(prctile(tmpTest(:),10),-0.05) max(prctile(tmpTest(:),90),0.05)];
ha(a_tmp+cc).Position = [0.12+(cc-1)*0.23 0.52 .3 0.05]; %+cc*

%         ax2=axes('Position', get(ha(a_tmp+cc),'outerpos'), 'XColor', rgb('White'), 'XTick', [], 'YTick', [], 'YColor', rgb('White'))
%         uistack(ax2,'down',1)

cc=2
axes(ha(a_tmp+cc));

xTmp = 1:length(zoneNames)+1;
p = boxplot2(dataTmp(:,4,:),xTmp, 'barwidth', 0.5);
ha(a_tmp+cc).XTickLabel = ['Global', zoneNames2];
ha(a_tmp+cc).XTick     = [1:1:length(zoneNames)+1];
ha(a_tmp+cc).XLim = [0.5 length(zoneNames)+1.5];
ha(a_tmp+cc).FontSize = 8;
ha(a_tmp+cc).Box     = 'on';
ha(a_tmp+cc).YTickLabelMode = 'auto';
structfun(@(x) set(x(1,:), 'color', colEXP(5,:), ...
    'markeredgecolor', colEXP(5,:)), p);
set([p.box],'linewidth', 1.5);
set([p.med],'linewidth', 1.5);
set([p.lwhis p.uwhis], 'linestyle', '-', 'linewidth', 1.5);
set(p.out, 'marker', 'none');
hold on,
rf=refline(0,0);
rf.Color =  rgb('Black');
ha(a_tmp+cc).YLim = colLim;
tmpTest = dataTmp(end,4,:);
ha(a_tmp+cc).YLim = [min(prctile(tmpTest(:),15),-0.05) prctile(tmpTest(:),85)]%colLim;
% ha(a_tmp+cc).YLim = [min(prctile(tmpTest(:),10),-0.05) max(prctile(tmpTest(:),90),0.05)];
ha(a_tmp+cc).Position = [0.12+(cc-1)*0.47 0.52 .3 0.05];      
%         ax2=axes('Position', get(ha(a_tmp+cc),'outerpos'), 'XColor', rgb('White'), 'XTick', [], 'YTick', [], 'YColor', rgb('White'))
%         uistack(ax2,'down',1)



% C) MSC for global + zones
rr=6;
a_tmp = 12;

varN  = varNames{1};
unitN = unitNames{1};
dObs = dataObs.(varN);
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
        if dN==1
            plot(xSeason,nanmeanArea(dXX(idxZ,:),pix_a(idxZ)), 'Linewidth', 1.5, 'LineStyle', '-.', 'Color', colEXP(dN+1,:)), hold on
        else
            plot(xSeason,nanmeanArea(dXX(idxZ,:),pix_a(idxZ)), 'Linewidth', 1.25, 'LineStyle', '-', 'Color', colEXP(dN+1,:)), hold on
        end
    end
    plot(xSeason,nanmeanArea(dObs(idxZ,:),pix_a(idxZ)), 'Linewidth', 1.75, 'LineStyle', ':', 'Color', colEXP(1,:)), hold on
    plot([0.5 12.5],[0 0], '-', 'color', rgb('Black'))
    set(gca, 'Fontsize', 9)
    title(zNames{cc}, 'FontSize',10)
    if cc==1; yl= ylabel(['TWS [' unitN ']'], 'Fontweight', 'b', 'Fontweight', 'b', 'FontSize', 10); yl.Position(1) = -1.5; end
    ha(cnt).XAxis.MinorTickValues =  [1:1:12];
    ha(a_tmp+cc).XLim = [0.5 12.5];
    
    ha(a_tmp+cc).Position(1) = ha(a_tmp+cc).Position(1)  + 0.03;
    ha(a_tmp+cc).Position(2) = 0.28;
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


rr=7;
a_tmp = 15;

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
        if dN==1
            plot(xSeason,nanmeanArea(dXX(idxZ,:),pix_a(idxZ)), 'Linewidth', 1.5, 'LineStyle', '-.', 'Color', colEXP(dN+1,:)), hold on
        else
            plot(xSeason,nanmeanArea(dXX(idxZ,:),pix_a(idxZ)), 'Linewidth', 1.25, 'LineStyle', '-', 'Color', colEXP(dN+1,:)), hold on
        end
    end
    plot(xSeason,nanmeanArea(dObs(idxZ,:),pix_a(idxZ)), 'Linewidth', 1.75, 'LineStyle', ':', 'Color', colEXP(1,:)), hold on
    plot([0.5 12.5],[0 0], '-', 'color', rgb('Black'))
    set(gca, 'Fontsize', 9)
    title(zNames{cc}, 'FontSize',10)
    if cc==4; yl= ylabel(['TWS [' unitN ']'], 'Fontweight', 'b', 'Fontweight', 'b', 'FontSize', 10); yl.Position(1) = -1.5; end
    ha(a_tmp+cc-3).XAxis.MinorTickValues =  [1:1:12];
    ha(a_tmp+cc-3).XLim = [0.5 12.5];
    
    ha(a_tmp+cc-3).Position(1) = ha(a_tmp+cc-3).Position(1)  + 0.03;
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
        elseif cc == 2 || cc == 3
            t1 = text(0.03,-0.05+0.1*dN, ['MEF = ' num2str(T_MEFMean.(varN){dN,cc})],'Units','Normalized', 'Fontsize', 8, 'Color', colEXP(dN+1,:), 'HorizontalAlignment', 'left');            
        else
            t1 = text(0.97,0.55+0.1*dN, ['MEF = ' num2str(T_MEFMean.(varN){dN,cc})],'Units','Normalized', 'Fontsize', 8, 'Color', colEXP(dN+1,:), 'HorizontalAlignment', 'right');
        end
    end
end

% set axis
set(ha(13:18), 'Xlim', [0.5 12.5], 'XTick', [3:3:12], 'GridAlpha', 0.25, 'XMinorGrid', 'on', 'MinorGridAlpha', 0.05, 'MinorGridLineStyle', '-',...
    'YGrid', 'on')

for tmp=13:18,ha(tmp).YLabel.FontSize = 10; ha(tmp).YAxis.FontSize = 8;ha(tmp).XAxis.FontSize = 8; end

tmpT = ['GRACE TWS' cellfun(@(c)[c '_b_e_s_t'],expNames, 'uni', false)'];
l  = legend(tmpT, 'Position', [.01 0.01 .99 .07], 'box', 'off', 'Orientation', 'Horizontal','Fontsize', 9);

an = annotation('textbox',[0.01 .905 1 .1],'String', 'a)' ,'LineStyle','none', 'Fontsize', 11, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
an = annotation('textbox',[0.12 .905 1 .1],'String', 'Model Effiecency (MEF)' ,'LineStyle','none', 'Fontsize', 11, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
an = annotation('textbox',[0.55 .905 1 .1],'String', ['Difference in MEF with ' dataNames{1} '_b_e_s_t'] ,'LineStyle','none', 'Fontsize', 11, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

% an = annotation('textbox',[0.01 .55 1 .1],'String', 'B)' ,'LineStyle','none', 'Fontsize', 11, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');

an = annotation('textbox',[0.01 .39 1 .1],'String', 'b)' ,'LineStyle','none', 'Fontsize', 11, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
an = annotation('textbox',[0.05 .39 1 .1],'String', 'Mean Seasonal Dynamics' ,'LineStyle','none', 'Fontsize', 11, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');



end
