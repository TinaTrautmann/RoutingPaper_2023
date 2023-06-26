%% Plot MSC of obs & 2 simulations for 3 variables globally & for different zones
function [sname, T_corr, sp, T_MEF, T_RMSE, T_corrMean, T_MEFMean] = plotMSCvarsZones_xExp_bounded(dataNames, dataObs, dataXX, dataUnc, zoneNames, zoneIdx, pix_a, unitNames, col)
% plots the global and zonal MSC of x experiments against observations
% for different variables
%
% returns the save name of the figure, the figure handle and a
% table-structure with the correlation of the 1-2 experiments with
% observations for all zones and variables
%
% Input:
%     dataNames     = experiment names, in the order dataObs, data1, data2
%     dataObs       = structure of observations, fields = variables with MSC for each grid; size(npix,12)
%     dataXX        = cellarray with structures of each experiment, fields = variables with MSC for each grid; size(npix,12)
%     zoneNames     = name of the zones defined, order should reflect the
%                       value in zoneIdx (i.e. first name refers to
%                       zoneIdx = 1
%     zoneIdx       = classifies the grids to different zones, continous
%                       numbers starting from 1; size(npix,1)
%     pix_a         = area of each grid; size(npix,1)
%     unitnames     = cellstring with units of each variable, can be empty
%     col           = line color for plotting obs, experiment1 and
%     experiment 2; default is red, blue, green
%
% Output:
%     sname     = save name (= title of the figure)
%     T_corr    = correlation table for each variable, combined as a
%                   structure
%     sp        = figure handle

% % ----EXAMPLE Input----
% dataNames = {'obs', 'exp1', 'exp2'};
% dataObs.A = randi(5,100,12);
% data1.A   = dataObs.A  ./2 ;
% data2.A   = dataObs.A  + 1.3;
% 
% dataObs.B = randi(5,100,12);
% data1.B   = (dataObs.A ./2 ) + 1;
% data2.B   =  data1.A + 1.2;
% 
% dataObs.C = randi(5,100,12);
% data1.C   = (dataObs.A ./2 ) + 1;
% data2.C   =  data1.A + 1.2;
% 
% zoneNames   = {'zone1','zone2','zone3','zone4','zone5'};
% zoneIdx     = randi(5,100,1);
% 
% pix_a       = ones(100,1);
% unitNames   = {'mm', 'mm/d', 'mm/d'};

%% checks?

%%
if isempty(col)
    col         = [rgb('DarkRed'); rgb('Blue'); rgb('Green')];
end

try
    varNames  = fieldnames(dataObs);
catch
    disp('no dataObs');
end

dataNames = strrep(dataNames,'_','-');
zoneNames = strrep(zoneNames,'_','-');

eN = numel(dataXX);


% varNames  = {'ET', 'Q'}
zNames    = ['Global', zoneNames];
zNames2   = strrep(zNames,'-','');
zNames2   = strrep(zNames2,' ','');

nrows   = numel(varNames);
ncols   = numel(zoneNames)+1;
xSeason = 1:1:12;

% correlations -> 3 tables, 1 per variable
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

sp  = figure; 
if nrows==1
    set(gcf, 'Position', [5 1 ncols*5 nrows*5+2]);
    ha  = tight_subplot(nrows,ncols,[.05 .02],[.15 .15],[.05 .02]);
else
    set(gcf, 'Position', [5 1 ncols*5 nrows*5]);
    ha  = tight_subplot(nrows,ncols,[.05 .02],[.08 .08],[.05 .02]);
end
% loop over rows (=variables), then cols (=global + zones)
cnt = 1;
for rr=1:nrows
    varN  = varNames{rr};
    if strcmp(varN,'wTWS')==1
        varN2 = 'TWS';
    elseif strcmp(varN,'Q')==1
        varN2 = 'Q_R';
    else
        varN2 = strrep(varN,'_','-');
    end
    unitN = unitNames{rr};
    dObs = dataObs.(varN);
    for dN = 1:eN
        tmp = dataXX{dN};
        eval(char(['d' num2str(dN) '_all = tmp.(varN).all;']));
        eval(char(['d' num2str(dN) '_m  = tmp.(varN).minTotal;']));
    end

    for cc=1:ncols
        if cc==1
            idxZ = 1:1:size(dObs,1);
        else
            idxZ = find(zoneIdx==cc-1);
        end
        axes(ha(cnt));
        
        for dN = 1:eN
            eval(char(['dXX = d' num2str(dN) '_m;']))
            eval(char(['dXX_all = d' num2str(dN) '_all;']))
            
            tmpAll = NaN(size(dXX_all,3),12);
            for oN=1:size(dXX_all,3)
                tmpAll(oN,:) = nanmeanArea(dXX_all(idxZ,:,oN),pix_a(idxZ));
            end
            dXX_lo = min(tmpAll,[],1);
            dXX_up = max(tmpAll,[],1);
%             plot(xSeason, nanmeanArea(dXX_lo(idxZ,:),pix_a(idxZ)), 'Color', rgb('LightGray')), hold on
%             plot(xSeason, nanmeanArea(dXX_up(idxZ,:),pix_a(idxZ)), 'Color', rgb('LightGray')), hold on
            x2 = [xSeason, fliplr(xSeason)];
            inBetween = [dXX_lo, fliplr(dXX_up)];
            fill(x2, inBetween, col(dN+1,:), 'FaceAlpha', 0.2, 'EdgeAlpha', 0)
            p=plot(xSeason,nanmeanArea(dXX(idxZ,:),pix_a(idxZ)), 'LineStyle', '-', 'LineWidth', 1.25, 'Color', col(dN+1,:)); hold on
            
            eval(char(['p_' num2str(dN) '= p;']));
        end
        if ~isempty(dataUnc) % && strcmp(varN, 'wTWS')
            tmpUnc =  nanmeanArea(dataUnc.(varN)(idxZ,:),pix_a(idxZ));
            dXX_lo = nanmeanArea(dObs(idxZ,:),pix_a(idxZ)) - tmpUnc;
            dXX_up = nanmeanArea(dObs(idxZ,:),pix_a(idxZ)) + tmpUnc;
            inBetween = [dXX_lo, fliplr(dXX_up)];
            fill(x2, inBetween, col(1,:), 'FaceAlpha', 0.1, 'EdgeAlpha', 0)
        end
        p_obs = plot(xSeason,nanmeanArea(dObs(idxZ,:),pix_a(idxZ)), 'Linewidth', 1.75, 'LineStyle', ':', 'Color', col(1,:)), hold on
        plot([0.5 12.5],[0 0], '-', 'color', [.5 .5 .5])
        if rr==1; title(zNames{cc}); end
        if cc==1; yl= ylabel([varN2 ' [' unitN ']'], 'Fontweight', 'b', 'FontSize', 10); yl.Position(1) = -1.5; end
        ha(cnt).XAxis.MinorTickValues =  [1:1:12];
        cnt = cnt+1;
        
        %calculate correlation
        tmpO = dObs(idxZ,:);
        for dN = 1:eN
            eval(char(['dXX = d' num2str(dN) '_m;']))
            tmp1 = dXX(idxZ,:);
            T_corr.(varN){dN,cc} = round(corr(tmpO(:),tmp1(:), 'rows', 'complete'),2);
            T_MEF.(varN){dN,cc}  = round(calcMEF(tmpO(:),tmp1(:),ones(size(tmp1(:)))),2);
            T_RMSE.(varN){dN,cc} = round(calcRMSE(tmpO(:),tmp1(:)),0);
            T_corrMean.(varN){dN,cc} = round(corr(nanmeanArea(dObs(idxZ,:),pix_a(idxZ))', nanmeanArea(dXX(idxZ,:),pix_a(idxZ))', 'rows', 'complete'),2);
            T_MEFMean.(varN){dN,cc} = round(calcMEF(nanmeanArea(dObs(idxZ,:),pix_a(idxZ)), nanmeanArea(dXX(idxZ,:),pix_a(idxZ)), ones(1,12)),2);
        end
        %add the correlation in the plot
        for dN = 1:eN
%             t1 = text(0.97,-0.025+0.1*dN, ['r^2 = ' num2str(T_corr.(varN){dN,cc})],'Units','Normalized', 'Fontsize', 6, 'Color', col(dN+1,:), 'HorizontalAlignment', 'right');
            if strcmp(varN, 'Q')
                if cc == 1 || cc == 6 
                    t1 = text(0.97,-0.025+0.1*dN, ['MEF = ' num2str(T_MEFMean.(varN){dN,cc})],'Units','Normalized', 'Fontsize', 9, 'Color', col(dN+1,:), 'HorizontalAlignment', 'right');
                elseif cc == 5
                    t1 = text(0.03,0.55+0.1*dN, ['MEF = ' num2str(T_MEFMean.(varN){dN,cc})],'Units','Normalized', 'Fontsize', 9, 'Color', col(dN+1,:), 'HorizontalAlignment', 'left');
                else
                    t1 = text(0.97,0.55+0.1*dN, ['MEF = ' num2str(T_MEFMean.(varN){dN,cc})],'Units','Normalized', 'Fontsize', 9, 'Color', col(dN+1,:), 'HorizontalAlignment', 'right');
                end
            elseif strcmp(varN,'ET')
                if cc == 1 || cc == 6 || cc == 5 || cc == 4
                    t1 = text(0.97,-0.025+0.1*dN, ['MEF = ' num2str(T_MEFMean.(varN){dN,cc})],'Units','Normalized', 'Fontsize', 9, 'Color', col(dN+1,:), 'HorizontalAlignment', 'right');
                else
                    t1 = text(0.97,0.55+0.1*dN, ['MEF = ' num2str(T_MEFMean.(varN){dN,cc})],'Units','Normalized', 'Fontsize', 9, 'Color', col(dN+1,:), 'HorizontalAlignment', 'right');
                end
            elseif strcmp(varN, 'wTWS') && cc == 5
                t1 = text(0.97,-0.025+0.1*dN, ['MEF = ' num2str(T_MEFMean.(varN){dN,cc})],'Units','Normalized', 'Fontsize', 9, 'Color', col(dN+1,:), 'HorizontalAlignment', 'right');
            elseif cc == 2 || cc == 3
                t1 = text(0.03,-0.025+0.1*dN,['MEF = ' num2str(T_MEFMean.(varN){dN,cc})],'Units','Normalized', 'Fontsize', 9, 'Color', col(dN+1,:), 'HorizontalAlignment', 'left');
            else
                t1 = text(0.97,0.55+0.1*dN, ['MEF = ' num2str(T_MEFMean.(varN){dN,cc})],'Units','Normalized', 'Fontsize', 9, 'Color', col(dN+1,:), 'HorizontalAlignment', 'right');
            end
        end
        
%         %correlation with the other experiment
%         tmp1 = d1(idxZ,:);
%         tmp2 = d2(idxZ,:);
%         tmpCorr = round(corr(tmp1(:),tmp2(:), 'rows', 'complete'),2);
%         t2 = text(0.97,-0.025+0.1*3, ['r^2 = ' num2str(tmpCorr)],'Units','Normalized', 'Fontsize', 6, 'Color', rgb('DarkMagenta'), 'HorizontalAlignment', 'right');
         
        
    end
end
% set axis
set(ha(:), 'Xlim', [0.5 12.5], 'XTick', [3:3:12], 'GridAlpha', 0.25, 'XMinorGrid', 'on', 'MinorGridAlpha', 0.05, 'MinorGridLineStyle', '-',...
    'YGrid', 'on', 'Fontsize', 7)
set(ha(1:end-ncols), 'XTickLabel', [])

for tmp=1:length(ha),ha(tmp).YLabel.FontSize = 10;end

% legend & overall title
tmpT = cellfun(@(c)[c '_b_e_s_t'],dataNames, 'uni', false);
tmpT{1} = 'obs';
if nrows==1
    l  = legend([p_obs, p_1, p_2, p_3, p_4], tmpT, 'Position', [.44 0 .12 .07], 'box', 'off', 'Orientation', 'Horizontal','Fontsize', 9);
else                    
    l  = legend([p_obs, p_1, p_2, p_3, p_4], tmpT, 'Position', [.44 -.02 .12 .07], 'box', 'off', 'Orientation', 'Horizontal','Fontsize', 9);
end
an = annotation('textbox',[0 .5 1 .5],'String', 'Mean Seasonal Dynamics' ,'LineStyle','none', 'Fontsize', 10, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

% savename
sname = ['Comparison of the Mean Seasonal Cycle - ' strjoin(tmpT, ' vs ') ];

end