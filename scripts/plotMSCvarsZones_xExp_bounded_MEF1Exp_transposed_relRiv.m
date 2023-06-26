%% Plot MSC of obs & 2 simulations for 3 variables globally & for different zones
function [sname, T_corr, sp, T_MEF, T_RMSE, T_corrMean, T_MEFMean] = plotMSCvarsZones_xExp_bounded_MEF1Exp_transposed_relRiv(dataNames, dataObs, dataXX, dataUnc, zoneNames, zoneIdx, pix_a, unitNames, col)
% plots the global and zonal MSC of x experiments against observations
% for different variables, calculates metrics not against observations, but
% against the 1st Experiment
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

dataNames2 = strrep(dataNames,'-','');

ncols   = numel(varNames);
nrows   = numel(zoneNames)+1;
xSeason = 1:1:12;

% correlations -> 3 tables, 1 per variable
tmpArray    = NaN(numel(zNames),numel(dataNames)-2);
varTypes    = cell(numel(zNames),1);
varTypes(:) = {'double'};
for vn=1:numel(varNames)
    T_corr.(varNames{vn}) = array2table(tmpArray, 'VariableNames', dataNames2(3:end),'RowNames', zNames2);
    T_MEF.(varNames{vn})  = array2table(tmpArray, 'VariableNames', dataNames2(3:end),'RowNames', zNames2); 
    T_RMSE.(varNames{vn}) = array2table(tmpArray, 'VariableNames', dataNames2(3:end),'RowNames', zNames2);
    T_corrMean.(varNames{vn}) = array2table(tmpArray, 'VariableNames', dataNames2(3:end),'RowNames', zNames2);
    T_MEFMean.(varNames{vn})  = array2table(tmpArray, 'VariableNames', dataNames2(3:end),'RowNames', zNames2);
end

sp  = figure; 
if nrows==1
    set(gcf, 'Position', [5 1 ncols*5 nrows*5+2]);
    ha  = tight_subplot(nrows,ncols,[.05 .02],[.15 .15],[.05 .02]);
else
    set(gcf, 'Position', [5 1 ncols*5 nrows*5]);
    ha  = tight_subplot(nrows,ncols,[.05 .02],[.08 .08],[.05 .02]);
end
% loop over columns (=variables), then rows (=global + zones)
for cc=1:ncols
    cnt = cc;
    varN  = varNames{cc};
    if strcmp(varN,'wTWS')==1
        varN2 = 'TWS';
    elseif strcmp(varN,'Q')==1
        varN2 = 'Q_R';
    else
        varN2 = strrep(varN,'_','-');
    end

    unitN = unitNames{cc};
    dObs = dataObs.(varN);
    for dN = 1:eN
        tmp = dataXX{dN};
        eval(char(['d' num2str(dN) ' = tmp.(varN);']));
    end
    
    for rr=1:nrows
       axes(ha(cnt));
        
       if rr==1
           idxZ = 1:1:size(dObs,1);
       else
           idxZ = find(zoneIdx==rr-1);
       end

        for dN = 1:eN
            eval(char(['dXX = d' num2str(dN) ';']))
            
%             plot(xSeason, nanmeanArea(dXX_lo(idxZ,:),pix_a(idxZ)), 'Color', rgb('LightGray')), hold on
%             plot(xSeason, nanmeanArea(dXX_up(idxZ,:),pix_a(idxZ)), 'Color', rgb('LightGray')), hold on
            if dN==1
                p=plot(xSeason,nanmeanArea(dXX(idxZ,:),pix_a(idxZ)), 'LineStyle', '-.', 'LineWidth', 1.75, 'Color', col(dN+1,:)); hold on
            else
                p=plot(xSeason,nanmeanArea(dXX(idxZ,:),pix_a(idxZ)), 'LineStyle', '-', 'LineWidth', 1.25, 'Color', col(dN+1,:)); hold on
            end
            eval(char(['p_' num2str(dN) '= p;']));
        end
        p_obs = plot(xSeason,nanmeanArea(dObs(idxZ,:),pix_a(idxZ)), 'Linewidth', 1.25, 'LineStyle', ':', 'Color', col(1,:)), hold on
        plot([0.5 12.5],[0 0], '-', 'color', [.5 .5 .5])
        if rr==1; title([varN2 ' [' unitN ']'], 'Fontsize', 10); end
        if cc==1; yl= ylabel(zNames{rr}, 'Fontsize', 10, 'Fontweight', 'b'); yl.Position(1) = -1.5; end
        ha(cnt).XAxis.MinorTickValues =  [1:1:12];
        
        ha(cnt).Position(1) = ha(cnt).Position(1)+0.02;
        ha(cnt).Position(3) = 0.28;
        ha(cnt).Position(4) = 0.12;
        
        %calculate correlation
        tmpO = d1(idxZ,:);
        
        for dN = 2:eN
            eval(char(['dXX = d' num2str(dN) ';']))
            tmp1 = dXX(idxZ,:);
            T_corr.(varN){dN,rr} = round(corr(tmpO(:),tmp1(:), 'rows', 'complete'),2);
            T_MEF.(varN){dN,rr}  = round(calcMEF(tmpO(:),tmp1(:),ones(size(tmp1(:)))),2);
            T_RMSE.(varN){dN,rr} = round(calcRMSE(tmpO(:),tmp1(:)),0);
            T_corrMean.(varN){dN,rr} = round(corr(nanmeanArea(d1(idxZ,:),pix_a(idxZ))', nanmeanArea(dXX(idxZ,:),pix_a(idxZ))', 'rows', 'complete'),2);
            T_MEFMean.(varN){dN,rr} = round(calcMEF(nanmeanArea(d1(idxZ,:),pix_a(idxZ)), nanmeanArea(dXX(idxZ,:),pix_a(idxZ)), ones(1,12)),2);
        end
        %add the correlation in the plot
        for dN = 2:eN
%             t1 = text(0.97,-0.025+0.1*dN, ['r^2 = ' num2str(T_corr.(varN){dN,cc})],'Units','Normalized', 'Fontsize', 6, 'Color', col(dN+1,:), 'HorizontalAlignment', 'right');
            if strcmp(varN, 'Q')
                if rr == 1 || rr == 6 
                    t1 = text(0.97,-0.125+0.1*dN, ['MEF = ' num2str(T_MEFMean.(varN){dN,rr})],'Units','Normalized', 'Fontsize', 9, 'Color', col(dN+1,:), 'HorizontalAlignment', 'right');
                elseif rr == 5
                    t1 = text(0.03,0.45+0.1*dN, ['MEF = ' num2str(T_MEFMean.(varN){dN,rr})],'Units','Normalized', 'Fontsize', 9, 'Color', col(dN+1,:), 'HorizontalAlignment', 'left');
                else
                    t1 = text(0.97,0.45+0.1*dN, ['MEF = ' num2str(T_MEFMean.(varN){dN,rr})],'Units','Normalized', 'Fontsize', 9, 'Color', col(dN+1,:), 'HorizontalAlignment', 'right');
                end
            elseif strcmp(varN,'ET')
                if rr == 1 || rr == 6 || rr == 5 || rr == 4
                    t1 = text(0.97,-0.125+0.1*dN, ['MEF = ' num2str(T_MEFMean.(varN){dN,rr})],'Units','Normalized', 'Fontsize', 9, 'Color', col(dN+1,:), 'HorizontalAlignment', 'right');
                else
                    t1 = text(0.97,0.45+0.1*dN, ['MEF = ' num2str(T_MEFMean.(varN){dN,rr})],'Units','Normalized', 'Fontsize', 9, 'Color', col(dN+1,:), 'HorizontalAlignment', 'right');
                end
            elseif strcmp(varN, 'wTWS') && rr == 5
                t1 = text(0.97,-0.125+0.1*dN, ['MEF = ' num2str(T_MEFMean.(varN){dN,rr})],'Units','Normalized', 'Fontsize', 9, 'Color', col(dN+1,:), 'HorizontalAlignment', 'right');
            elseif rr == 2 || rr == 3
                t1 = text(0.03,-0.125+0.1*dN,['MEF = ' num2str(T_MEFMean.(varN){dN,rr})],'Units','Normalized', 'Fontsize', 9, 'Color', col(dN+1,:), 'HorizontalAlignment', 'left');
            else
                t1 = text(0.97,0.45+0.1*dN, ['MEF = ' num2str(T_MEFMean.(varN){dN,rr})],'Units','Normalized', 'Fontsize', 9, 'Color', col(dN+1,:), 'HorizontalAlignment', 'right');
            end
        end
        
        cnt = cnt+ncols;
        
    end
    
end
% set axis
set(ha(:), 'Xlim', [0.5 12.5], 'XTick', [3:3:12], 'GridAlpha', 0.25, 'XMinorGrid', 'on', 'MinorGridAlpha', 0.05, 'MinorGridLineStyle', '-')
set(ha(1:end-ncols), 'XTickLabel', [])

for tmp=1:length(ha),ha(tmp).YLabel.FontSize = 10;end

% legend & overall title
tmpT = dataNames;
tmpT{1} = 'obs';
tmpT{2} = [dataNames{2} '_b_e_s_t'];
if nrows==1
    l  = legend([p_obs, p_1, p_2, p_3, p_4, p_5], tmpT, 'Position', [.44 0 .12 .07], 'box', 'off','Fontsize', 9, 'NumColumns', 3);
else                    
    l  = legend([p_obs, p_1, p_2, p_3, p_4, p_5], tmpT, 'Position', [.44 0 .12 .05], 'box', 'off','Fontsize', 9, 'NumColumns', 3);
end
an = annotation('textbox',[0 .5 1 .5],'String', 'Mean Seasonal Dynamics' ,'LineStyle','none', 'Fontsize', 11, 'Fontweight', 'b','fitboxtotext', 'on', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center');

% savename
sname = ['Comparison of the Mean Seasonal Cycle - ' strjoin(tmpT, ' vs ') ];

end