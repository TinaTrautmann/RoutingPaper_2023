%% S1 Comparison of Calibration Runs for different Experiments with CaMa
% by Tina Trautmann, Aug. 2022

%% set save pth
spth = [pwd '/figures/supplement/'];
if ~exist(spth, 'dir'), mkdir(spth),end


%% load the table with
T_param = readtable([pwd '/data/model_output/Compare_optimizedParams_withCaMa.xls'],'ReadRowNames',true)

expNamesNEW     = T_param.Properties.VariableNames(4:end);
% scale it by max
A_scale = NaN(length(T_param.p_upBound),numel(expNamesNEW));
for pN=1:length(T_param.p_upBound)
    tmp     = T_param{pN,4:end};
    tmpMax  = max(T_param{pN,4:end});
    A_scale(pN,:) = tmp./tmpMax;
end

% assign only a value when the parameter variation is large
tmpVar = var(A_scale(:,1:41),0,2)
A_scale_class = A_scale;
A_scale_class(tmpVar<=0.025,:) = 0.5;

rowNames_Fig = {'s_E_V_I', 'p_i_n_t', 's_R_D_1', 's_R_D_2', 's_R_D_3', 's_R_D_4', 'wSoil_m_a_x_R_D_4', 'p_b_e_r_g', ...
    'rf_S_l_o_w', 'dc_S_l_o_w', 'k_S_o_i_l', '\alpha_V_e_g', 'k_T_r_a_n_s_p', 's_D_e_e_p', 'f_m_a_x', 'd_D_e_e_p'}
col          = othercolor('PRGn5',10);

figure, set(gcf, 'Position', [12 8 20 8])
p=imagesc(A_scale_class)
set(gca, 'XAxisLocation', 'top', 'TickLength', [0.001 0.001], 'XTick', [6, 11+5 , 11+15, 11+25, 11+35], 'XTickLabel', {'MOD', 'MOD-R25', 'MOD-R05', 'MOD-R01', 'MOD-CaMa'})
set(get(gca, 'XAxis'), 'FontWeight', 'bold')
set(gca, 'YTick', [1:1:16], 'YTickLabel', rowNames_Fig)
hold on, plot([11.5 11.5],[0.5 16.5], 'k-')
hold on, plot([21.5 21.5],[0.5 16.5], 'k-')
hold on, plot([31.5 31.5],[0.5 16.5], 'k-')
hold on, plot([41.5 41.5],[0.5 16.5], 'k-')

colormap(col)
cb=colorbar;
cb.Limits = [0 1];
cb.Label.String  = ['calibrated parameter value' char(10) 'normalized by maximum'];
set(gca, 'Position',  [0.2 0.05 0.65 0.8])
anSM = annotation('textbox',[0.05 0.5 0.8 0.25], 'String', ' ')
tSM  = text(-8,7, 'soil water', 'FontWeight', 'b', 'Rotation',90)
anET = annotation('textbox',[0.05 0.20 0.8 0.15], 'String', ' ')
tET  = text(-8,12.5, 'ET', 'FontWeight', 'b', 'Rotation',90)
tCBmax = text(56,0.5, 'high', 'FontWeight', 'b', 'Color', col(end,:))
tCBmin = text(56,16, 'low', 'FontWeight', 'b', 'Color', col(1,:))

% print(gcf,[spth 'S2_optimizedParams_scaled_high_vs_low.png'],'-dpng','-r300');

% calculate variance of param values overall vs all experiments
A_var = NaN(length(T_param.p_upBound),5);
for pN=1:length(T_param.p_upBound)
%     tmpP = T_param{pN,4:end};
    tmpP = A_scale(pN,:);
    %overall
    A_var(pN,1) = var(tmpP,0,2);
    %MOD
    A_var(pN,2) = var(tmpP(1:11),0,2);
    %MOD-25
    A_var(pN,3) = var(tmpP(12:21),0,2);
    %MOD-05
    A_var(pN,4) = var(tmpP(22:31),0,2);
    %MOD-01
    A_var(pN,5) = var(tmpP(32:41),0,2);   
    %MOD-CaMa
    A_var(pN,6) = var(tmpP(42:51),0,2);   
end

% %% look at final costs - infoIn!!!
% allExp = fieldnames(infoIn);
% for eN=1:numel(allExp)
%     expNameNEW  = allExp{eN};
%     pthOptiFull = infoIn.(expNameNEW).info.opti.paths.outFullPath;
%     if ispc==1, pthOptiFull=strrep(pthOptiFull,'/Net/Groups/BGI/','M:\'); end
%     data.(expNameNEW) = load(pthOptiFull);
% end
% 
% data.VEG = data.E_B_bL_RD4;
% data = rmfield(data,'E_B_bL_RD4')


%% cost components
%get the cost components for the minimum total costs
T_cost = readtable([pwd '/data/model_output/Compare_finalCosts_withCaMa.xls'],'ReadRowNames',true)

costCompNames = T_cost.Properties.VariableNames;

%% ... boxplots of costs + param plot combined 
colXX  = [rgb('green'); rgb('Chocolate'); rgb('DarkMagenta'); rgb('Blue'); rgb('Crimson')];

sp  = figure; set(gcf, 'Position', [2 2 18 20]);
ha  = tight_subplot(6,1,[.05 .05],[.05 .05],[.05 .02]);

axes(ha(1))
p=imagesc(A_scale_class)
set(gca, 'XAxisLocation', 'top', 'TickLength', [0.001 0.001], 'XTick', [6, 11+5 , 11+15, 11+25, 11+35], 'XTickLabel', {'MOD', 'MOD-R25', 'MOD-R05', 'MOD-R01', 'MOD-CaMa'})
set(gca, 'YTick', [1:1:16], 'YTickLabel', rowNames_Fig, 'Fontsize', 8)
set(get(gca, 'XAxis'), 'FontWeight', 'bold', 'Fontsize', 10)
hold on, plot([11.5 11.5],[0.5 16.5], 'k-')
hold on, plot([21.5 21.5],[0.5 16.5], 'k-')
hold on, plot([31.5 31.5],[0.5 16.5], 'k-')
hold on, plot([41.5 41.5],[0.5 16.5], 'k-')
colormap(col)
cb=colorbar;
cb.Limits = [0 1];
cb.Label.String  = ['calibrated parameter value' char(10) 'normalized by maximum'];
set(gca, 'Position',  [0.15 0.6 0.7 0.35])
anSM = annotation('rectangle',[0.03 0.7958 0.8202 0.1096])
tSM  = text(-5,6.5, 'soil water', 'FontWeight', 'b', 'FontSize', 9, 'Rotation',90)
anET = annotation('rectangle',[0.03 0.6648 0.8202 0.0661])
tET  = text(-5,12.5, 'ET', 'FontWeight', 'b', 'FontSize', 9, 'Rotation',90)
tCBmax = text(56,0.5, 'high', 'FontWeight', 'b', 'Color', col(end,:))
tCBmin = text(56,16, 'low', 'FontWeight', 'b', 'Color', col(1,:))

% the boxplots of costs
for nC = 1:numel(costCompNames)
    axes(ha(nC+1))
    tmpE = NaN(5,11);
    
    %MOD
    tmpE(1,:) = T_cost{1:11,nC}'
    %MOD-25
    tmpE(2,1:10) = T_cost{12:21,nC}'
    %MOD-05
    tmpE(3,1:10) = T_cost{22:31,nC}'
    %MOD-01
    tmpE(4,1:10) = T_cost{32:41,nC}'
    %MOD-CaMa
    tmpE(5,1:10) = T_cost{42:51,nC}'

    
    p=boxplot(tmpE', 'Widths', 0.25, 'Symbol', 'k.', 'Factorseparator', 1);
    
    a = get(get(gca,'children'),'children');   % Get the handles of all the objects
    t = get(a,'tag');   % List the names of all the objects
    idx=strcmpi(t,'box');  % Find Box objects
    boxes=a(idx);          % Get the children you need
    set(boxes,'linewidth',1, 'Color', rgb('Black')); % Set width
    yL = ylabel([costCompNames{nC}], 'FontWeight', 'b', 'FontSize', 9)
    yL.Position(1) = -0.1
    set(gca, 'Position',  [0.15 0.45-(nC-1)*0.11 0.7 0.09])   
    set(gca, 'XTickLabel', '', 'XLim', [0.4 5.5]) 
    set(gca, 'XMinorGrid', 'off', 'MinorGridLineStyle', '-', 'MinorGridColor', rgb('Black'), 'MinorGridAlpha', 1)
 
end

%adjust YLim
ha(3).YLim(2) = ha(3).YLim(1) + 0.1
ha(4).YLim(2) = ha(4).YLim(1) + 0.02
ha(5).YLim(2) = ha(5).YLim(1) + 0.02
ha(6).YLim(2) = ha(6).YLim(1) + 0.04
set(ha(2:6), 'YTickMode', 'manual')

annotation('textbox',[0.001 0.975 0.5 0.03], 'String', 'a) Calibrated Model Parameters', 'Fontsize', 10', 'Fontweight', 'bold', 'LineStyle', 'none')
annotation('textbox',[0.001 0.55 0.5 0.03], 'String', 'b) Final Calibration Costs', 'Fontsize', 10', 'Fontweight', 'bold', 'LineStyle', 'none')

print(gcf,[spth 'S1_FinalCosts_Parameters_CaMa.png'],'-dpng','-r300');

