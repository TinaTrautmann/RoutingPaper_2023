function [sname] = plotTimeSeries2(name, tSteps, time, legendtext, In1, In2, In3, In4, In5, In6, In7)
% Plots time series and MSC of up to 5 variables and returns save name
% name          - for title
% tSteps        - time steps of input data ('d' or 'm')
% time          - startDate & endDate
% legendtext    - cell structure with strings for legend
% In1  - data 1 (red dotted)
% In2  - data 2 (blue line)
% In3  - data 3 (green line)
% In4  - data 4 (orange line)
% In5  - data 5 (black line)
% In6  - data 6 (violet line)
% In7  - data 7 (turquoise line)

name        = strrep(name,'_',' ');
startDate   = time{1};
endDate     = time{2};
if tSteps == 'd'
    xTime        = createDateVector(startDate,endDate , 'd');
    xMonth       = createDateVector(startDate,endDate ,'m');
    [~,vecMonth]    = datevec(xMonth);
    try
        dataM1           = aggDay2Mon(In1,startDate,endDate);
    end
    try
        dataM2           = aggDay2Mon(In2,startDate,endDate);
    end
    try
        dataM3           = aggDay2Mon(In3,startDate,endDate);
    end
    try
        dataM4           = aggDay2Mon(In4,startDate,endDate);
    end
    try
        dataM5           = aggDay2Mon(In5,startDate,endDate);
    end
    try
        dataM6           = aggDay2Mon(In6,startDate,endDate);
    end
    try
        dataM7           = aggDay2Mon(In7,startDate,endDate);
    end

elseif tSteps=='m'
    xTime      = createDateVector(startDate,endDate ,'m');
    [~,vecMonth]    = datevec(xTime);
        try
        dataM1           = In1;
    end
    try
        dataM2           = In2;
    end
    try
        dataM3           = In3;
    end
    try
        dataM4           = In4;
    end
    try
        dataM5           = In5;
    end
    try
        dataM6           = In6;
    end
    try
        dataM7           = In7;
    end

else
    error('no valid timestep given!')
end

xSeason = [1:1:12];

% calc MSC
try
    [~, dataSeason1] = calcMSC(dataM1,vecMonth);
end
try
    [~, dataSeason2] = calcMSC(dataM2,vecMonth);
end
try
    [~, dataSeason3] = calcMSC(dataM3,vecMonth);
end
try
    [~, dataSeason4] = calcMSC(dataM4,vecMonth);
end
try
    [~, dataSeason5] = calcMSC(dataM5,vecMonth);
end
try
    [~, dataSeason6] = calcMSC(dataM6,vecMonth);
end
try
    [~, dataSeason7] = calcMSC(dataM7,vecMonth);
end
% correlation with first input data (daily values)
cor     = NaN(1,7);
try
    cor(1,2) = round(corr(In1(:),In2(:), 'rows', 'complete'),2);
    cor(1,3) = round(corr(In1(:),In3(:), 'rows', 'complete'),2);
    cor(1,4) = round(corr(In1(:),In4(:), 'rows', 'complete'),2);
    cor(1,5) = round(corr(In1(:),In5(:), 'rows', 'complete'),2);
    cor(1,6) = round(corr(In1(:),In6(:), 'rows', 'complete'),2);
    cor(1,7) = round(corr(In1(:),In7(:), 'rows', 'complete'),2);
end

% include correlation in legendtext if In1 exists
if isempty(In1)==1
    leg = legendtext;
else
    leg = [legendtext(1)];
    for i=2:numel(legendtext)
        tmp  = [char(legendtext{i}) ' - cor=' num2str(cor(1,i))];
        leg  = [leg , tmp];
    end
end

leg = strrep(leg, '_','-');

%% Plot
figure('Color', [1 1 1]);
set(gcf, 'Position', [5 5 20 7])

%% Time Series
sp1 = subplot(1,2,1);
try
    plot(xTime,In1, 'r:' , 'Linewidth', 1.5)
end
try
    hold on
    plot(xTime,In2, 'b-' , 'Linewidth', 1)
end
try
    hold on
    plot(xTime,In3, '-' , 'Linewidth', 1, 'color', [0 .7 0])
end
try
    hold on
    plot(xTime,In4, '-' , 'Linewidth', 1, 'color', [.6 0 .4] )
end
try
    hold on
    plot(xTime,In5, 'k-' , 'Linewidth', 1)
end
try
    hold on
    plot(xTime,In6, '-' , 'Linewidth', 1, 'color', [.8 .4 0] )
end
try
    hold on
    plot(xTime,In7, '-' , 'Linewidth', 1, 'color', rgb('MediumTurquoise'))
end
hold on, plot([xTime(1) xTime(end)],[0 0], '-', 'color', [.5 .5 .5]) 
% set axis
box on, grid on
title(['time series (' tSteps ')'], 'fontsize', 9, 'fontweight', 'b');

%% MSC
sp2 = subplot(1,2,2);
try
    plot(xSeason,dataSeason1, 'r:' , 'Linewidth', 1.5)
end
try
    hold on
    plot(xSeason,dataSeason2, 'b-' , 'Linewidth', 1)
end
try
    hold on
    plot(xSeason,dataSeason3, '-' , 'Linewidth', 1, 'color', [0 .7 0])
end
try
    hold on
    plot(xSeason,dataSeason4, '-' , 'Linewidth', 1, 'color', [.6 0 .4] )
end
try
    hold on
    plot(xSeason,dataSeason5, 'k-' , 'Linewidth', 1)
end
try
    hold on
    plot(xSeason,dataSeason6, '-' , 'Linewidth', 1, 'color', [.8 .4 0])
end
try
    hold on
    plot(xSeason,dataSeason7, '-' , 'Linewidth', 1, 'color', rgb('MediumTurquoise'))
end
hold on, plot([0.5 12.5],[0 0], '-', 'color', [.5 .5 .5]) 

% set axis
set(gca, 'Xlim', [0.5 12.5], 'XTick', [1:1:12])
box on, grid on
title('mean seasonal' , 'fontsize', 9, 'fontweight', 'b');

legend(leg, 'Position', [.44 .025 .12 .07], 'box', 'off', 'Orientation', 'Horizontal');

% add overall title

set(sp1,'fontsize',8,'Position', [.04 .22 .58 .64])
set(sp2,'fontsize',8,'Position', [.68 .22 .29 .64])

annotation('textbox', 'String', name, 'fontsize', 10, 'fontweight', 'bold',...
    'edgecolor', 'none', 'Position', [0 0.94 1 0.07], 'HorizontalAlignment', 'center' );


% save name
name        = strrep(name, ',', '_');
name        = strrep(name, '&', '_');
name        = regexprep(name,'\W','');
sname       = name;

end