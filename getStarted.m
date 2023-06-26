 ARA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    START     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Workspace
addpath(genpath(pwd))

% Figure Properties
set(0,'DefaultAxesFontName','Arial', 'Units', 'centimeters');
set(0,'DefaultTextFontName','Arial');
set(0,'DefaultFigurePaperUnits','centimeters');
set(0,'DefaultFigureUnits', 'centimeters');
set(0,'DefaultFigurePaperPositionMode','auto');
set(0,'DefaultFigurePosition', [5 5 16 8]);

%%%%%%%%%%%%%%%%%%%%%%%%%   ANALYSIS IN MAIN MANUSCRIPT  %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 0) STUDY AREA
run plot_Study_Area.m
gone

%% I) EFFECT ON MODEL CALIBRATION
% A) effect of river storage on TWS pattern
run run_Effect_TWSpattern_withCaMa.m
gone

% B) effect of river storage on model calibration
run run_Effect_Calibration_withCaMa.m
gone

%% II) EFFECT ON MODEL VALIDATION
run run_Effect_Validation_withCaMa.m
gone


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      SUPPLEMENTS      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% S0) Regional flow accumulation
run S0_regionalFlowAccumulation.m
gone

%% S1) Comparison of calibration runs in terms of calibrated parameters and final costs
run S1_Comparison_CalibrationRuns_withCaMa.m
gone

%% S2) Analysis of selection of calibration grid cells
run S2_Analysis_RelRiver_CalibrationGrids_withCaMa.m
gone

%% S3-4) Influence of selection of calibration grid cells
run S3_Analysis_RelRiver_CalibrationGrids.m
gone

%% S5) Comparison with EartH2Observe model ensemble
run S5_Comparison_EartH2Observe.m
gone

%% S6) Evaluation of discharge at GRDC stations
run S6_Evaluation_Discharge_withCaMa.m
gone

