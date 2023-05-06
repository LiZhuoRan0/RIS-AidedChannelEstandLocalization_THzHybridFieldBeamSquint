%% theta
clear
clc
set(0,'defaultfigurecolor','w') 
figure; hold on; box on; grid on;
set(gca,'FontName','Times New Roman','FontSize',12);

load('.\DataFinal\Loc_CDF\main_Base_and_NumP_CDL_16.mat')
p1  = cdfplot(RMSE.theta_co);
    set(p1,'LineStyle','--', 'Color', 'r', 'LineWidth', 3)    
p2  = cdfplot(RMSE.theta);
    set(p2,'LineStyle','--', 'Color', 'g', 'LineWidth', 3)
% load('.\DataFinal\Loc_CDF\main_Base_PDL.mat')
% p3  = cdfplot(RMSE.theta);
%     set(p3,'LineStyle','--', 'Color', 'm', 'LineWidth', 3)
load('.\DataFinal\Loc_CDF\NewPDL\main_Base_PDL.mat')
p3  = cdfplot(RMSE.theta);
    set(p3,'LineStyle','--', 'Color', 'm', 'LineWidth', 3)

load('.\DataFinal\Loc_CDF\NewPDL\main_Base_PDL_Baseline.mat')
p4  = cdfplot(RMSE.theta_full);
    set(p4,'LineStyle','--', 'Color', 'b', 'LineWidth', 3)
p5  = cdfplot(RMSE.theta_full_ESPRIT);
    set(p5,'LineStyle','--', 'Color', 'c', 'LineWidth', 3)

xlabel('RMSE_{\vartheta}/rad');
ylabel('CDF');
title('')
set(gca, 'XScale', 'log');
% xlim([1e-7 1e-2])
l1 = legend([p1, p2, p3, p4, p5],{  
        'GMMV-OMP',...
        'CDL',...
        'PDL',...
        'MUSIC',...
        'ESPRIT'},...
        'Interpreter','latex', 'Box','off');
l1.FontSize = 12;
%% r
set(0,'defaultfigurecolor','w') 
figure; hold on; box on; grid on;
set(gca,'FontName','Times New Roman','FontSize',12);

load('.\DataFinal\Loc_CDF\main_Base_and_NumP_CDL_16.mat')
p6  = cdfplot(RMSE.r_co);
    set(p6,'LineStyle','--', 'Color', 'r', 'LineWidth', 3)    
p7  = cdfplot(RMSE.r);
    set(p7,'LineStyle','--', 'Color', 'g', 'LineWidth', 3)
% load('.\DataFinal\Loc_CDF\main_Base_PDL.mat')
% p8  = cdfplot(RMSE.r);
%     set(p8,'LineStyle','--', 'Color', 'm', 'LineWidth', 3)
load('.\DataFinal\Loc_CDF\NewPDL\main_Base_PDL.mat')
p8  = cdfplot(RMSE.r);
    set(p8,'LineStyle','--', 'Color', 'm', 'LineWidth', 3)

load('.\DataFinal\Loc_CDF\NewPDL\main_Base_PDL_Baseline.mat')
p9  = cdfplot(RMSE.r_full);
    set(p9,'LineStyle','--', 'Color', 'b', 'LineWidth', 3)
p10 = cdfplot(RMSE.r_full_ESPRIT);
    set(p10,'LineStyle','--', 'Color', 'c', 'LineWidth', 3)

xlabel('RMSE_{r}/m');
ylabel('CDF');
title('')
set(gca, 'XScale', 'log');
% xlim([1e-4 0.3])
l2 = legend([p6, p7, p8, p9, p10],{  
        'GMMV-OMP',...
        'CDL',...
        'PDL',...
        'MUSIC',...
        'ESPRIT'},...
        'Interpreter','latex', 'Box','off');
l2.FontSize = 12;