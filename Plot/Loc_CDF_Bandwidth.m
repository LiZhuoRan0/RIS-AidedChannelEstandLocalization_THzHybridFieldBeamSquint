%% theta
clear
clc
set(0,'defaultfigurecolor','w') 
figure; hold on; box on; grid on;
set(gca,'FontName','Times New Roman','FontSize',12);
%% CDL
load('.\DataFinal\Loc_CDF\main_Bandwidth_CDL_100e6_500e6.mat')
p1  = cdfplot(RMSE.theta(1,:));
    set(p1,'LineStyle','--', 'Color', 'r', 'LineWidth', 3)    
p2  = cdfplot(RMSE.theta(2,:));
    set(p2,'LineStyle','--', 'Color', 'g', 'LineWidth', 3)

load('.\DataFinal\Loc_CDF\main_Bandwidth_CDL_1e9_5e9_10e9.mat')
p11  = cdfplot(RMSE.theta(1,:));
    set(p11,'LineStyle','--', 'Color', 'c', 'LineWidth', 3)    

load('.\DataFinal\Loc_CDF\main_Bandwidth_CDL_25e8_75e8.mat')
p22  = cdfplot(RMSE.theta(1,:));
    set(p22,'LineStyle','--', 'Color', 'm', 'LineWidth', 3)      

load('.\DataFinal\Loc_CDF\main_Bandwidth_CDL_1e9_5e9_10e9.mat')
p33  = cdfplot(RMSE.theta(2,:));
    set(p33,'LineStyle','--', 'Color', 'k', 'LineWidth', 3)

load('.\DataFinal\Loc_CDF\main_Bandwidth_CDL_25e8_75e8.mat')  
p44  = cdfplot(RMSE.theta(2,:));
    set(p44,'LineStyle',':', 'Color', 'k', 'LineWidth', 3)        

% load('.\DataFinal\Loc_CDF\main_Base_and_NumP_CDL_16.mat')
% p3  = cdfplot(RMSE.theta);
%     set(p3,'LineStyle','--', 'Color', 'b', 'LineWidth', 3)
load('.\DataFinal\Loc_CDF\main_Bandwidth_CDL_1e9_5e9_10e9.mat')
p3  = cdfplot(RMSE.theta(3,:));
    set(p3,'LineStyle','--', 'Color', 'b', 'LineWidth', 3)    
%% PDL
load('.\DataFinal\Loc_CDF\NewPDL\main_Bandwidth_PDL.mat')
p4  = cdfplot(RMSE.theta(1,:));
    set(p4,'LineStyle','-', 'Color', 'r', 'LineWidth', 3)
p5  = cdfplot(RMSE.theta(2,:));
    set(p5,'LineStyle','-', 'Color', 'g', 'LineWidth', 3)
p6  = cdfplot(RMSE.theta(3,:));
    set(p6,'LineStyle','-', 'Color', 'k', 'LineWidth', 3)
p7  = cdfplot(RMSE.theta(4,:));
    set(p7,'LineStyle','-', 'Color', 'c', 'LineWidth', 3)    

load('.\DataFinal\Loc_CDF\main_Base_PDL.mat')
p8  = cdfplot(RMSE.theta);
    set(p8,'LineStyle','-', 'Color', 'b', 'LineWidth', 3)    

xlabel('RMSE_{\vartheta}/rad');
ylabel('CDF');
title('')
set(gca, 'XScale', 'log');
% xlim([1e-8 1e-2])
l1 = legend([p1, p2, p11, p22, p33, p44, p3, p4, p5, p6, p7, p8],{  
        'CDL B=100MHz',...
        'CDL B=500MHz',...
        'CDL B=1GHz',...
        'CDL B=2.5GHz',...
        'CDL B=5GHz',...
        'CDL B=7.5GHz',...
        'CDL B=10GHz',...
        'PDL B=100MHz',...
        'PDL B=500MHz',...
        'PDL B=1GHz',...
        'PDL B=5GHz',...
        'PDL B=10GHz'},...
        'Interpreter','latex', 'Box','off');
l1.FontSize = 12;
%% Range
clear
clc
set(0,'defaultfigurecolor','w') 
figure; hold on; box on; grid on;
set(gca,'FontName','Times New Roman','FontSize',12);
% CDL
load('.\DataFinal\Loc_CDF\main_Bandwidth_CDL.mat')
pp1  = cdfplot(RMSE.r(1,:));
    set(pp1,'LineStyle','--', 'Color', 'r', 'LineWidth', 3)    
pp2  = cdfplot(RMSE.r(2,:));
    set(pp2,'LineStyle','--', 'Color', 'g', 'LineWidth', 3)

load('.\DataFinal\Loc_CDF\main_Bandwidth_CDL_1e9_5e9_10e9.mat')
pp11  = cdfplot(RMSE.r(1,:));
    set(pp11,'LineStyle','--', 'Color', 'c', 'LineWidth', 3)    

load('.\DataFinal\Loc_CDF\main_Bandwidth_CDL_25e8_75e8.mat')
pp22  = cdfplot(RMSE.r(1,:));    
    set(pp22,'LineStyle','--', 'Color', 'm', 'LineWidth', 3)      

load('.\DataFinal\Loc_CDF\main_Bandwidth_CDL_1e9_5e9_10e9.mat')
pp33  = cdfplot(RMSE.r(2,:));
    set(pp33,'LineStyle','--', 'Color', 'k', 'LineWidth', 3)

load('.\DataFinal\Loc_CDF\main_Bandwidth_CDL_25e8_75e8.mat')  
pp44  = cdfplot(RMSE.r(2,:));
    set(pp44,'LineStyle',':', 'Color', 'k', 'LineWidth', 3)    

% load('.\DataFinal\Loc_CDF\main_Base_and_NumP_CDL_16.mat')
% pp3  = cdfplot(RMSE.r);
%     set(pp3,'LineStyle','--', 'Color', 'b', 'LineWidth', 3)
load('.\DataFinal\Loc_CDF\main_Bandwidth_CDL_1e9_5e9_10e9.mat')
pp3  = cdfplot(RMSE.r(3,:));
    set(pp3,'LineStyle','--', 'Color', 'b', 'LineWidth', 3)      
% PDL
load('.\DataFinal\Loc_CDF\NewPDL\main_Bandwidth_PDL.mat')
pp4  = cdfplot(RMSE.r(1,:));
    set(pp4,'LineStyle','-', 'Color', 'r', 'LineWidth', 3)
pp5  = cdfplot(RMSE.r(2,:));
    set(pp5,'LineStyle','-', 'Color', 'g', 'LineWidth', 3)
pp6  = cdfplot(RMSE.r(3,:));
    set(pp6,'LineStyle','-', 'Color', 'k', 'LineWidth', 3)
pp7  = cdfplot(RMSE.r(4,:));
    set(pp7,'LineStyle','-', 'Color', 'm', 'LineWidth', 3)    

load('.\DataFinal\Loc_CDF\main_Base_PDL.mat')
pp8  = cdfplot(RMSE.r);
    set(pp8,'LineStyle','-', 'Color', 'b', 'LineWidth', 3)    

xlabel('RMSE_r/m');
ylabel('CDF');
title('')
set(gca, 'XScale', 'log');
xlim([1e-5 1e-1])
l2 = legend([pp1, pp2, pp11, pp22, pp33, pp44, pp3, pp4, pp5, pp6, pp7, pp8],{  
        'CDL B=100MHz',...
        'CDL B=500MHz',...
        'CDL B=1GHz',...
        'CDL B=2.5GHz',...
        'CDL B=5GHz',...
        'CDL B=7.5GHz',...
        'CDL B=10GHz',...
        'PDL B=100MHz',...
        'PDL B=500MHz',...
        'PDL B=1GHz',...
        'PDL B=5GHz',...
        'PDL B=10GHz'},...
        'Interpreter','latex', 'Box','off');
l2.FontSize = 12;