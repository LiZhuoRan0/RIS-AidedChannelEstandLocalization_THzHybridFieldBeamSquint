%% theta
clear
clc
set(0,'defaultfigurecolor','w') 
figure; hold on; box on; grid on;
set(gca,'FontName','Times New Roman','FontSize',12);
% CDL
load('.\DataFinal\Loc_CDF\main_UEDistance_CDL.mat')
p1  = cdfplot(RMSE.theta(1,:));
    set(p1,'LineStyle','--', 'Color', 'r', 'LineWidth', 3)    
p2  = cdfplot(RMSE.theta(2,:));
    set(p2,'LineStyle','--', 'Color', 'g', 'LineWidth', 3)

load('.\DataFinal\Loc_CDF\main_Base_and_NumP_CDL_16.mat')
p3  = cdfplot(RMSE.theta);
    set(p3,'LineStyle','--', 'Color', 'b', 'LineWidth', 3)
% PDL
load('.\DataFinal\Loc_CDF\NewPDL\main_UEDistance_PDL.mat')
p4  = cdfplot(RMSE.theta(1,:));
    set(p4,'LineStyle','-', 'Color', 'r', 'LineWidth', 3)
p5  = cdfplot(RMSE.theta(2,:));
    set(p5,'LineStyle','-', 'Color', 'g', 'LineWidth', 3) 

load('.\DataFinal\Loc_CDF\main_Base_PDL.mat')
p6  = cdfplot(RMSE.theta);
    set(p6,'LineStyle','-', 'Color', 'b', 'LineWidth', 3)    

xlabel('RMSE_{\vartheta}/rad');
ylabel('CDF');
title('')
set(gca, 'XScale', 'log');
% xlim([1e-8 1e-2])
l1 = legend([p1, p2, p3, p4, p5, p6],{  
        'CDL $r_0^{BU}$=11.2m',...
        'CDL $r_0^{BU}$=16.8m',...
        'CDL $r_0^{BU}$=22.4m',...
        'PDL $r_0^{BU}$=11.2m',...
        'PDL $r_0^{BU}$=16.8m',...
        'PDL $r_0^{BU}$=22.4m'},...
        'Interpreter','latex', 'Box','off');
l1.FontSize = 12;
%% Range
clear
clc
set(0,'defaultfigurecolor','w') 
figure; hold on; box on; grid on;
set(gca,'FontName','Times New Roman','FontSize',12);
% CDL
load('.\DataFinal\Loc_CDF\main_UEDistance_CDL.mat')
pp1  = cdfplot(RMSE.r(1,:));
    set(pp1,'LineStyle','--', 'Color', 'r', 'LineWidth', 3)    
pp2  = cdfplot(RMSE.r(2,:));
    set(pp2,'LineStyle','--', 'Color', 'g', 'LineWidth', 3)

load('.\DataFinal\Loc_CDF\main_Base_and_NumP_CDL_16.mat')
pp3  = cdfplot(RMSE.r);
    set(pp3,'LineStyle','--', 'Color', 'b', 'LineWidth', 3)
% PDL
load('.\DataFinal\Loc_CDF\NewPDL\main_UEDistance_PDL.mat')
pp4  = cdfplot(RMSE.r(1,:));
    set(pp4,'LineStyle','-', 'Color', 'r', 'LineWidth', 3)
pp5  = cdfplot(RMSE.r(2,:));
    set(pp5,'LineStyle','-', 'Color', 'g', 'LineWidth', 3)

load('.\DataFinal\Loc_CDF\main_Base_PDL.mat')
pp6  = cdfplot(RMSE.r);
    set(pp6,'LineStyle','-', 'Color', 'b', 'LineWidth', 3)    

xlabel('RMSE_r/m');
ylabel('CDF');
title('')
set(gca, 'XScale', 'log');
% xlim([1e-4 1e-1])
l2 = legend([pp1, pp2, pp3, pp4, pp5, pp6],{  
        'CDL $r_0^{BU}$=11.2m',...
        'CDL $r_0^{BU}$=16.8m',...
        'CDL $r_0^{BU}$=22.4m',...
        'PDL $r_0^{BU}$=11.2m',...
        'PDL $r_0^{BU}$=16.8m',...
        'PDL $r_0^{BU}$=22.4m'},...
        'Interpreter','latex', 'Box','off');
l2.FontSize = 12;