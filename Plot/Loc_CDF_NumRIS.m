%% theta
clear
clc
set(0,'defaultfigurecolor','w') 
figure; hold on; box on; grid on;
set(gca,'FontName','Times New Roman','FontSize',12);
%% CDL
load('.\DataFinal\Loc_CDF\main_NumRIS_CDL.mat')
p1  = cdfplot(RMSE.theta(1,:));
    set(p1,'LineStyle','--', 'Color', 'r', 'LineWidth', 3)    
p2  = cdfplot(RMSE.theta(2,:));
    set(p2,'LineStyle','--', 'Color', 'g', 'LineWidth', 3)
p3  = cdfplot(RMSE.theta(3,:));
    set(p3,'LineStyle','--', 'Color', 'm', 'LineWidth', 3)
p4  = cdfplot(RMSE.theta(4,:));
    set(p4,'LineStyle','--', 'Color', 'b', 'LineWidth', 3)      
%% PDL
load('.\DataFinal\Loc_CDF\NewPDL\main_NumRIS_PDL.mat')
p5  = cdfplot(RMSE.theta(1,:));
    set(p5,'LineStyle','-', 'Color', 'r', 'LineWidth', 3)
p6  = cdfplot(RMSE.theta(2,:));
    set(p6,'LineStyle','-', 'Color', 'g', 'LineWidth', 3)
p7  = cdfplot(RMSE.theta(3,:));
    set(p7,'LineStyle','-', 'Color', 'm', 'LineWidth', 3)    
p8  = cdfplot(RMSE.theta(4,:));
    set(p8,'LineStyle','-', 'Color', 'b', 'LineWidth', 3)

xlabel('RMSE_{\vartheta}/rad');
ylabel('CDF');
title('')
set(gca, 'XScale', 'log');
% xlim([1e-6 2e-2])
l1 = legend([p1, p2, p3, p4, p5, p6, p7, p8],{  
        'CDL $N_{RIS}$=64',...
        'CDL $N_{RIS}$=128',...
        'CDL $N_{RIS}$=256',...
        'CDL $N_{RIS}$=512',...
        'PDL $N_{RIS}$=64',...
        'PDL $N_{RIS}$=128',...
        'PDL $N_{RIS}$=256',...
        'PDL $N_{RIS}$=512'},...
        'Interpreter','latex', 'Box','off');
l1.FontSize = 12;
%% Range
clear
clc
set(0,'defaultfigurecolor','w') 
figure; hold on; box on; grid on;
set(gca,'FontName','Times New Roman','FontSize',12);
%% CDL
load('.\DataFinal\Loc_CDF\main_NumRIS_CDL.mat')
pp1  = cdfplot(RMSE.r(1,:));
    set(pp1,'LineStyle','--', 'Color', 'r', 'LineWidth', 3)    
pp2  = cdfplot(RMSE.r(2,:));
    set(pp2,'LineStyle','--', 'Color', 'g', 'LineWidth', 3)
pp3  = cdfplot(RMSE.r(3,:));
    set(pp3,'LineStyle','--', 'Color', 'm', 'LineWidth', 3)
pp4  = cdfplot(RMSE.r(4,:));
    set(pp4,'LineStyle','--', 'Color', 'b', 'LineWidth', 3)      
%% PDL
load('.\DataFinal\Loc_CDF\NewPDL\main_NumRIS_PDL.mat')
pp5  = cdfplot(RMSE.r(1,:));
    set(pp5,'LineStyle','-', 'Color', 'r', 'LineWidth', 3)
pp6  = cdfplot(RMSE.r(2,:));
    set(pp6,'LineStyle','-', 'Color', 'g', 'LineWidth', 3)
pp7  = cdfplot(RMSE.r(3,:));
    set(pp7,'LineStyle','-', 'Color', 'm', 'LineWidth', 3)    
pp8  = cdfplot(RMSE.r(4,:));
    set(pp8,'LineStyle','-', 'Color', 'b', 'LineWidth', 3)

xlabel('RMSE_r/m');
ylabel('CDF');
title('')
set(gca, 'XScale', 'log');
% xlim([1e-5 2e-1])
l2 = legend([pp1, pp2, pp3, pp4, pp5, pp6, pp7, pp8],{  
        'CDL $N_{RIS}$=64',...
        'CDL $N_{RIS}$=128',...
        'CDL $N_{RIS}$=256',...
        'CDL $N_{RIS}$=512',...
        'PDL $N_{RIS}$=64',...
        'PDL $N_{RIS}$=128',...
        'PDL $N_{RIS}$=256',...
        'PDL $N_{RIS}$=512'},...
        'Interpreter','latex', 'Box','off');
l2.FontSize = 12;