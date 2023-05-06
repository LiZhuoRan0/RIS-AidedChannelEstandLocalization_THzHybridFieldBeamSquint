%% theta
clear
clc
set(0,'defaultfigurecolor','w') 
figure; hold on; box on; grid on;
set(gca,'FontName','Times New Roman','FontSize',12);
%% CDL
load('.\DataFinal\Loc_CDF\main_NumP_CDL.mat')
p1  = cdfplot(RMSE.theta(1,:));
    set(p1,'LineStyle','--', 'Color', 'b', 'LineWidth', 3)    
load('.\DataFinal\Loc_CDF\main_Base_and_NumP_CDL_16.mat')
p2  = cdfplot(RMSE.theta);
    set(p2,'LineStyle','--', 'Color', 'm', 'LineWidth', 3)  
load('.\DataFinal\Loc_CDF\main_NumP_CDL_32.mat')    
p3  = cdfplot(RMSE.theta(1,:));
    set(p3,'LineStyle','--', 'Color', 'g', 'LineWidth', 3)

 
%% PDL
load('.\DataFinal\Loc_CDF\main_Base_PDL.mat')
p4  = cdfplot(RMSE.theta);
    set(p4,'LineStyle','-', 'Color', 'b', 'LineWidth', 3)    
        
load('.\DataFinal\Loc_CDF\NewPDL\main_NumP_PDL.mat')
p5  = cdfplot(RMSE.theta(1,:));
    set(p5,'LineStyle','-', 'Color', 'm', 'LineWidth', 3)

load('.\DataFinal\Loc_CDF\NewPDL\main_NumP_PDL.mat')
p6  = cdfplot(RMSE.theta(2,:));
    set(p6,'LineStyle','-', 'Color', 'g', 'LineWidth', 3)
xlabel('RMSE_{\vartheta}/rad');
ylabel('CDF');
title('')
set(gca, 'XScale', 'log');
% xlim([1e-8 1e-2])
l1 = legend([p1, p2, p3, p4, p5, p6],{  
        'CDL $P^{NRIS}$=8',...
        'CDL $P^{NRIS}$=16',...
        'CDL $P^{NRIS}$=32',...
        'PDL $P^{NRIS}$=8',...
        'PDL $P^{NRIS}$=16',...
        'PDL $P^{NRIS}$=32'},...
        'Interpreter','latex', 'Box','off');
l1.FontSize = 12;
%% Range
clear
clc
set(0,'defaultfigurecolor','w') 
figure; hold on; box on; grid on;
set(gca,'FontName','Times New Roman','FontSize',12);
%% CDL
load('.\DataFinal\Loc_CDF\main_NumP_CDL.mat')
p1  = cdfplot(RMSE.r(1,:));
    set(p1,'LineStyle','--', 'Color', 'b', 'LineWidth', 3)    
load('.\DataFinal\Loc_CDF\main_Base_and_NumP_CDL_16.mat')
p2  = cdfplot(RMSE.r);
    set(p2,'LineStyle','--', 'Color', 'm', 'LineWidth', 3)  
load('.\DataFinal\Loc_CDF\main_NumP_CDL_32.mat')    
p3  = cdfplot(RMSE.r(1,:));
    set(p3,'LineStyle','--', 'Color', 'g', 'LineWidth', 3)

 
%% PDL
load('.\DataFinal\Loc_CDF\main_Base_PDL.mat')
p4  = cdfplot(RMSE.r);
    set(p4,'LineStyle','-', 'Color', 'b', 'LineWidth', 3)    

load('.\DataFinal\Loc_CDF\NewPDL\main_NumP_PDL.mat')
p5  = cdfplot(RMSE.r(1,:));
    set(p5,'LineStyle','-', 'Color', 'm', 'LineWidth', 3)

load('.\DataFinal\Loc_CDF\NewPDL\main_NumP_PDL.mat')
p6  = cdfplot(RMSE.r(2,:));
    set(p6,'LineStyle','-', 'Color', 'g', 'LineWidth', 3)
xlabel('RMSE_r/m');
ylabel('CDF');
title('')
set(gca, 'XScale', 'log');
% xlim([1e-8 1e-2])
l1 = legend([p1, p2, p3, p4, p5, p6],{  
        'CDL $P^{NRIS}$=8',...
        'CDL $P^{NRIS}$=16',...
        'CDL $P^{NRIS}$=32',...
        'PDL $P^{NRIS}$=8',...
        'PDL $P^{NRIS}$=16',...
        'PDL $P^{NRIS}$=32'},...
        'Interpreter','latex', 'Box','off');
l1.FontSize = 12;