%% theta
clear
clc
set(0,'defaultfigurecolor','w') 
figure; hold on; box on; grid on;
set(gca,'FontName','Times New Roman','FontSize',12);
%% CDL
load('.\DataFinal\Loc_CDF\main_NumBS_CDL_64_128_512.mat')
p1  = cdfplot(RMSE.theta(1,:));
    set(p1,'LineStyle','--', 'Color', 'r', 'LineWidth', 3)    
p2  = cdfplot(RMSE.theta(2,:));
    set(p2,'LineStyle','--', 'Color', 'g', 'LineWidth', 3)

% load('.\DataFinal\Loc_CDF\main_Base_and_NumP_CDL_16.mat')
% p3  = cdfplot(RMSE.theta);
%     set(p3,'LineStyle','--', 'Color', 'm', 'LineWidth', 3)
% 192 256 384
load('.\DataFinal\Loc_CDF\main_NumBS_CDL_192_256_384.mat')
p31  = cdfplot(RMSE.theta(1,:));
    set(p31,'LineStyle','-.', 'Color', 'm', 'LineWidth', 3)    
p32  = cdfplot(RMSE.theta(2,:));
    set(p32,'LineStyle','--', 'Color', 'm', 'LineWidth', 3)
p33  = cdfplot(RMSE.theta(3,:));
    set(p33,'LineStyle',':', 'Color', 'm', 'LineWidth', 3)    

load('.\DataFinal\Loc_CDF\main_NumBS_CDL_64_128_512.mat')
p4  = cdfplot(RMSE.theta(3,:));
    set(p4,'LineStyle','--', 'Color', 'b', 'LineWidth', 3)      
%% PDL
load('.\DataFinal\Loc_CDF\NewPDL\main_NumBS_PDL.mat')
p5  = cdfplot(RMSE.theta(1,:));
    set(p5,'LineStyle','-', 'Color', 'r', 'LineWidth', 3)
p6  = cdfplot(RMSE.theta(2,:));
    set(p6,'LineStyle','-', 'Color', 'g', 'LineWidth', 3)

load('.\DataFinal\Loc_CDF\main_Base_PDL.mat')
p7  = cdfplot(RMSE.theta);
    set(p7,'LineStyle','-', 'Color', 'm', 'LineWidth', 3)    

load('.\DataFinal\Loc_CDF\NewPDL\main_NumBS_PDL.mat')
p8  = cdfplot(RMSE.theta(3,:));
    set(p8,'LineStyle','-', 'Color', 'b', 'LineWidth', 3)

xlabel('RMSE_{\vartheta}/rad');
ylabel('CDF');
title('')
set(gca, 'XScale', 'log');
xlim([1e-9 1e-0])
l1 = legend([p1, p2, p31, p32, p33, p4, p5, p6, p7, p8],{  
        'CDL $N$=64',...
        'CDL $N$=128',...
        'CDL $N$=192',...
        'CDL $N$=256',...
        'CDL $N$=384',...
        'CDL $N$=512',...
        'PDL $N$=64',...
        'PDL $N$=128',...
        'PDL $N$=256',...
        'PDL $N$=512'},...
        'Interpreter','latex', 'Box','off');
l1.FontSize = 12;
%% Range
clear
clc
set(0,'defaultfigurecolor','w') 
figure; hold on; box on; grid on;
set(gca,'FontName','Times New Roman','FontSize',12);
%% CDL
load('.\DataFinal\Loc_CDF\main_NumBS_CDL_64_128_512.mat')
pp1  = cdfplot(RMSE.r(1,:));
    set(pp1,'LineStyle','--', 'Color', 'r', 'LineWidth', 3)    
pp2  = cdfplot(RMSE.r(2,:));
    set(pp2,'LineStyle','--', 'Color', 'g', 'LineWidth', 3)

% load('.\DataFinal\Loc_CDF\main_Base_and_NumP_CDL_16.mat')
% pp3  = cdfplot(RMSE.r);
%     set(pp3,'LineStyle','--', 'Color', 'm', 'LineWidth', 3)
load('.\DataFinal\Loc_CDF\main_NumBS_CDL_192_256_384.mat')
pp31  = cdfplot(RMSE.r(1,:));
    set(pp31,'LineStyle','-.', 'Color', 'm', 'LineWidth', 3)    
pp32  = cdfplot(RMSE.r(2,:));
    set(pp32,'LineStyle','--', 'Color', 'm', 'LineWidth', 3)
pp33  = cdfplot(RMSE.r(3,:));
    set(pp33,'LineStyle',':', 'Color', 'm', 'LineWidth', 3)        

load('.\DataFinal\Loc_CDF\main_NumBS_CDL_64_128_512.mat')
pp4  = cdfplot(RMSE.r(3,:));
    set(pp4,'LineStyle','--', 'Color', 'b', 'LineWidth', 3)      
%% PDL
load('.\DataFinal\Loc_CDF\NewPDL\main_NumBS_PDL.mat')
pp5  = cdfplot(RMSE.r(1,:));
    set(pp5,'LineStyle','-', 'Color', 'r', 'LineWidth', 3)
pp6  = cdfplot(RMSE.r(2,:));
    set(pp6,'LineStyle','-', 'Color', 'g', 'LineWidth', 3)

load('.\DataFinal\Loc_CDF\main_Base_PDL.mat')
pp7  = cdfplot(RMSE.r);
    set(pp7,'LineStyle','-', 'Color', 'm', 'LineWidth', 3)    

load('.\DataFinal\Loc_CDF\NewPDL\main_NumBS_PDL.mat')
pp8  = cdfplot(RMSE.r(3,:));
    set(pp8,'LineStyle','-', 'Color', 'b', 'LineWidth', 3)

xlabel('RMSE_r/m');
ylabel('CDF');
title('')
set(gca, 'XScale', 'log');
% xlim([1e-8 1e-2])
l2 = legend([pp1, pp2, pp31, pp32, pp33, pp4, pp5, pp6, pp7, pp8],{  
        'CDL $N$=64',...
        'CDL $N$=128',...
        'CDL $N$=192',...
        'CDL $N$=256',...
        'CDL $N$=384',...
        'CDL $N$=512',...
        'PDL $N$=64',...
        'PDL $N$=128',...
        'PDL $N$=256',...
        'PDL $N$=512'},...
        'Interpreter','latex', 'Box','off');
l2.FontSize = 12;