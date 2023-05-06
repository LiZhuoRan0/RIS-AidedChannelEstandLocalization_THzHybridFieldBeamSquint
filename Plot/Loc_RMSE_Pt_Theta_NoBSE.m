%% theta
clear
clc
set(0,'defaultfigurecolor','w') 
figure; hold on; box on; grid on;
set(gca,'FontName','Times New Roman','FontSize',12);
xlabel('P/dBm');
ylabel('RMSE_{\vartheta}/rad');
set(gca, 'YScale', 'log');
ylim([1e-6 1e1])
P_t_dBm_ALL = [15 20 25 30 35 40 45 50 60];
%% -    without BSE Near
load('.\DataFinal\Loc_RMSE_Pt\RMSE_Pt_MUSICESPRIT_Near_NoBSE_60dBm.mat')
theta_full = RMSE.theta_full;
theta_full_ESPRIT = RMSE.theta_full_ESPRIT;
load('.\DataFinal\Loc_RMSE_Pt\RMSE_Pt_MUSICESPRIT_Near_NoBSE.mat')
p1 = plot(P_t_dBm_ALL(1:2:end), [RMSE.theta_full(1:2:end) theta_full],'r:o', 'LineWidth', 1.5,'MarkerSize',10);
p2 = plot(P_t_dBm_ALL(1:2:end), [RMSE.theta_full_ESPRIT(1:2:end) theta_full_ESPRIT],'r:x', 'LineWidth', 1.5,'MarkerSize',10);
%% -    without BSE Far
load('.\DataFinal\Loc_RMSE_Pt\RMSE_Pt_MUSICESPRIT_Far_NoBSE_60dBm.mat')
theta_full = RMSE.theta_full;
theta_full_ESPRIT = RMSE.theta_full_ESPRIT;
load('.\DataFinal\Loc_RMSE_Pt\RMSE_Pt_MUSICESPRIT_Far_NoBSE.mat')
p11 = plot(P_t_dBm_ALL(1:2:end), [RMSE.theta_full(1:2:end) theta_full],'g:o', 'LineWidth', 1.5,'MarkerSize',10);
p22 = plot(P_t_dBm_ALL(1:2:end), [RMSE.theta_full_ESPRIT(1:2:end) theta_full_ESPRIT],'g:x', 'LineWidth', 1.5,'MarkerSize',10);
%% -    with BSE Near
load('.\DataFinal\Loc_RMSE_Pt\RMSE_Pt_MUSICESPRIT_Near_55dBm.mat')
theta_full = RMSE.theta_full;
theta_full_ESPRIT = RMSE.theta_full_ESPRIT;
load('.\DataFinal\Loc_RMSE_Pt\RMSE_Pt_MUSICESPRIT_Near.mat')
p3 = plot(P_t_dBm_ALL(1:2:end), [RMSE.theta_full(1:2:end) theta_full],'b:o', 'LineWidth', 1.5,'MarkerSize',10);
p4 = plot(P_t_dBm_ALL(1:2:end), [RMSE.theta_full_ESPRIT(1:2:end) theta_full_ESPRIT],'b:x', 'LineWidth', 1.5,'MarkerSize',10);
%% -    with BSE Far
load('.\DataFinal\Loc_RMSE_Pt\RMSE_Pt_MUSICESPRIT_Far_55dBm.mat')
theta_full = RMSE.theta_full;
theta_full_ESPRIT = RMSE.theta_full_ESPRIT;
load('.\DataFinal\Loc_RMSE_Pt\RMSE_Pt_MUSICESPRIT_Far.mat')
p33 = plot(P_t_dBm_ALL(1:2:end), [RMSE.theta_full(1:2:end) theta_full],'m:o', 'LineWidth', 1.5,'MarkerSize',10);
p44 = plot(P_t_dBm_ALL(1:2:end), [RMSE.theta_full_ESPRIT(1:2:end) theta_full_ESPRIT],'m:x', 'LineWidth', 1.5,'MarkerSize',10);
%% -    legend
l1 = legend([p1, p2],{  
        'Near-Field~ MUSIC~~B=500~MHz',...
        'Near-Field~ ESPRIT~B=500~MHz'},...
        'Interpreter','latex', 'Box','off');
l1.FontSize = 12;    
ah1=axes('position',get(gca,'position'),'visible','off');    
l2 = legend(ah1, [p3, p4],{  
        'Near-Field~ MUSIC~~B=10~GHz',...
        'Near-Field~ ESPRIT~B=10~GHz'},...    
        'Interpreter','latex', 'Box','off');    
l2.FontSize = 12; 
ah2=axes('position',get(gca,'position'),'visible','off'); 
l3 = legend(ah2, [p11, p22],{  
        'Far-Field~ MUSIC~~B=500~MHz',...
        'Far-Field~ ESPRIT~B=500~MHz'},...    
        'Interpreter','latex', 'Box','off');    
l3.FontSize = 12; 
ah3=axes('position',get(gca,'position'),'visible','off'); 
l4 = legend(ah3, [p33, p44],{  
        'Far-Field~ MUSIC~~B=10~GHz',...
        'Far-Field~ ESPRIT~B=10~GHz'},...    
        'Interpreter','latex', 'Box','off');    
l4.FontSize = 12; 