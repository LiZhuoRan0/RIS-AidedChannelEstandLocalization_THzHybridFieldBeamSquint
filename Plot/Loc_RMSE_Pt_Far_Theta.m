%% theta
clear
clc
set(0,'defaultfigurecolor','w') 
figure; hold on; box on; grid on;
set(gca,'FontName','Times New Roman','FontSize',12);
xlabel('P/dBm');
ylabel('RMSE_{\vartheta}/rad');
set(gca, 'YScale', 'log');
ylim([1e-6 1e0])
P_t_dBm_All = [15 20 25 30 35 40 45];
%% -    OMP
load('.\DataFinal\Loc_RMSE_Pt\RMSE_Pt_far_CDL.mat')
p1 = plot(P_t_dBm_All(1:2:end), RMSE.theta_co(1:2:end),'r:s', 'LineWidth', 1.5);
p2 = plot(P_t_dBm_All(1:2:end), RMSE.theta(1:2:end),'g-o', 'LineWidth', 1.5); 
%% -    TDoA
load('.\DataFinal\Loc_RMSE_Pt\RMSE_Pt_PDL_Far.mat')
p3 = plot(P_t_dBm_All(1:2:end), RMSE.theta(1:2:end),'m:o', 'LineWidth', 1.5, 'MarkerSize',10 );
load('.\DataFinal\Loc_RMSE_Pt\RMSE_Pt_MUSICESPRIT_Far.mat')
p4 = plot(P_t_dBm_All(1:2:end), RMSE.theta_full(1:2:end),'b:o', 'LineWidth', 1.5,'MarkerSize',10);
p5 = plot(P_t_dBm_All(1:2:end), RMSE.theta_full_ESPRIT(1:2:end),'b:x', 'LineWidth', 1.5,'MarkerSize',10);
p6 = plot(P_t_dBm_All(1:2:end), RMSE.theta_hybrid(1:2:end),'b:s', 'LineWidth', 1.5,'MarkerSize',10);
%% -    legend
l1 = legend([p1, p2],{  
        'GMMV-OMP',...
        'CDL'},...
        'Interpreter','latex', 'Box','off');
l1.FontSize = 12;    
ah=axes('position',get(gca,'position'),'visible','off');    
l2 = legend(ah, [p3, p4, p5, p6],{  
        'PDL',...    
        'Full-digital MUSIC',...
        'Full-digital ESPRIT',...
        'Hybrid MUSIC'},...    
        'Interpreter','latex', 'Box','off');    
l2.FontSize = 12; 