%%
clear
clc
set(0,'defaultfigurecolor','w') 
figure; hold on; box on; grid on;
set(gca,'FontName','Times New Roman','FontSize',12);
xlabel('P/dBm');
ylabel('NMSE/dB');
ylim([-55 -12.5])
%% My Algorithm
load('.\DataFinal\CE_OMP_BS_UE\far_LA_PgOMP.mat')
p1 = plot(P_t_dBm+15, 10*log10(NMSE(:,1)),'m-o', 'LineWidth', 1.5 ,'MarkerSize',10);
p2 = plot(P_t_dBm+15, 10*log10(NMSE(:,2)),'m:s', 'LineWidth', 1.5 ,'MarkerSize',10);
load('.\DataFinal\CE_OMP_BS_UE\far_PgOMP.mat')
p3 = plot(P_t_dBm+15, 10*log10(NMSE(:,1)),'b:o', 'LineWidth', 1.5 ,'MarkerSize',10);
p4 = plot(P_t_dBm+15, 10*log10(NMSE(:,2)),'b:s', 'LineWidth', 1.5 ,'MarkerSize',10);
load('.\DataFinal\CE_OMP_BS_UE\far_LA_PgOMP_PathInCluster_1.mat')
p5 = plot(P_t_dBm+15, 10*log10(NMSE(:,1)),'g:o', 'LineWidth', 1.5 ,'MarkerSize',10);
p6 = plot(P_t_dBm+15, 10*log10(NMSE(:,2)),'g:s', 'LineWidth', 1.5 ,'MarkerSize',10);
load('.\DataFinal\CE_OMP_BS_UE\far_PgOMP_PathInCluster_1.mat')
p7 = plot(P_t_dBm+15, 10*log10(NMSE(:,1)),'r:o', 'LineWidth', 1.5 ,'MarkerSize',10);
p8 = plot(P_t_dBm+15, 10*log10(NMSE(:,2)),'r:s', 'LineWidth', 1.5 ,'MarkerSize',10);
%% Baseline
load('.\DataFinal\CE_OMP_BS_UE\far_LS.mat')
p9 = plot(P_t_dBm+15, 10*log10(NMSE_GenieLS),'c:o', 'LineWidth', 1.5 ,'MarkerSize',10);
load('.\DataFinal\CE_OMP_BS_UE\far_LS_PathInCluster_1.mat')
p10 = plot(P_t_dBm+15, 10*log10(NMSE_GenieLS),'c:s', 'LineWidth', 1.5 ,'MarkerSize',10);
load('.\DataFinal\CE_OMP_BS_UE\far_PSOMP.mat')
p11 = plot(P_t_dBm+15, 10*log10(NMSE),'c:<', 'LineWidth', 1.5 ,'MarkerSize',10);
load('.\DataFinal\CE_OMP_BS_UE\far_PSOMP_PathInCluster_1.mat')
p12 = plot(P_t_dBm+15, 10*log10(NMSE),'c:>', 'LineWidth', 1.5 ,'MarkerSize',10);
%% legend    
l1 = legend([p1, p2, p3, p4, p5, p6, p7, p8],{  
        'LA-GMMV-OMP~$G_l=6$ $N_s=6$',...
        'LA-GMMV-OMP~$G_l=6$ $N_s=1$',...
        'GMMV-OMP~~~~~~$G_l=6$ $N_s=6$',...
        'GMMV-OMP~~~~~~$G_l=6$ $N_s=1$',...
        'LA-GMMV-OMP~$G_l=1$ $N_s=6$',...
        'LA-GMMV-OMP~$G_l=1$ $N_s=1$',...
        'GMMV-OMP~~~~~~$G_l=1$ $N_s=6$',...
        'GMMV-OMP~~~~~~$G_l=1$ $N_s=1$'},...
        'Interpreter','latex', 'Box','off');
l1.FontSize = 10;    
ah=axes('position',get(gca,'position'),'visible','off'); 

l2 = legend(ah, [p9, p10, p11, p12],{  
        'Genie-LS~$G_l=6$',...
        'Genie-LS~$G_l=1$',...
        'PSOMP~~$G_l=6$',...
        'PSOMP~~$G_l=1$'},...
        'Interpreter','latex', 'Box','off');
l2.FontSize = 10; 