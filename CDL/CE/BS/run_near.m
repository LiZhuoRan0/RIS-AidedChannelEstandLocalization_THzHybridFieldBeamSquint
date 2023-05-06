clear
clc
main_near_LA_PgOMP
save('.\data\near_LA_PgOMP',...
    'P_t_dBm','NMSE', 'Times')
%%
clear
clc
main_near_PgOMP
save('.\data\near_PgOMP',...
    'P_t_dBm','NMSE', 'Times')
%%
clear
clc
main_near_LA_PgOMP_PathInCluster_1
save('.\data\near_LA_PgOMP_PathInCluster_1',...
    'P_t_dBm','NMSE', 'Times')
%%
clear
clc
main_near_PgOMP_PathInCluster_1
save('.\data\near_PgOMP_PathInCluster_1',...
    'P_t_dBm','NMSE', 'Times')
%%
clear; clc
main_near_GenieLS_PathInCluster_1
save('.\data\near_LS_PathInCluster_1',...
    'P_t_dBm','NMSE_GenieLS')
%%
clear; clc
main_near_GenieLS
save('.\data\near_LS',...
    'P_t_dBm','NMSE_GenieLS')
%%
clear; clc
main_near_PSOMP
save('.\data\near_PSOMP',...
    'P_t_dBm','NMSE', 'Times')
%%
clear; clc
main_near_PSOMP_PathInCluster_1
save('.\data\near_PSOMP_PathInCluster_1',...
    'P_t_dBm','NMSE', 'Times')