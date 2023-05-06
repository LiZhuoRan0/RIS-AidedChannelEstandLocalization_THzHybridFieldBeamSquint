clear
clc
Main_Near_LA_PgOMP
save('.\data\near_LA_PgOMP_New',...
    'P_t_dBm','NMSE', 'Times')
%%
clear
clc
Main_Near_PgOMP
save('.\data\near_PgOMP_New',...
    'P_t_dBm','NMSE', 'Times')
%%
clear
clc
Main_Near_LA_PgOMP_PathInCluster_1
save('.\data\near_LA_PgOMP_PathInCluster_1',...
    'P_t_dBm','NMSE', 'Times')
%%
clear
clc
Main_Near_PgOMP_PathInCluster_1
save('.\data\near_PgOMP_PathInCluster_1',...
    'P_t_dBm','NMSE', 'Times')
%%
clear; clc
Main_Near_GenieLS_PathInCluster_1
save('.\data\near_LS_PathInCluster_1',...
    'P_t_dBm','NMSE_GenieLS')
%%
clear; clc
Main_Near_GenieLS
save('.\data\near_LS',...
    'P_t_dBm','NMSE_GenieLS')
%%
clear; clc
Main_Near_PSOMP
save('.\data\near_PSOMP',...
    'P_t_dBm','NMSE', 'Times')