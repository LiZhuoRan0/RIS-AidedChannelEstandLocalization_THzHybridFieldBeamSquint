clear
clc
main_far_LA_PgOMP
save('.\data\far_LA_PgOMP',...
    'P_t_dBm','NMSE', 'Times')
%%
clear
clc
main_far_PgOMP
save('.\data\far_PgOMP',...
    'P_t_dBm','NMSE', 'Times')
%%
clear
clc
main_far_LA_PgOMP_PathInCluster_1
save('.\data\far_LA_PgOMP_PathInCluster_1',...
    'P_t_dBm','NMSE', 'Times')
%%
clear
clc
main_far_PgOMP_PathInCluster_1
save('.\data\far_PgOMP_PathInCluster_1',...
    'P_t_dBm','NMSE', 'Times')
%%
clear; clc
main_far_GenieLS_PathInCluster_1
save('.\data\far_LS_PathInCluster_1',...
    'P_t_dBm','NMSE_GenieLS')
%%
clear; clc
main_far_GenieLS
save('.\data\far_LS',...
    'P_t_dBm','NMSE_GenieLS')
%%
clear; clc
main_far_PSOMP
save('.\data\far_PSOMP',...
    'P_t_dBm','NMSE', 'Times')
%%
clear; clc
main_far_PSOMP_PathInCluster_1
save('.\data\far_PSOMP_PathInCluster_1',...
    'P_t_dBm','NMSE', 'Times')