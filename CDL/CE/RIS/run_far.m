clear
clc
Main_Far_LA_PgOMP
save('.\data\far_LA_PgOMP',...
    'P_t_dBm','NMSE', 'Times')
%%
clear
clc
Main_Far_PgOMP
save('.\data\far_PgOMP',...
    'P_t_dBm','NMSE', 'Times')
%%
clear
clc
Main_Far_LA_PgOMP_PathInCluster_1
save('.\data\far_LA_PgOMP_PathInCluster_1',...
    'P_t_dBm','NMSE', 'Times')
%%
clear
clc
Main_Far_PgOMP_PathInCluster_1
save('.\data\far_PgOMP_PathInCluster_1',...
    'P_t_dBm','NMSE', 'Times')
%%
clear; clc
Main_Far_GenieLS_PathInCluster_1
save('.\data\far_LS_PathInCluster_1',...
    'P_t_dBm','NMSE_GenieLS')
%%
clear; clc
Main_Far_GenieLS
save('.\data\far_LS',...
    'P_t_dBm','NMSE_GenieLS')
%%
clear; clc
Main_Far_PSOMP
save('.\data\far_PSOMP',...
    'P_t_dBm','NMSE', 'Times')
%% 
clear; clc
Main_Far_PSOMP_PathInCluster_1
save('.\data\far_PSOMP_PathInCluster_1',...
    'P_t_dBm','NMSE', 'Times')