%% Add Path
addpath('.\CDF')
% addpath('.\RMSE_power')
% addpath('.\RMSE_P')
%% CDF
clear; clc
main_Base
save('.\Data\CDF\main_Base_CDL', 'P_t_dBm', 'RMSE')

clear; clc
main_UEDistance
save('.\Data\CDF\main_UEDistance_CDL', 'P_t_dBm', 'UEx_All', 'UEy_All', 'RMSE')

clear; clc
main_Bandwidth
save('.\Data\CDF\main_Bandwidth_CDL', 'P_t_dBm', 'B_All', 'RMSE')

clear; clc
main_NumBS
save('.\Data\CDF\main_NumBS_CDL_512', 'P_t_dBm', 'NAll', 'RMSE')

clear; clc
main_NumRIS
save('.\Data\CDF\main_NumRIS_CDL', 'P_t_dBm', 'NRIS_All', 'RMSE')
%% RMSE-Pt
clear; clc
main_Armijo
save('.\Data\Pt\RMSE_Pt_near_CDL', 'P_t_dBm', 'RMSE')

clear; clc
main_farfield_Armijo
save('.\Data\Pt\RMSE_Pt_far_CDL', 'P_t_dBm', 'RMSE')
%% RMSE-P
% clear; clc
% main
% save('.\Data\P\RMSE_P_CDL_CDF_32', 'P_BS_All', 'P_RIS_All', 'RMSE')