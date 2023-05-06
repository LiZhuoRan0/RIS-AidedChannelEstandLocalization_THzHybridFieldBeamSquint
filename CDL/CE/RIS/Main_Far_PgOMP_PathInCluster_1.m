clear
clc
addpath('.\function')
fprintf('==============================================\n');
rng(666,"twister");
%% parameter setting
N               =           256;                % number of BS antennas
NRIS            =           256;                % number of RIS elements
NRF             =           4;
fc              =           100e9;              %100GHz
lambda_c        =           3e8/fc;
d               =           lambda_c/2;
B               =           10e9;               %10GHz
M               =           64;                 %number of carrier frequency
K               =           -0.45/1000;% -0.45dB/km
S               =           sqrt(4*pi);
angleOverSmp    =           2;
S_eff           =           (NRIS*d)^2;%RIS scattering area
AntGain         =           19.1+1.34;%dB 
L_2             =           [6 1];
% L_2             =           [6];
threshold       =           [0.85 0.98];

P_RIS                           = 32;% number observed time slots

Nit             =           2e1;
% Nit             =           5;
L_max           =           20;% maximum iteration number of the OMP
%% -    BS£¨RIS£¨UE, scatter location setting
% RISŒª÷√
RIS.x                          = 40*sqrt(2);
RIS.y                          = 0;

angle_RIS                    = -pi/4;% BS ULA orientation
RIS.x1                         = RIS.x - d*(NRIS-1)/2*abs(sin(angle_RIS));
RIS.y1                         = RIS.y + d*(NRIS-1)/2*abs(sin(angle_RIS));
RIS.xN                         = RIS.x + d*(NRIS-1)/2*abs(sin(angle_RIS));
RIS.yN                         = RIS.y - d*(NRIS-1)/2*abs(sin(angle_RIS));
% location of each RIS element
RIS.xAll = linspace(RIS.x1, RIS.xN, N).';
RIS.yAll = linspace(RIS.y1, RIS.yN, N).';

% BS location
BS.x                          = 0;% x-coordinate of the BS center
BS.y                          = 0;
angle_BS                    = pi/4;% BS ULA orientation
BS.x1                         = BS.x - d*(N-1)/2*abs(sin(angle_BS));
BS.xN                         = BS.x + d*(N-1)/2*abs(sin(angle_BS));
BS.y1                         = BS.y - d*(N-1)/2*abs(sin(angle_BS));
BS.yN                         = BS.y + d*(N-1)/2*abs(sin(angle_BS));
% location of each BS antenna element
BS.xAll = linspace(BS.x1, BS.xN, N).';
BS.yAll = linspace(BS.y1, BS.yN, N).';

num_scatter_UE2RIS= 3;
PathInCluster = 1;
SizeOfScatter = 1;

%% -    transmit power, noise power
P_noise = 10^(-17.4)*(B/2048)*(M);

P_t_dBm = (15:5:45);%dBmW

NMSE = zeros(length(P_t_dBm), length(L_2));
Times = zeros(length(P_t_dBm), length(L_2));

%% parameter generation
% UE location
UE.x                         = 40.2 ;
UE.y                         = -20.2;
 
r_UE2RIS                    = sqrt((UE.x-RIS.x).^2 +...
                                    (UE.y-RIS.y).^2); 
r_UE2RIS_All                 = sqrt((UE.x-RIS.xAll).^2 +...
                            (UE.y-RIS.yAll).^2);  
r_BS2RIS                    = sqrt((BS.x - RIS.x).^2 + (BS.y - RIS.y).^2);
GainAnt                     = 10^(AntGain/10);%4*pi/(atan(NRIS*d/2/r_BS2RIS)*atan(NRIS*d/2/r_BS2RIS)*2);% THz Antenna Gain                  
theta_UE2RIS                = atan((UE.y-RIS.y)/(UE.x-RIS.x));
%% -    -   scatterers
Scatter_UE2RIS = genScatterUEBS(K, S, N, RIS, UE, num_scatter_UE2RIS,...                        
                PathInCluster,...
                SizeOfScatter);
%% parameter output                         
fprintf('true UE2RIS distance = %10.6f\n', r_UE2RIS);
fprintf('true UE2RIS angle (rad) = %10.8f\n', theta_UE2RIS); 
fprintf('true UE2RIS angle (no unit) = %10.8f\n', sin(pi/2+angle_RIS-theta_UE2RIS));  
fprintf('P_RIS = %d\n', P_RIS);
fprintf('num_scatter_UE2RIS = %d\n', num_scatter_UE2RIS);
%% localization
for i_P_t_dBm = 1:length(P_t_dBm)
    for i_it = 1:Nit
        fprintf('==============================================\n');                  
        fprintf('P_t_dBm = %3d\n',P_t_dBm(i_P_t_dBm));
        fprintf('i_it/Nit = %3d/%3d\n', i_it, Nit);
%% -    generate Y, A, H according to the parameters defined before
%% -    -   UE2RIS2BS
        H_UE2RIS_absolute = NearFieldH(K, NRIS, sin(pi/2+angle_RIS-[theta_UE2RIS, Scatter_UE2RIS.theta])...            
                    , [r_UE2RIS, Scatter_UE2RIS.r], [r_UE2RIS, Scatter_UE2RIS.r_LastHop], [1./(r_UE2RIS_All*4*pi).*(10.^(K.*r_UE2RIS_All/10)), Scatter_UE2RIS.alpha], fc, B, M, 'absolute');
        H_BR = zeros(N, NRIS, M);
        for i_M = 1:M
           H_BR(:, :, i_M) = genRIS2BS_H(RIS.x1, RIS.y1,...
                                RIS.xN, RIS.yN,...
                                BS.x1, BS.y1,...
                                BS.xN, BS.yN,...
                                NRIS, N, fc-B/2+(B/M)*(i_M-1), K); %fc-B/2+(B/M)*(i_M-1)
        end
        
        Y_RIS = zeros(P_RIS*NRF, M);
        A_RIS = zeros(P_RIS*NRF, NRIS, M);
        A_RIS_1 = zeros(NRF, NRIS, M);        
       
        for i_P = 1:P_RIS
            A_RIS_tmp = genA_RIS(sin(-angle_BS), sqrt((BS.x-RIS.x)^2+(BS.y-RIS.y)^2), fc, B, N, NRF, i_P, P_RIS)/sqrt(NRF)*sqrt(N);
            Phi = diag( 2*randi([0 1], NRIS, 1) -1 ); 
            Y_RIS_absolute_tmp = zeros(N, M);
            for i_M = 1:M
                  Y_RIS_absolute_tmp(:, i_M) = H_BR(:, :, i_M)*Phi*H_UE2RIS_absolute(:, i_M)*sqrt(S_eff);
                  A_RIS_1(:, :, i_M) = A_RIS_tmp*H_BR(:, :, i_M)*Phi*sqrt(S_eff);
            end
            Noise = (1/sqrt(2)*randn(NRF, M) + 1j/sqrt(2)*randn(NRF, M));
            Y_RIS_absolute = A_RIS_tmp*(Y_RIS_absolute_tmp*sqrt(GainAnt))+ sqrt(P_noise)/sqrt(10^(P_t_dBm(i_P_t_dBm)/10))*Noise;
            Y_RIS((i_P-1)*NRF + (1:NRF),:) = Y_RIS_absolute;
            A_RIS((i_P-1)*NRF + (1:NRF), :, :) = A_RIS_1*sqrt(GainAnt);
        end 
%% channel estimation and calculate NMSE
        for i_L = 1:length(L_2)
            [H_hat, it_num] = PgOMP_CE_RIS(L_2(i_L), threshold(i_L), Y_RIS, A_RIS, L_max,...
                       NRIS, fc, NRF, P_RIS, angleOverSmp, M, B);
            NMSE(i_P_t_dBm, i_L) = NMSE(i_P_t_dBm, i_L) + (norm(H_UE2RIS_absolute - H_hat, 'fro')/norm(H_UE2RIS_absolute, 'fro')).^2;
            Times(i_P_t_dBm, i_L) = Times(i_P_t_dBm, i_L) + it_num;
            fprintf("L_2 = %d\n", L_2(i_L));
            fprintf('Average NMSE=%ddB; Current NMSE=%ddB\n',10*log10(NMSE(i_P_t_dBm, i_L)/i_it), 20*log10(norm(H_UE2RIS_absolute - H_hat, 'fro')/norm(H_UE2RIS_absolute, 'fro')));
            fprintf("Average Itetation Number=%d; Current Iteration Number =%d\n", Times(i_P_t_dBm, i_L)/i_it, it_num);
        end
    end
    NMSE(i_P_t_dBm,:) = NMSE(i_P_t_dBm,:)/Nit;
    Times(i_P_t_dBm, :) = Times(i_P_t_dBm, :)/Nit;
end