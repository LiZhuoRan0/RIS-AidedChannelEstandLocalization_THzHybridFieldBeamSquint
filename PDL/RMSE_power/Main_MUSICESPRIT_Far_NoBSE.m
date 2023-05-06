clear
clc
addpath('..')
fprintf('==============================================\n');
rng(666, 'twister');
%% parameter setting
N               =           256;                % number of BS antennas
NRIS            =           256;                % number of RIS elements
NRF             =           4;
fc              =           100e9;              %100GHz
lambda_c        =           3e8/fc;
d               =           lambda_c/2;
B               =           500e6;               %10GHz
M               =           2048;               %number of carrier frequency
K               =           -0.45/1000;% -0.45dB/km
S               =           sqrt(4*pi);
angleOverSmp    =           2;
% S_eff           =           2;
S_eff           =           (NRIS*d)^2;
AntGain         =           19.1+1.34;%dB  

P_BS                           = 8;% number observed time slots
P_RIS                           = 16;

Nit             =           2e1;
% Nit             =           10;
%% -    BS，RIS，UE, scatter location setting
% RIS location
RIS.x                          = 40*sqrt(2);
RIS.y                          = 0;

angle_RIS                    = -pi/4;% BS ULA 朝向
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
angle_BS                      = pi/4;% BS ULA orientation
BS.x1                         = BS.x - d*(N-1)/2*abs(sin(angle_BS));
BS.xN                         = BS.x + d*(N-1)/2*abs(sin(angle_BS));
BS.y1                         = BS.y - d*(N-1)/2*abs(sin(angle_BS));
BS.yN                         = BS.y + d*(N-1)/2*abs(sin(angle_BS));
% location of each BS antenna element
BS.xAll = linspace(BS.x1, BS.xN, N).';
BS.yAll = linspace(BS.y1, BS.yN, N).';
%% -    transmit power, noise power
P_noise = 10^(-17.4)*(B/2048)*(M);

P_t_dBm_All = (15:5:45);%dBmW
%%
RMSE.x_full                      = zeros(1, length(P_t_dBm_All));
RMSE.y_full                      = zeros(1, length(P_t_dBm_All));
RMSE.r_full                      = zeros(1, length(P_t_dBm_All));
RMSE.theta_full                  = zeros(1, length(P_t_dBm_All));

RMSE.x_full_ESPRIT                      = zeros(1, length(P_t_dBm_All));
RMSE.y_full_ESPRIT                      = zeros(1, length(P_t_dBm_All));
RMSE.r_full_ESPRIT                      = zeros(1, length(P_t_dBm_All));
RMSE.theta_full_ESPRIT                  = zeros(1, length(P_t_dBm_All));

RMSE.x_hybrid                      = zeros(1, length(P_t_dBm_All));
RMSE.y_hybrid                       = zeros(1, length(P_t_dBm_All));
RMSE.r_hybrid                       = zeros(1, length(P_t_dBm_All));
RMSE.theta_hybrid                  = zeros(1, length(P_t_dBm_All));
%% -    -   scatterers
num_scatter_UE2BS = 3;
num_scatter_UE2RIS= 3;
PathInCluster = 6;
SizeOfScatter = 1;
%% parameter generation
% UE location
UE.x                         = 40.2 ;
UE.y                         = -20.1;

r_UE2BS                     = sqrt((UE.x-BS.x).^2 +...
                                    (UE.y-BS.y).^2);   
r_UE2BS_All                 = sqrt((UE.x-BS.xAll).^2 +...
                            (UE.y-BS.yAll).^2);  
r_UE2RIS                    = sqrt((UE.x-RIS.x).^2 +...
                                    (UE.y-RIS.y).^2); 
r_UE2RIS_All                 = sqrt((UE.x-RIS.xAll).^2 +...
                            (UE.y-RIS.yAll).^2);  
r_BS2RIS                    = sqrt((BS.x - RIS.x).^2 + (BS.y - RIS.y).^2);
GainAnt                     = 10^(AntGain/10);%4*pi/(atan(NRIS*d/2/r_BS2RIS)*atan(NRIS*d/2/r_BS2RIS)*2);% THz Antenna Gain
theta_UE2BS                 = atan((UE.y-BS.y)/(UE.x-BS.x));         
theta_UE2RIS                = atan((UE.y-RIS.y)/(UE.x-RIS.x));
%% -    -   scatterers
Scatter_UE2BS = genScatterUEBS(K, S, N, BS, UE, num_scatter_UE2BS,...                        
                PathInCluster,...
                SizeOfScatter);
Scatter_UE2RIS = genScatterUEBS(K, S, N, RIS, UE, num_scatter_UE2RIS,...                        
                PathInCluster,...
                SizeOfScatter);
%% parameter output
fprintf('true UE2BS distance = %10.6f\n', r_UE2BS);                            
fprintf('true UE2RIS distance = %10.6f\n', r_UE2RIS);
fprintf('true UE2BS angle (rad) = %10.8f\n', theta_UE2BS); 
fprintf('true UE2BS angle (no unit) = %10.8f\n', sin(pi/2-angle_BS+theta_UE2BS)); 
fprintf('true UE2RIS angle (rad) = %10.8f\n', theta_UE2RIS); 
fprintf('true UE2RIS angle (no unit) = %10.8f\n', sin(pi/2+angle_RIS-theta_UE2RIS));  
fprintf('P_BS = %d\n', P_BS);
fprintf('P_RIS = %d\n', P_RIS);
fprintf('num_scatter_UE2BS = %d\n', num_scatter_UE2BS);
fprintf('num_scatter_UE2RIS = %d\n', num_scatter_UE2RIS); 
%% localization
H_UE2BS_absolute = NearFieldH(K, N, sin(pi/2-angle_BS+[theta_UE2BS, Scatter_UE2BS.theta])...
    , [r_UE2BS, Scatter_UE2BS.r], [r_UE2BS, Scatter_UE2BS.r_LastHop], [1./(r_UE2BS_All*4*pi).*(10.^(K.*r_UE2BS_All/10)) Scatter_UE2BS.alpha], fc, B, M, 'absolute');


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
    A_RIS_tmp = genA_RIS(sin(-angle_BS), sqrt((BS.x-RIS.x)^2+(BS.y-RIS.y)^2), fc, B, N, NRF, i_P, P_RIS)/sqrt(NRF)*sqrt(NRIS);
    Phi = diag( 2*randi([0 1], NRIS, 1) -1 );
    Y_RIS_absolute_tmp = zeros(N, M);
    %             Y_RIS_relative = zeros(N, M);
    for i_M = 1:M
        Y_RIS_absolute_tmp(:, i_M) = H_BR(:, :, i_M)*Phi*H_UE2RIS_absolute(:, i_M)*sqrt(S_eff);
        A_RIS_1(:, :, i_M) = A_RIS_tmp*H_BR(:, :, i_M)*Phi*sqrt(S_eff);
        %                   Y_RIS_relative(:, i_M) = H_BR(:, :, i_M)*Phi*H_UE2RIS_relative(:, i_M);
    end    
    Y_RIS_absolute = A_RIS_tmp*(Y_RIS_absolute_tmp*sqrt(GainAnt));
    Y_RIS((i_P-1)*NRF + (1:NRF),:) = Y_RIS_absolute;
    A_RIS((i_P-1)*NRF + (1:NRF), :, :) = A_RIS_1*sqrt(GainAnt);
end
for i_P_t_dBm = 1:length(P_t_dBm_All)   
    P_t_dBm = P_t_dBm_All(i_P_t_dBm);
    for i_it = 1:Nit
        fprintf('==============================================\n');   
        fprintf('P_t_dBm = %3d\n',P_t_dBm);
        fprintf('i_it/Nit = %3d/%3d\n', i_it, Nit);
        %% -    generate Y, A, H according to the parameters defined before
        %% -    -   UE2BS
        A = exp(1j*2*pi*rand(NRF*P_BS, N))/sqrt(NRF);   
        A_full = eye(N);
        A(1,:) = zeros(1, N);
        if mod(N,2) == 0%N is even
            A(1, ceil(N/2)) = 1/sqrt(NRF);
            A(1, ceil((N+1)/2)) = 1/sqrt(NRF);
        else
            A(1, ceil((N+1)/2)) = 1/sqrt(NRF);
        end
        

        Noise = (1/sqrt(2)*randn(size(A, 1), M) + 1j/sqrt(2)*randn(size(A, 1), M));
        Noise_Full = (1/sqrt(2)*randn(size(A_full, 1), M) + 1j/sqrt(2)*randn(size(A_full, 1), M));
%         Noise(1,:) = 0*Noise(1,:);
        Y_absolute_all_full = A_full*H_UE2BS_absolute + sqrt(P_noise)/sqrt(10^(P_t_dBm/10))*Noise_Full;
        Y_absolute_all = (A*H_UE2BS_absolute + sqrt(P_noise)/sqrt(10^(P_t_dBm/10))/sqrt(10^(AntGain/10))*Noise);
        r_est = MUSIC_noEVD(M, Y_absolute_all, fc, B, 'B');
        fprintf('UE2BS distance estimated by MUSIC = %10.6f\n', r_est);
%% -    -   UE2RIS2BS
        Noise = 1/sqrt(2)*randn(size(Y_RIS)) + 1j/sqrt(2)*randn(size(Y_RIS));
        Y_RIS_Noise = Y_RIS + sqrt(P_noise)/sqrt(10^(P_t_dBm/10))*Noise;
        %             10*log10(10.^(P_t_dBm/10).*norm(Y_RIS,'fro')^2./norm(sqrt(P_noise)*Noise,'fro')^2)
        r_UE2RIS_est = MUSIC_noEVD(M, Y_RIS_Noise, fc, B, 'R');
        r_RIS_est = r_UE2RIS_est - sqrt((BS.x - RIS.x).^2 + (BS.y - RIS.y).^2)  + 3e8/B*M;  
        fprintf('UE2RIS distance estimated by MUSIC = %10.6f\n', r_RIS_est);    
        %% Obain hyperbola through 2 delays         
        a = (r_est-r_RIS_est)/2;
        c = (RIS.x - BS.x)/2;
        b = sqrt(c*c-a*a);
        
%% - Estimate the AoA
        
        theta_UE2BS_cor_est_full = MUSIC_BS(Y_absolute_all_full, A_full);
        theta_UE2BS_cor_est_hybrid = MUSIC_BS(Y_absolute_all, A);
        psi = TLS_ESPRIT_Algorithm(Y_absolute_all_full, 1);
        theta_UE2BS_cor_est_full_ESPRIT = log(psi)/(1j*pi);
       %% - localization through the AoA and hyperbola
        k_B2U_MUSIC_full = tan(asin(theta_UE2BS_cor_est_full)+angle_BS-pi/2);
        k_B2U_MUSIC_hybrid = tan(asin(theta_UE2BS_cor_est_hybrid)+angle_BS-pi/2);
        k_B2U_ESPRIT = tan(asin(theta_UE2BS_cor_est_full_ESPRIT)+angle_BS-pi/2);
        
        m_1_MUSIC_full = 1-(a/b.*k_B2U_MUSIC_full).^2;
        m_1_MUSIC_hybrid = 1-(a/b.*k_B2U_MUSIC_hybrid).^2;
        m_1_ESPRIT = 1-(a/b.*k_B2U_ESPRIT).^2;              
        
        m_2 = -2*c;
        m_3 = b*b;
        
        UE.x_1_angle_MUSIC_full = (-m_2+sqrt(m_2.^2-4.*m_1_MUSIC_full*m_3))./(2.*m_1_MUSIC_full);
        UE.y_1_angle_MUSIC_full = UE.x_1_angle_MUSIC_full.*k_B2U_MUSIC_full;
        UE.x_1_angle_MUSIC_hybrid = (-m_2+sqrt(m_2.^2-4.*m_1_MUSIC_hybrid*m_3))./(2.*m_1_MUSIC_hybrid);
        UE.y_1_angle_MUSIC_hybrid = UE.x_1_angle_MUSIC_hybrid.*k_B2U_MUSIC_hybrid; 
        UE.x_1_angle_ESPRIT = (-m_2+sqrt(m_2.^2-4.*m_1_ESPRIT*m_3))./(2.*m_1_ESPRIT);
        UE.y_1_angle_ESPRIT = UE.x_1_angle_ESPRIT.*k_B2U_ESPRIT; 

        %% Output
        theta_UE2BS_final_full = theta_UE2BS_cor_est_full;
        theta_UE2BS_final_full_ESPRIT = theta_UE2BS_cor_est_full_ESPRIT;
        theta_UE2BS_final_hybrid = theta_UE2BS_cor_est_hybrid;

        fprintf('True UE2BS angle (in unit 1) = %10.8f\n', sin(pi/2-angle_BS+atan((UE.y-BS.y)/(UE.x-BS.x))));    
        fprintf('UE2BS angle estimated by the FullDigtal MUSIC (in unit 1) = %10.8f\n', theta_UE2BS_final_full);    
        fprintf('UE2BS angle estimated by the FullDigtal ESPRIT (in unit 1) = %10.8f\n', theta_UE2BS_final_full_ESPRIT);    
        fprintf('UE2BS angle estimated by the Hybrid MUSIC (in unit 1) = %10.8f\n', theta_UE2BS_final_hybrid);
        %% calculate the angle RMSE
        RMSE.theta_full(i_P_t_dBm) = RMSE.theta_full(i_P_t_dBm) + norm(asin(theta_UE2BS_final_full)+angle_BS-pi/2 - theta_UE2BS).^2;
        RMSE.theta_full_ESPRIT(i_P_t_dBm) = RMSE.theta_full_ESPRIT(i_P_t_dBm) + norm(asin(theta_UE2BS_final_full_ESPRIT)+angle_BS-pi/2 - theta_UE2BS).^2;
        RMSE.theta_hybrid(i_P_t_dBm) = RMSE.theta_hybrid(i_P_t_dBm) + norm(asin(theta_UE2BS_final_hybrid)+angle_BS-pi/2 - theta_UE2BS).^2;
        %% calculate the RMSE of x,y,r
        RMSE.x_full(i_P_t_dBm) = RMSE.x_full(i_P_t_dBm) + norm(UE.x_1_angle_MUSIC_full - UE.x)^2;
        RMSE.y_full(i_P_t_dBm) = RMSE.y_full(i_P_t_dBm) + norm(UE.y_1_angle_MUSIC_full - UE.y)^2;
        RMSE.r_full(i_P_t_dBm) = RMSE.r_full(i_P_t_dBm) + norm(sqrt(UE.x^2 + UE.y^2) - ...
                        sqrt(UE.x_1_angle_MUSIC_full^2 + UE.y_1_angle_MUSIC_full^2))^2;

        RMSE.x_full_ESPRIT(i_P_t_dBm) = RMSE.x_full_ESPRIT(i_P_t_dBm) + norm(UE.x_1_angle_ESPRIT - UE.x)^2;
        RMSE.y_full_ESPRIT(i_P_t_dBm) = RMSE.y_full_ESPRIT(i_P_t_dBm) + norm(UE.y_1_angle_ESPRIT - UE.y)^2;
        RMSE.r_full_ESPRIT(i_P_t_dBm) = RMSE.r_full_ESPRIT(i_P_t_dBm) + norm(sqrt(UE.x^2 + UE.y^2) - ...
                        sqrt(UE.x_1_angle_ESPRIT^2 + UE.y_1_angle_ESPRIT^2))^2;
                    
        RMSE.x_hybrid(i_P_t_dBm) = RMSE.x_hybrid(i_P_t_dBm) + norm(UE.x_1_angle_MUSIC_hybrid - UE.x)^2;
        RMSE.y_hybrid(i_P_t_dBm) = RMSE.y_hybrid(i_P_t_dBm) + norm(UE.y_1_angle_MUSIC_hybrid - UE.y)^2;
        RMSE.r_hybrid(i_P_t_dBm) = RMSE.r_hybrid(i_P_t_dBm) + norm(sqrt(UE.x^2 + UE.y^2) - ...
                        sqrt(UE.x_1_angle_MUSIC_hybrid^2 + UE.y_1_angle_MUSIC_hybrid^2))^2;
        
        fprintf('RMSE.x_full = %10.8f\n', sqrt(RMSE.x_full(1, i_P_t_dBm)/i_it));
        fprintf('RMSE.y_full = %10.8f\n', sqrt(RMSE.y_full(1, i_P_t_dBm)/i_it));
        fprintf('RMSE.r_full = %10.8f\n', sqrt(RMSE.r_full(1, i_P_t_dBm)/i_it));
        fprintf('RMSE.theta_full = %10.8f\n', sqrt(RMSE.theta_full(1, i_P_t_dBm)/i_it));
        
        fprintf('RMSE.x_full_ESPRIT = %10.8f\n', sqrt(RMSE.x_full_ESPRIT(1, i_P_t_dBm)/i_it));
        fprintf('RMSE.y_full_ESPRIT = %10.8f\n', sqrt(RMSE.y_full_ESPRIT(1, i_P_t_dBm)/i_it));
        fprintf('RMSE.r_full_ESPRIT = %10.8f\n', sqrt(RMSE.r_full_ESPRIT(1, i_P_t_dBm)/i_it));
        fprintf('RMSE.theta_full_ESPRIT = %10.8f\n', sqrt(RMSE.theta_full_ESPRIT(1, i_P_t_dBm)/i_it));        
                    
        fprintf('RMSE.x_hybrid = %10.8f\n', sqrt(RMSE.x_hybrid(1, i_P_t_dBm)/i_it));
        fprintf('RMSE.y_hybrid = %10.8f\n', sqrt(RMSE.y_hybrid(1, i_P_t_dBm)/i_it));
        fprintf('RMSE.r_hybrid = %10.8f\n', sqrt(RMSE.r_hybrid(1, i_P_t_dBm)/i_it));
        fprintf('RMSE.theta_hybrid = %10.8f\n', sqrt(RMSE.theta_hybrid(1, i_P_t_dBm)/i_it));        
    end
end
RMSE.x_full = sqrt(RMSE.x_full/Nit);
RMSE.y_full = sqrt(RMSE.y_full/Nit);
RMSE.r_full = sqrt(RMSE.r_full/Nit);
RMSE.theta_full = sqrt(RMSE.theta_full/Nit);

RMSE.x_full_ESPRIT = sqrt(RMSE.x_full_ESPRIT/Nit);
RMSE.y_full_ESPRIT = sqrt(RMSE.y_full_ESPRIT/Nit);
RMSE.r_full_ESPRIT = sqrt(RMSE.r_full_ESPRIT/Nit);
RMSE.theta_full_ESPRIT = sqrt(RMSE.theta_full_ESPRIT/Nit);

RMSE.x_hybrid = sqrt(RMSE.x_hybrid/Nit);
RMSE.y_hybrid = sqrt(RMSE.y_hybrid/Nit);
RMSE.r_hybrid = sqrt(RMSE.r_hybrid/Nit);
RMSE.theta_hybrid = sqrt(RMSE.theta_hybrid/Nit);
save('RMSE_Pt_Far', 'RMSE', 'P_t_dBm')