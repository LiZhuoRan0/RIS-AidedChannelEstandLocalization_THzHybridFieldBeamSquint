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
M               =           2048;                 %number of carrier frequency
K               =           -0.45/1000;% -0.45dB/km
S               =           sqrt(4*pi);
angleOverSmp    =           2;
% S_eff           =           2;
S_eff           =           (NRIS*d)^2;
AntGain         =           19.1+1.34;%dB 

P_BS                           = 8;% number observed time slots
P_RIS                           = 16;

Nit             =           2e2;
% Nit             =           10;
%% -    BS，RIS，UE, scatter location setting
% RIS location
RIS.x                          = 20*sqrt(2);
RIS.y                          = 0;

angle_RIS                    = -pi/4;% BS ULA orientation
RIS.x1                         = RIS.x - d*(NRIS-1)/2*abs(sin(angle_RIS));
RIS.y1                         = RIS.y + d*(NRIS-1)/2*abs(sin(angle_RIS));
RIS.xN                         = RIS.x + d*(NRIS-1)/2*abs(sin(angle_RIS));
RIS.yN                         = RIS.y - d*(NRIS-1)/2*abs(sin(angle_RIS));
% location of each RIS element
RIS.xAll = linspace(RIS.x1, RIS.xN, NRIS).';
RIS.yAll = linspace(RIS.y1, RIS.yN, NRIS).';

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
P_noise = 10^(-17.4)*(10e9/2048)*(M);

P_t_dBm = 30;%dBmW
%% -    -   scatterers
num_scatter_UE2BS = 3;
num_scatter_UE2RIS= 3;
PathInCluster = 6;
SizeOfScatter = 1;
%% parameter generation
% UE location
UE.x                         = 20.1 ;
UE.y                         = -10.1;

r_UE2BS                     = sqrt((UE.x-BS.x).^2 +...
                                    (UE.y-BS.y).^2);   
r_UE2BS_All                 = sqrt((UE.x-BS.xAll).^2 +...
                            (UE.y-BS.yAll).^2);  
r_UE2RIS                    = sqrt((UE.x-RIS.x).^2 +...
                                    (UE.y-RIS.y).^2); 
r_UE2RIS_All                 = sqrt((UE.x-RIS.xAll).^2 +...
                            (UE.y-RIS.yAll).^2);  
r_BS2RIS                    = sqrt((BS.x - RIS.x).^2 + (BS.y - RIS.y).^2);
GainAnt                     = 10^(AntGain/10);%4*pi/(atan(256*d/2/r_BS2RIS)*atan(256*d/2/r_BS2RIS)*2);% THz Antenna Gain
theta_UE2BS                 = atan((UE.y-BS.y)/(UE.x-BS.x));         
theta_UE2RIS                = atan((UE.y-RIS.y)/(UE.x-RIS.x));
%% -    -   scatterers
Scatter_UE2BS = genScatterUEBS(K, S, N, BS, UE, num_scatter_UE2BS,...                        
                PathInCluster,...
                SizeOfScatter);
Scatter_UE2RIS = genScatterUEBS(K, S, NRIS, RIS, UE, num_scatter_UE2RIS,...                        
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
B_All               =           [100e6 500e6 1e9 5e9];               %10GHz
RMSE.x                      = zeros(length(B_All), Nit);
RMSE.y                      = zeros(length(B_All), Nit);
RMSE.r                      = zeros(length(B_All), Nit);
RMSE.theta                  = zeros(length(B_All), Nit);
%%
for i_B = 1:length(B_All)
    B = B_All(i_B);
%     P_noise = 10^(-17.4)*(B);

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
    for i_it = 1:Nit
        fprintf('==============================================\n');
        fprintf('P_t_dBm = %3d\n',P_t_dBm);        
        fprintf('B = %3d 1e9\n',B/1e9);
        fprintf('i_it/Nit = %3d/%3d\n', i_it, Nit);
        %% -    generate Y, A, H according to the parameters defined before
        %% -    -   UE2BS
        A = exp(1j*2*pi*rand(NRF*P_BS, N))/sqrt(NRF);
        A(1,:) = zeros(1, N);
        if mod(N,2) == 0%N is even
            A(1, ceil(N/2)) = 1/sqrt(NRF);
            A(1, ceil((N+1)/2)) = 1/sqrt(NRF);
        else
            A(1, ceil((N+1)/2)) = 1/sqrt(NRF);
        end
    
    
        Noise = (1/sqrt(2)*randn(size(A, 1), M) + 1j/sqrt(2)*randn(size(A, 1), M));
        %         Noise(1,:) = 0*Noise(1,:);
        Y_absolute_all = (A*H_UE2BS_absolute + sqrt(P_noise)/sqrt(10^(P_t_dBm/10))/sqrt(10^(AntGain/10))*Noise);
        r_est = MUSIC_noEVD(M, Y_absolute_all, fc, B, 'B');
        fprintf('UE2BS distance estimated by MUSIC = %10.6f\n', r_est);        
        %% -    -   -   convert absolute phase to the relative phase
        Y_relative_all = zeros(NRF*P_BS, M);
        for i = 2:NRF*P_BS
            Y_relative_all(i,:) = Y_absolute_all(i,:)./Y_absolute_all(1,:);
            Y_relative_all(i,:) = Y_relative_all(i,:)...
                /norm(Y_relative_all(i,:), 'fro')...
                *norm(Y_absolute_all(i,:), 'fro');
        end
        Y_relative_all(1,:) = Y_absolute_all(1,:);
        %% -    -   UE2RIS2BS        
        Noise = 1/sqrt(2)*randn(size(Y_RIS)) + 1j/sqrt(2)*randn(size(Y_RIS));
        Y_RIS_Noise = Y_RIS + sqrt(P_noise)/sqrt(10^(P_t_dBm/10))*Noise;
        %             10*log10(10.^(P_t_dBm/10).*norm(Y_RIS,'fro')^2./norm(sqrt(P_noise)*Noise,'fro')^2)
        r_UE2RIS_est = MUSIC_noEVD(M, Y_RIS_Noise, fc, B, 'R');
        r_RIS_est = r_UE2RIS_est - sqrt((BS.x - RIS.x).^2 + (BS.y - RIS.y).^2);
        fprintf('UE2RIS distance estimated by MUSIC = %10.6f\n', r_RIS_est);    
        %% Obain hyperbola through 2 delays
        a = (r_est-r_RIS_est)/2;
        c = (RIS.x - BS.x)/2;
        b = sqrt(c*c-a*a);
    
        theta_B2U_all = linspace(-pi/3, -pi/16, 200);
        k_B2U_all = tan(theta_B2U_all);
    
        m_1 = 1-(a/b.*k_B2U_all).^2;
        m_2 = -2*c;
        m_3 = b*b;
        UE.x_est = (-m_2+sqrt(m_2.^2-4.*m_1*m_3))./(2.*m_1);
        UE.y_est = UE.x_est.*k_B2U_all;        
        %% Obtain the coarse angle estimation on the hyperbola through correlation
        r_B2U_h = sqrt(UE.x_est.*UE.x_est + UE.y_est.*UE.y_est);
        theta_B2U_h = atan(UE.y_est./UE.x_est);
    
        r_est_All = sqrt((UE.x_est-BS.xAll).^2 +...
            (UE.y_est-BS.yAll).^2);
        theta_UE2BS_cor_est = find_max_cor_theta(K, r_est_All, r_B2U_h, theta_B2U_h, Y_absolute_all, A, B, N, M, fc, angle_BS);
        %% localization through the AoA and hyperbola
        k_B2U = tan(asin(theta_UE2BS_cor_est)+angle_BS-pi/2);
        m_1 = 1-(a/b.*k_B2U).^2;
        m_2 = -2*c;
        m_3 = b*b;
        UE.x_1_angle = (-m_2+sqrt(m_2.^2-4.*m_1*m_3))./(2.*m_1);
        UE.y_1_angle = UE.x_1_angle.*k_B2U;
    
        %         r_est_2 = sqrt(UE.x_1_angle^2 + UE.y_1_angle^2);
        r_est_2_All = sqrt((UE.x_1_angle-BS.xAll).^2 +...
            (UE.y_1_angle-BS.yAll).^2);
        %% Refine the angle through Gradient Descent
        [theta_UE2BS_final] = Gradient_descent_rel_Armijo(...
            K, Y_relative_all, A, B, N, M, fc, r_est, sin(pi/2-angle_BS+theta_UE2BS_cor_est), r_est_2_All);        
        fprintf('UE2BS angle estimated by correlation (in unit 1) = %10.8f\n', sin(pi/2-angle_BS+theta_UE2BS_cor_est));
        fprintf('UE2BS angle estimated by Gradient Descent (in unit 1) = %10.8f\n', theta_UE2BS_final);
        % if the estimated results estimated by the Gradient Descent is bad, then use the result estimated by the correlation
        H_1 = NearFieldH(K, N, sin(pi/2-angle_BS+theta_UE2BS_cor_est), r_est, r_est, 1./(r_est_2_All*4*pi).*(10.^(K.*r_est_2_All/10)), fc, B, M, 'absolute');
        H_2 = NearFieldH(K, N, theta_UE2BS_final, r_est, r_est, 1./(r_est_2_All*4*pi).*(10.^(K.*r_est_2_All/10)), fc, B, M, 'absolute');
        if norm(Y_absolute_all-A*H_1, 'fro') < norm(Y_absolute_all-A*H_2, 'fro')
            theta_UE2BS_final = sin(pi/2-angle_BS+theta_UE2BS_cor_est);
        end        
        fprintf('True UE2BS angle (in unit 1) = %10.8f\n', sin(pi/2-angle_BS+atan((UE.y-BS.y)/(UE.x-BS.x))));        
        fprintf('Final estimation of the UE2BS angle (in unit 1) = %10.8f\n', theta_UE2BS_final);
        %% localization through the AoA and hyperbola
        k_B2U = tan(asin(theta_UE2BS_final)+angle_BS-pi/2);
        m_1 = 1-(a/b.*k_B2U).^2;
        m_2 = -2*c;
        m_3 = b*b;
        UE.x_1_angle = (-m_2+sqrt(m_2.^2-4.*m_1*m_3))./(2.*m_1);
        UE.y_1_angle = UE.x_1_angle.*k_B2U;
        %% calculate RMSE
        RMSE.theta(i_B, i_it) = norm(asin(theta_UE2BS_final)+angle_BS-pi/2 - theta_UE2BS);
        RMSE.x(i_B, i_it) = norm(UE.x_1_angle - UE.x);
        RMSE.y(i_B, i_it) = norm(UE.y_1_angle - UE.y);
        RMSE.r(i_B, i_it) = norm(sqrt(UE.x^2 + UE.y^2) - ...
            sqrt(UE.x_1_angle^2 + UE.y_1_angle^2));
        fprintf('RMSE.x = %10.6f\n', RMSE.x(i_B, i_it));
        fprintf('RMSE.y = %10.6f\n', RMSE.y(i_B, i_it));
        fprintf('RMSE.r = %10.6f\n', RMSE.r(i_B, i_it));
        fprintf('RMSE.theta = %10.8f\n', RMSE.theta(i_B, i_it));
    end
end