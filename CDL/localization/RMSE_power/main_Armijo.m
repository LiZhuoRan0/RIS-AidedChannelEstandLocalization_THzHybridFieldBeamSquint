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
B               =           10e9;               %10GHz
M               =           64;                 %number of carrier frequency
K               =           -0.45/1000;% -0.45dB/km
S               =           sqrt(4*pi);
AntGain         =           19.1+1.34;%dB  
angleOverSmp    =           2;
% S_eff           =           2;
S_eff           =           (NRIS*d)^2;

P_BS                           = 16;% number observed time slots
P_RIS                           = 32;

Nit             =           4e1;
% Nit             =           10;
%% -    BS£¬RIS£¬UE, scatter location setting
% RIS location
RIS.x                          = 20*sqrt(2);
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

P_t_dBm_All = (15:10:45);%dBmW
%%
RMSE.x_co                      = zeros(1, length(P_t_dBm_All));
RMSE.y_co                      = zeros(1, length(P_t_dBm_All));
RMSE.r_co                      = zeros(1, length(P_t_dBm_All));
RMSE.theta_co                  = zeros(1, length(P_t_dBm_All));

RMSE.x                      = zeros(1, length(P_t_dBm_All));
RMSE.y                      = zeros(1, length(P_t_dBm_All));
RMSE.r                      = zeros(1, length(P_t_dBm_All));
RMSE.theta                  = zeros(1, length(P_t_dBm_All));
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
fprintf('P_BS = %d\n', P_BS);
fprintf('P_RIS = %d\n', P_RIS);
fprintf('num_scatter_UE2BS = %d\n', num_scatter_UE2BS);
fprintf('num_scatter_UE2RIS = %d\n', num_scatter_UE2RIS);
fprintf('true UE2BS distance = %10.6f\n', r_UE2BS);                            
fprintf('true UE2RIS distance = %10.6f\n', r_UE2RIS);
fprintf('true UE2BS angle (rad) = %10.8f\n', theta_UE2BS); 
fprintf('true UE2BS angle (no unit) = %10.8f\n', sin(pi/2-angle_BS+theta_UE2BS)); 
fprintf('true UE2RIS angle (rad) = %10.8f\n', theta_UE2RIS); 
fprintf('true UE2RIS angle (no unit) = %10.8f\n', sin(pi/2+angle_RIS-theta_UE2RIS));   
%% localization
for i_P_t_dBm = 1:length(P_t_dBm_All)  

    P_t_dBm = P_t_dBm_All(i_P_t_dBm);

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
        fprintf('i_it/Nit = %3d/%3d\n', i_it, Nit);
%% -    generate Y, A, H according to the parameters defined before
%% -    -   UE2BS 
        H_UE2BS_absolute = NearFieldH(K, N, sin(pi/2-angle_BS+[theta_UE2BS, Scatter_UE2BS.theta])...
            , [r_UE2BS, Scatter_UE2BS.r], [r_UE2BS, Scatter_UE2BS.r_LastHop], [1./(r_UE2BS_All*4*pi).*(10.^(K.*r_UE2BS_All/10)) Scatter_UE2BS.alpha], fc, B, M, 'absolute');
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
%% -    generate the polar-domain dictionary, obtain the initial angle estimation (after sin)
        theta_UE2BS_coarse = PSOMP(Y_absolute_all, A, 1,...
                             N, fc, NRF, P_BS, angleOverSmp, M, B);
        theta_UE2RIS_coarse = PSOMP_RIS(Y_RIS_Noise, A_RIS, fc, B, N, NRIS, M, ...
                            1, NRF, P_RIS, angleOverSmp);
        fprintf('UE2BS angle estimation using correlation (no unit) = %10.8f\n', theta_UE2BS_coarse);
        fprintf('UE2RIS angle estimation using correlation (no unit) = %10.8f\n', theta_UE2RIS_coarse);
        UE_est_co = jointLoc(asin(theta_UE2BS_coarse)+angle_BS-pi/2, -asin(theta_UE2RIS_coarse)+pi/2+angle_RIS, BS, RIS);
        r_est_co = sqrt((UE_est_co.x - BS.x).^2+(UE_est_co.y - BS.y).^2);
        r_est_RIS_co = sqrt((UE_est_co.x - RIS.x).^2+(UE_est_co.y - RIS.y).^2);
        r_est_All = sqrt((UE_est_co.x-BS.xAll).^2 +...
                                    (UE_est_co.y-BS.yAll).^2);
        r_est_RIS_All = sqrt((UE_est_co.x-RIS.xAll).^2 +...
                                    (UE_est_co.y-RIS.yAll).^2); 
%% -    fine angle estimation
%% -    -   UE2BS, Gradient Descent 
        [theta_UE2BS_fine] = Gradient_descent_rel_Armijo(...
            K, Y_relative_all, A, B, N, M, fc, r_est_co, theta_UE2BS_coarse, r_est_All);       
        fprintf('UE2BS angle estimation using Gradient Descent (no unit) = %10.8f\n', theta_UE2BS_fine);
%% -    -   UE2RIS, Hierarchical search
        theta_UE2RIS_fine = HieDic(K, Y_RIS_Noise, A_RIS, fc, B, r_est_RIS_co, r_est_RIS_All, theta_UE2RIS_coarse);                               
        fprintf('UE2RIS angle estimation using hierarchical dictionary (no unit) = %10.8f\n', theta_UE2RIS_fine);
%% -    calculate the coordinate  
        UE_est = jointLoc(asin(theta_UE2BS_fine)+angle_BS-pi/2, -asin(theta_UE2RIS_fine)+pi/2+angle_RIS, BS, RIS);
        % relative to the BS
        r_est = sqrt((UE_est.x - BS.x).^2+(UE_est.y-BS.y).^2);
        theta_est = -sin(atan((UE_est.y - BS.y)/(UE_est.x - BS.x)));
%% calculate RMSE
        RMSE.theta_co(i_P_t_dBm) = RMSE.theta_co(i_P_t_dBm) + norm(asin(theta_UE2BS_coarse)+angle_BS-pi/2 - theta_UE2BS).^2;
        RMSE.x_co(i_P_t_dBm) = RMSE.x_co(i_P_t_dBm) + norm(UE_est_co.x - UE.x)^2;
        RMSE.y_co(i_P_t_dBm) = RMSE.y_co(i_P_t_dBm) + norm(UE_est_co.y - UE.y)^2;
        RMSE.r_co(i_P_t_dBm) = RMSE.r_co(i_P_t_dBm) + norm(sqrt(UE.x^2 + UE.y^2) - ...
                        sqrt(UE_est_co.x^2 + UE_est_co.y^2))^2;
        fprintf('RMSE.x_co = %10.8f\n', sqrt(RMSE.x_co(1, i_P_t_dBm)/i_it));
        fprintf('RMSE.y_co = %10.8f\n', sqrt(RMSE.y_co(1, i_P_t_dBm)/i_it));
        fprintf('RMSE.r_co = %10.8f\n', sqrt(RMSE.r_co(1, i_P_t_dBm)/i_it));        
        fprintf('RMSE.theta_co = %10.8f\n', sqrt(RMSE.theta_co(1, i_P_t_dBm)/i_it));
        
        RMSE.theta(i_P_t_dBm) = RMSE.theta(i_P_t_dBm) + norm(asin(theta_UE2BS_fine)+angle_BS-pi/2 - theta_UE2BS).^2;
        RMSE.x(i_P_t_dBm) = RMSE.x(i_P_t_dBm) + norm(UE_est.x - UE.x)^2;
        RMSE.y(i_P_t_dBm) = RMSE.y(i_P_t_dBm) + norm(UE_est.y - UE.y)^2;
        RMSE.r(i_P_t_dBm) = RMSE.r(i_P_t_dBm) + norm(sqrt(UE.x^2 + UE.y^2) - ...
                        sqrt(UE_est.x^2 + UE_est.y^2))^2;
        fprintf('RMSE.x = %10.8f\n', sqrt(RMSE.x(1, i_P_t_dBm)/i_it));
        fprintf('RMSE.y = %10.8f\n', sqrt(RMSE.y(1, i_P_t_dBm)/i_it));
        fprintf('RMSE.r = %10.8f\n', sqrt(RMSE.r(1, i_P_t_dBm)/i_it));        
        fprintf('RMSE.theta = %10.8f\n', sqrt(RMSE.theta(1, i_P_t_dBm)/i_it));        
        fprintf('current RMSE.theta = %10.8f\n', norm(asin(theta_UE2BS_fine)+angle_BS-pi/2 - theta_UE2BS));
    
   end
end
RMSE.x_co = sqrt(RMSE.x_co/Nit);
RMSE.y_co = sqrt(RMSE.y_co/Nit);
RMSE.r_co = sqrt(RMSE.r_co/Nit);
RMSE.theta_co = sqrt(RMSE.theta_co/Nit);

RMSE.x = sqrt(RMSE.x/Nit);
RMSE.y = sqrt(RMSE.y/Nit);
RMSE.r = sqrt(RMSE.r/Nit);
RMSE.theta = sqrt(RMSE.theta/Nit);