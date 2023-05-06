clear
clc
addpath('.\function')
fprintf('==============================================\n');
rng(666, 'twister');
%% parameter setting
N               =           256;                % number of BS antennas
NRF             =           4;
fc              =           100e9;              %100GHz
lambda_c        =           3e8/fc;
d               =           lambda_c/2;
B               =           10e9;               %10GHz
M               =           64;                 %number of carrier frequency
K               =           -0.45/1000;% -0.45dB/km
S               =           sqrt(4*pi);
% S               =           18;
AntGain         =           19.1+1.34;%dB  
angleOverSmp    =           2;

P_BS                           = 16;% number observed time slots

% Nit             =           5e1;
Nit             =           1;
%% -    BS£¬RIS£¬UE, scatter location setting
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
% location of the UE
UE.x                         = 40.1;
UE.y                         = -20.1;

r_UE2BS                     = sqrt((UE.x-BS.x).^2 +...
                                    (UE.y-BS.y).^2);   
r_UE2BS_All                     = sqrt((UE.x-BS.xAll).^2 +...
                                    (UE.y-BS.yAll).^2);  

theta_UE2BS                 = atan((UE.y-BS.y)/(UE.x-BS.x));         

%% -    transmit power, noise power
P_noise = 10^(-17.4)*(B/2048)*(M);

P_t_dBm = (15:5:45);%dBmW

NMSE_GenieLS = zeros(length(P_t_dBm), 1);
%% -    -   scatterer
num_scatter_UE2BS = 3;
PathInCluster = 3;
SizeOfScatter = 1;
Scatter_UE2BS = genScatterUEBS(K, S, N, BS, UE, num_scatter_UE2BS,...                
                PathInCluster,...
                SizeOfScatter);
%% parameter output
fprintf('P_BS = %d\n', P_BS);
fprintf('num_scatter_UE2BS = %d\n', num_scatter_UE2BS);
fprintf('true UE2BS distance = %10.6f\n', r_UE2BS);                            
fprintf('true UE2BS angle (rad) = %10.8f\n', theta_UE2BS); 
% (with respect to the normal line of the BS)
fprintf('true UE2BS angle (has no unit) = %10.8f\n', sin(pi/2-angle_BS+theta_UE2BS)); 
%% localization
for i_P_t_dBm = 1:length(P_t_dBm)
    for i_it = 1:Nit
        fprintf('==============================================\n');        
        fprintf('P_t_dBm = %3d\n',P_t_dBm(i_P_t_dBm));
        fprintf('i_it/Nit = %3d/%3d\n', i_it, Nit);
%% -    generate Y, A, H according to the parameters defined before
%% -    -   UE2BS
        H_UE2BS_absolute = NearFieldH(K, N, sin(pi/2-angle_BS+[theta_UE2BS, Scatter_UE2BS.theta])...
            , [r_UE2BS, Scatter_UE2BS.r], [r_UE2BS, Scatter_UE2BS.r_LastHop], [1./(r_UE2BS_All*4*pi).*(10.^(K.*r_UE2BS_All/10)) Scatter_UE2BS.alpha], fc, B, M, 'absolute');
%         Y_absolute_all = zeros(NRF*P_BS, M);
        A = exp(1j*2*pi*rand(NRF*P_BS, N))/sqrt(NRF);
        A(1,:) = zeros(1, N);
        if mod(N,2) == 0%N is even
            A(1, ceil(N/2)) = 1/sqrt(NRF);
            A(1, ceil((N+1)/2)) = 1/sqrt(NRF);
        else
            A(1, ceil((N+1)/2)) = 1/sqrt(NRF);
        end
        
        Noise = (1/sqrt(2)*randn(size(A, 1), M) + 1j/sqrt(2)*randn(size(A, 1), M));
        Y_absolute_all = (A*H_UE2BS_absolute + sqrt(P_noise)/sqrt(10^(P_t_dBm(i_P_t_dBm)/10))/sqrt(10^(AntGain/10))*Noise);
%         10*log10(10^(AntGain/10)*10.^(P_t_dBm/10).*norm(A*H_UE2BS_absolute,'fro')^2./norm(sqrt(P_noise)*Noise,'fro')^2)
%         Y_absolute_all = A*H_UE2BS_absolute;
%% channel estimation and calulate the NMSE                
        H_hat_GenieLS = CE_GenieLS(Y_absolute_all, A,...
                    sin(pi/2-angle_BS+[theta_UE2BS, Scatter_UE2BS.theta]), ...
                    [r_UE2BS, Scatter_UE2BS.r],...
                    [r_UE2BS, Scatter_UE2BS.r_LastHop],...        
                    [1./(r_UE2BS_All*4*pi) Scatter_UE2BS.alpha],...
                    fc, B, K);                           
        NMSE_GenieLS(i_P_t_dBm) =  NMSE_GenieLS(i_P_t_dBm) + ...
                    (norm(H_hat_GenieLS - H_UE2BS_absolute, 'fro')/norm(H_hat_GenieLS, 'fro')).^2;                               
        fprintf("GenieLS average NMSE  =%ddB\n", 10*log10(NMSE_GenieLS(i_P_t_dBm)/i_it));
    end
    NMSE_GenieLS(i_P_t_dBm, :) = NMSE_GenieLS(i_P_t_dBm, :)/Nit;    
end
% figure; plot(P_t_dBm, 10*log10(NMSE_GenieLS(:,1)));