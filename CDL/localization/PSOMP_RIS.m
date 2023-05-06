function [angle_est, flag] = PSOMP_RIS(Y, A_cmb, fc, B, N, NRIS, M, ...
                            L_hat, NRF, P, angleOverSmp)
%% angle estimation by the correlation of the OMP
% Inputs:
%       L_hat                           =           ~
%       Y                               =           Received pilot
%       A_cmb                           =           sensing matrix
%       N                               =           number of BS array elements
%       NRIS                            =           number of RIS array elements
%       fc                              =           center carrier frequency
%       NRF                             =           number of RF chains
%       P                               =           number of time slots
%       angleOverSmp                	=           angle oversampling coefficient
%       M                               =           number of carriers
%       B                               =           bandwidth

% Outputs:
%       angle_est                       =           estimated angle
%       flag                            =           ~
%% 
lambda_c = 3e8/fc;
d = lambda_c/2;
AngleRange = 1;
% W(NRIS¡ÁQ¡ÁM)

[W, s_final, Z_Delta] = PolarCodebook_wideband_angleBlock(AngleRange, NRIS,...
                            d, lambda_c, angleOverSmp, M, B);

[NRIS, Q, M] = size(W);

% Phi_bar = zeros(P, Q, M);
Phi_bar = zeros(P*NRF, Q, M);
for i_M = 1:M
    Phi_bar(:, :, i_M) = A_cmb(:, :, i_M)*W(:, :, i_M);
end    
%% 
[sup_max, flag] = gOMP_MMV(Phi_bar, Y, L_hat);


angle_est =(2*ceil(sup_max/s_final)-NRIS*angleOverSmp*AngleRange-1)/(NRIS*angleOverSmp);
end