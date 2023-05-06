function [angle_est, flag] = PSOMP(Y, A, L_max,...
                        N, fc, NRF, P, angleOverSmp, M, B)
%% angle estimation by the correlation of the OMP
% Inputs:
%       L_max                           =           maximum iteration number
%       Y                               =           Received pilot
%       A                               =           sensing matrix
%       N                               =           number of array elements
%       fc                              =           center carrier frequency
%       NRF                             =           number of RF chains
%       P                               =           number of time slots
%       angleOverSmp                	=           angle oversampling coefficient
%       M                               =           number of carriers
%       B                               =           bandwidth

% Outputs:
%       angle_est                       =           estimated angle
%%
lambda_c = 3e8/fc;
d = lambda_c/2;
AngleRange = 1;

[W, s_final, Z_Delta] = PolarCodebook_wideband_angleBlock(AngleRange, N,...
                            d, lambda_c, angleOverSmp, M, B);

[N, Q, M] = size(W);

Y_bar = Y;
Phi_bar = zeros(P*NRF, Q, M);
for i_M = 1:M
    Phi_bar(:, :, i_M) = A*W(:, :, i_M);
end    


% estimate the index
[sup_max] = gOMP_MMV(Phi_bar, Y_bar, L_max);

% Get the Angle by the index of the atom in the support set
angle_est =(2*ceil(sup_max/s_final)-N*angleOverSmp*AngleRange-1)/(N*angleOverSmp);                  
end