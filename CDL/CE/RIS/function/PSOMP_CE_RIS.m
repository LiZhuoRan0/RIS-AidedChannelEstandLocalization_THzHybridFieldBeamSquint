function  [H_hat, it_num] = PSOMP_CE_RIS(L_2, threshold, Y, A, L_max,...
                        N, fc, NRF, P, angleOverSmp, M, B)                  
%% channel estimation without the UE's location
% Inputs:
%       L_max                           =           maximum iteration number
%       L_2                             =           number of atoms added in the support set in estimating the NLoS paths
%       threshold                       =           adaptive itration stop condition
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
%       H_hat                           =           estimated near-field Channel
%       times                           =           iteration number
%%
lambda_c = 3e8/fc;
d = lambda_c/2;
AngleRange = 1;

[W, s_final, Z_Delta] = PolarCodebook_narrowband_angleBlock(AngleRange, N,...
                            d, lambda_c, angleOverSmp, M, B);

[N, Q, M] = size(W);

Y_bar = Y;
Phi_bar = zeros(P*NRF, Q, M);
for i_M = 1:M
    Phi_bar(:, :, i_M) = A(:, :, i_M)*W(:, :, i_M);
end    

[hat_H, it_num] = gOMP_MMV_CE(L_2, threshold, Phi_bar, Y_bar, L_max);

H_hat = zeros(N, M);
for i_M = 1:M
    H_hat(:, i_M) = W(:, :, i_M)*hat_H(:, i_M);
end    
                    
end