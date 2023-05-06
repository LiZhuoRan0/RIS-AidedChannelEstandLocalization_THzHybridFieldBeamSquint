function  [H_hat, it_num] = PgOMP_CE_RIS_AS(L_2, threshold, Y, A, L_max,...
                        N, fc, NRF, P, angleOverSmp, M, B, theta_est, r_est)             
%% channel estimation with the UE's location              
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
%       it_num                           =           iteration number
%%
lambda_c = 3e8/fc;
d = lambda_c/2;
AngleRange = 1;

[W, ~, ~] = PolarCodebook_wideband_angleBlock(AngleRange, N,...
                            d, lambda_c, angleOverSmp, M, B);

[N, Q, M] = size(W);

Y_bar = Y;
% Phi_bar = zeros(P, Q, M);
Phi_bar = zeros(P*NRF, Q, M);
for i_M = 1:M
    Phi_bar(:, :, i_M) = A(:, :, i_M)*W(:, :, i_M);
end    
%% update the dictionary according to the UE's locaiton
Phi_bar_new = zeros(P*NRF, Q+1, M);
Phi_bar_new(:,1:Q,:) = Phi_bar;
for i_M = 1:M
    Phi_bar_new(:, Q+1, i_M) = A(:,:,i_M)*genb(theta_est, r_est, N, fc-B/2+(i_M-1)*B/M, d);
end

W_new = zeros(N, Q+1, M);
W_new(:,1:Q,:) = W;
for i_M = 1:M
    W_new(:, Q+1, i_M) = genb(theta_est, r_est, N, fc-B/2+(i_M-1)*B/M, d);
end
%%
[hat_H, it_num] = gOMP_MMV_CE_2(L_2, threshold, Phi_bar_new, Y_bar, L_max);

H_hat = zeros(N, M);
for i_M = 1:M
    H_hat(:, i_M) = W_new(:, :, i_M)*hat_H(:, i_M);
end    
                 
end