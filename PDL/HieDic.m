function theta_out = HieDic(K, Y, A, fc, B, r_co, r_est_RIS_All, theta_in)
%% refine on the basis of the coarse estimation
% inputs
%       r_co        coarse estimation of the distance
%       theta_in    coarse estimation of the angle

% outputs
%       theta_out   the refined angle after hierarchical search (after action of the sin)
%%
lambda_c = 3e8/fc;
d = lambda_c/2;
range = 0.002*[1 0.1 0.01];%1/512=0.002
number = 41;
theta_out = theta_in;
NRIS = size(A, 2);
M = size(A, 3);
W = zeros(NRIS, number, M);
Phi_bar = zeros(size(A, 1), number, M);
product = zeros(number, M);
for i = 1:length(range)
    theta = linspace(theta_out-range(i), theta_out+range(i), number);
    for i_n = 1:number 
        W(:, i_n, :) = NearFieldH(K, NRIS, theta(i_n)...            
                , r_co, r_co, 1./(r_est_RIS_All*4*pi).*(10.^(K.*r_est_RIS_All/10)), fc, B, M, 'absolute');
    end    
    for i_M = 1:M
        Phi_bar(:, :, i_M) = A(:, :, i_M)*W(:, :, i_M);
    end     
    for i_M = 1:M
        product(:, i_M) = abs(Phi_bar(:, :, i_M)'*Y(:, i_M));
    end

    Cor_tmp = sum(product.^2,2);
    [~, sup_max] = max(Cor_tmp);

    theta_out = 2*range(i)/(number-1)*(sup_max-1) + theta_out - range(i);
end
end