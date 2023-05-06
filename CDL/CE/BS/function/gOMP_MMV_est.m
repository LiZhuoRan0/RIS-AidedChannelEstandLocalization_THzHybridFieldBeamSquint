function [hat_x] = gOMP_MMV_est(A, y, k)
%% generalized Orthogonal Matching Pursuit

% Inputs:
%       A         = sensing matrix£¨P*NRF¡ÁQ¡ÁM£©
%       y         = received signal£¨P*NRF¡ÁM£©
%       k         = maximum iteration number

% output:
%      hat_x          = sparse approximation£¨N¡ÁM£©
%% 
L = 6;
[P_NRF, Q, M] = size(A);
[P_NRF, M] = size(y);   
A_copy = A;
r0=y;                       
hat_x = zeros(Q, M);
product = zeros(Q, M);
pos = [];
for times = 1:k
    for i_M = 1:M
        product(:, i_M) = abs(A(:, :, i_M)'*r0(:, i_M));
    end 
    [~, pos_new] = sort(sum(abs(product).^2,2), 'descend');

    pos = union(pos, pos_new(1:L));

    Aug_t = A_copy(:, pos, :);                     

    for i_M = 1:M
        LS=pinv(Aug_t(:, :, i_M))*y(:, i_M);          
        r0(:, i_M)=y(:, i_M)-Aug_t(:, :, i_M)*LS;                         
    end
end    
for i_M = 1:M
    hat_x(pos, i_M) = pinv(Aug_t(:, :, i_M))*y(:, i_M);           
end
end

