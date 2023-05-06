function [hat_x, it_num] = gOMP_MMV_CE(L_2, threshold, A, y, L_max)
%% generalized Orthogonal Matching Pursuit

% Inputs:
%       A         = sensing matrix£¨P*NRF¡ÁQ¡ÁM£©
%       y         = received signal£¨P*NRF¡ÁM£©
%       L_max     = maximum iteration number
%       L_2       = number of atoms added in the support set in estimating the NLoS paths
%       threshold = adaptive itration stop condition

% output:
%      hat_x          = sparse approximation£¨N¡ÁM£©
%      it_num         = iteration number
%%
L_1 = 1;
[P_NRF, Q, M] = size(A);
[P_NRF, M] = size(y);   
r0=y;         
hat_x = zeros(Q, M);
product = zeros(Q, M);
pos = [];
pos_pre = [];
r_norm = realmax;
ratio_pre = realmax;
for times = 1:L_max
% for times = 1:1
    for i_M = 1:M
        product(:, i_M) = abs(A(:, :, i_M)'*r0(:, i_M));
    end 
    [~, pos_new] = sort(sum(abs(product).^2,2), 'descend');
    if times == 1
        pos = union(pos, pos_new(1:L_1));
    else
        pos = union(pos, pos_new(1:L_2));
    end
    Aug_t = A(:, pos, :);                       
    for i_M = 1:M
        LS=Aug_t(:, :, i_M)\y(:, i_M);          
        r0(:, i_M)=y(:, i_M)-Aug_t(:, :, i_M)*LS;                   
    end
    
    if norm(r0, 'fro')/r_norm > threshold || (L_2 ~= 1 && norm(r0, 'fro')/r_norm < 0.8 && times>=8)
%     if times > 2 && norm(r0, 'fro')/r_norm < ratio_pre 
        fprintf('current residual/previous residual = %10.5f\n',norm(r0, 'fro')/r_norm);
        pos = pos_pre;
        it_num = times -1;
        break;
    end
    ratio_pre = norm(r0, 'fro')/r_norm;
    fprintf('current residual/previous residual = %10.5f\n',ratio_pre);
    r_norm = norm(r0, 'fro');
    pos_pre = pos;
    it_num = times;
end    
fprintf('iteration number           = %3d\n',it_num);
Aug_t = A(:, pos, :);
for i_M = 1:M
    hat_x(pos, i_M) = Aug_t(:, :, i_M)\y(:, i_M);       
end
end

