function [sup_max, flag] = gOMP_MMV(A, y, k)
%% generalized Orthogonal Matching Pursuit

% Inputs:
%       A         = sensing matrix£¨P*NRF¡ÁQ¡ÁM£©
%       y         = received signal£¨P*NRF¡ÁM£©
%       k         = maximum iteration number

% output:
%      sup_max          = index

%% 
flag = 1;
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
    if times == 1
        Cor_tmp = sum(abs(product).^2,2);
        [val_max, sup_max] = max(Cor_tmp);
        if val_max < 2*mean(Cor_tmp)
            flag = 0;
        end
        break;
    end
    [~, pos_new] = sort(sum(abs(product).^2,2), 'descend');

    pos = union(pos, pos_new(1:L));

    Aug_t = A_copy(:, pos, :);                      

    for i_M = 1:M
        LS=pinv(Aug_t(:, :, i_M))*y(:, i_M);          
        r0(:, i_M)=y(:, i_M)-Aug_t(:, :, i_M)*LS;                             
    end
end    
end

