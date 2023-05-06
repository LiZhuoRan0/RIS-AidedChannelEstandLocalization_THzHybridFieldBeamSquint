function [hat_x] = SOMP_MMV_CE(A, y, k)
%% Orthogonal Matching Pursuit
% Inputs:
%       A                       =       sensing matrix£¨P*NRF¡ÁQ¡ÁM£©
%       y                       =       received signal£¨P*NRF¡ÁM£©
%       k                       =       maximum iteration number

% output:
%       hat_x                   =       estimated Channel
%%
L = 1 ;
[P_NRF, Q, M] = size(A);
[P_NRF, M] = size(y);   
r0=y;                                         
hat_x = zeros(Q, M);
product = zeros(Q, M);
pos = [];
% for times = 1:k*L_1
for times = 1:k
    for i_M = 1:M
        product(:, i_M) = abs(A(:, :, i_M)'*r0(:, i_M));
    end 
    [~, pos_new] = sort(sum(abs(product).^2,2), 'descend');
    pos = union(pos, pos_new(1:L));
    Aug_t = A(:, pos, :);                    
    for i_M = 1:M
        LS=Aug_t(:, :, i_M)\y(:, i_M);       
        r0(:, i_M)=y(:, i_M)-Aug_t(:, :, i_M)*LS;                      
    end
end    
for i_M = 1:M
    hat_x(pos, i_M) = Aug_t(:, :, i_M)\y(:, i_M);          
end

end

