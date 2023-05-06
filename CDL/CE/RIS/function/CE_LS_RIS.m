function [H_hat] = CE_LS_RIS(Y, A)
M = size(A, 3);
N = size(A, 2);
H_hat = zeros(N, M);
for m = 1:M
    H_hat(:, m) = A(:,:,m)\Y(:, m);
end
end

