function [H_hat] = CE_GenieLS(Y, A, sin_theta, r, r_LastHop,  alpha,  fc, B, K)
Nsrs = length(r);
M = size(Y, 2);
N = size(A, 2);

H = zeros(N, M, Nsrs);
LS = zeros(Nsrs, M);
H_hat = zeros(N, M);

for i_Nsrs = 1:Nsrs
    H(:, :, i_Nsrs) = NearFieldH(K, N, sin_theta(i_Nsrs), r(i_Nsrs), r_LastHop(i_Nsrs), alpha(:,i_Nsrs), fc, B, M, 'absolute');
end
for m = 1:M
    LS(:, m) = (A*squeeze(H(:,m,:)))\Y(:,m);
    H_hat(:,m) = squeeze(H(:,m,:))*LS(:, m);
end

end

