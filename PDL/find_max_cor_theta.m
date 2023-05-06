function theta_UE2BS_cor_est = find_max_cor_theta(K, r_est_All, r_B2U_h, theta_B2U_h, Y_tmp_absolute, A, B, N, M, fc, angle_BS)
%% coarse angle estimation by correlation on the hyperbola

dict = zeros(size(Y_tmp_absolute, 1), length(theta_B2U_h), size(Y_tmp_absolute, 2));
cor = zeros(1, length(theta_B2U_h));
Y_tmp_absolute = Y_tmp_absolute/norm(Y_tmp_absolute, 'fro');
for i = 1:length(theta_B2U_h)
    theta = theta_B2U_h(i);
    r_est = r_B2U_h(i);
    dict(:, i, :) = A*NearFieldH(K, N, sin(pi/2-angle_BS+theta), r_est, r_est, 1./(r_est_All*4*pi).*(10.^(K.*r_est_All/10)), fc, B, M, 'absolute');
    dict(:, i, :) = dict(:, i, :)/norm(squeeze(dict(:, i, :)), 'fro');
    cor(i) = reshape(Y_tmp_absolute, [], 1)'*reshape(squeeze(dict(:, i, :)), [], 1);
end
[~, theta_UE2BS_cor_index_est] = max(cor);
theta_UE2BS_cor_est = theta_B2U_h(theta_UE2BS_cor_index_est);

fprintf('estimated UE2BS angle (rad)  = %10.6f\n', theta_UE2BS_cor_est); 
end
