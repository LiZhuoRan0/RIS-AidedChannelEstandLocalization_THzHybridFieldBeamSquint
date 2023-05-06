function [theta_final] = Gradient_descent_rel_Armijo(...
    K, Y, A, B, N, M, fc, r_est, theta_est, r_UE2BS_All)
%% parameter setting
alpha_ini       = 1e8;% initial step size
precision       = 1e-7;% estimation accuracy
% delta_precision = 1e-4;
Nit             = 1e2;
% c               = 1e5/7.5;
c               = 1e-4;

lambda_c = 3e8/fc;
d = lambda_c/2;
delta = ((2*(0:N-1)-N+1)/2).';
d_h_d_theta = zeros(N, M);
f = fc + (-B/2 : B/M :(B/2 - B/M));
lambda = 3e8./f;

H = NearFieldH(K, N, theta_est, r_est, r_est, 1./(r_UE2BS_All*4*pi).*(10.^(K.*r_UE2BS_All/10)), fc, B, M, 'relative');
theta_final = theta_est;
theta_pre = 0;% previous estimated angle

Y_vec = reshape(Y, [], 1);

%% 
for num_it =  1:Nit
    %% -    derivative
    alpha = alpha_ini;
    for i_M = 1:M
        d_h_d_theta(:, i_M) = 2j*pi/(lambda(i_M))...
                *exp(-1j*2*pi/lambda(i_M)*(sqrt(r_est.^2+(delta*d).^2-2*r_est*theta_final*delta*d)-r_est))...
                *r_est.*delta*d...
                ./sqrt(r_est.^2+(delta*d).^2-2*r_est*theta_final*delta*d)...
                .*lambda(i_M)./(r_UE2BS_All*4*pi).*(10.^(K.*r_UE2BS_All/10))/sqrt(N);
    end    
    A_d_h_d_theta_vec = reshape(A*d_h_d_theta, [], 1);   
    A_H_vec = reshape(A*H, [], 1);

    delta_lzr = -2*real(Y_vec'*A_d_h_d_theta_vec)...
            +2*real(A_d_h_d_theta_vec'*A_H_vec);
%% Armijo
    H = NearFieldH(K, N, theta_final, r_est, r_est, 1./(r_UE2BS_All*4*pi).*(10.^(K.*r_UE2BS_All/10)), fc, B, M, 'relative');
    norm_pre = norm(Y - A*H, 'fro');
    % constraint on alpha*rho*delta_lzr
    while abs(alpha*delta_lzr) > 1e-3
        alpha = alpha/2;
    end
    while abs(alpha * delta_lzr) > precision
        H = NearFieldH(K, N, theta_final-alpha*delta_lzr, r_est, r_est, 1./(r_UE2BS_All*4*pi).*(10.^(K.*r_UE2BS_All/10)), fc, B, M, 'relative');
        norm_tmp = norm(Y - A*H, 'fro');
        if norm_tmp^2 < norm_pre^2 - c*abs(alpha * delta_lzr^2)
            break;
        else
            alpha = alpha * 0.5;
        end
    end
    theta_final = theta_final-alpha*delta_lzr;
    H = NearFieldH(K, N, theta_final, r_est, r_est, 1./(r_UE2BS_All*4*pi).*(10.^(K.*r_UE2BS_All/10)), fc, B, M, 'relative');
        
%     fprintf('it     = %4d\n', num_it);
%     fprintf('theta  = %10.9f\n', theta_final);
%     fprintf('Loss   = %20.19f\n', norm(Y - A*H, 'fro'));
%     fprintf('**********************************\n');
    
    if abs(theta_pre - theta_final) < precision
%     if abs(delta_lzr) < delta_precision
        break
    else
        theta_pre = theta_final;
    end
end
end