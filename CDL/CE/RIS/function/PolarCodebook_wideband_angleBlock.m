function [Wfinal, s_final, Z_Delta] = PolarCodebook_wideband_angleBlock(AngleRange, N,...
                            d, lambda_c, angleOverSmp, M, B)
%% 宽带码本，每一个角度块有不同距离

% Inputs:
%       rho_min                     =       minimum allowable distance
%       beta_Delta                  =       threshold                   
%       N                           =       阵元数              
%       d                           =       antenna spacing             
%       lambda_c                    =       中心wavelength                  
%       angleOverSmp                =       角度域过采样系数
%       M                           =       子载波数
%       B                           =       带宽
% Outputs:
%       W                           =       polar-domain transform matrix
%       s_final                     =       number of distance grids
%       Z_Delata                    =       maximum distance in the near-field dictionary
%%
fc = 3e8/lambda_c;
s_final = 10;
Q = s_final*ceil(N*angleOverSmp*AngleRange);
Wfinal = zeros(N, Q, M);

Z_Delta = 50;
for i_M = 1:M%carrier frequency
    theta = zeros(ceil(N*angleOverSmp*AngleRange), 1);
    r = zeros(N, s_final);
    W = [];
    for n = 1:ceil(N*angleOverSmp*AngleRange)%angle
        theta(n) = (2*n - angleOverSmp*N - 1)/N/angleOverSmp;
         W_angle = [];
        for s = 0:s_final-1%distance           
            if s == 0%  far-field case
                a = genb(theta(n), 100, N, fc-B/2+(i_M-1)*B/M, d);
                W_angle = [W_angle a];              
            else
                r(n, s+1) = 1/s*Z_Delta*(1 - theta(n)*theta(n));              
                b = genb(theta(n), r(n, s+1), N, fc-B/2+(i_M-1)*B/M, d);
                W_angle = [W_angle b];
            end        
        end
        W = [W W_angle];
    end
    Wfinal(:, :, i_M) = W;
end

end