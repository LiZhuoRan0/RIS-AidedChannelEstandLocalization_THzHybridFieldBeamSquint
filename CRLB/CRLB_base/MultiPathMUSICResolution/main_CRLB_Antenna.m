%% 
clear
clc
% addpath('.\func')
rng(666);
%% 参数设置
Nall               =           [4 8 16 32 64 128 256 512 1024];                % 基站处天线
fc              =           100e9;              %100GHz
lambda_c        =           3e8/fc;
d               =           lambda_c/2;
T               =           200;
Nit             =           1e2;
% SNRdBs          =           0:5:20;
% SNRdBs          =           100:5:105;
SNRdB           =   20;
NumSrs          =           2;
lzr = zeros(NumSrs, Nit);
%% 
% theta   = 2*rand(1, NumSrs) - 1;
theta = [0.00100 0];
theta = sort(theta);
alpha = ones(NumSrs,1);
% theta   = 0.5;%pi/6

CRLB        = zeros(length(Nall), NumSrs);


for i_N = 1:length(Nall)
    N = Nall(i_N);
    H   = genSteerVector(theta, N, d, lambda_c);
    Y   = zeros(N, T);
%     MseMUSIC    = zeros(length(SNRdBs), NumSrs);
%     MseESPRIT    = zeros(length(SNRdBs), NumSrs);
%     MseCorr    = zeros(length(SNRdBs), 1);
    
    a1  = genPartialSteerVector(theta(1), N, d, lambda_c, 1);
    a2  = genPartialSteerVector(theta(2), N, d, lambda_c, 1);
    
    D = [a1 a2];
    Cst = D'*(eye(N)-H*inv(H'*H)*H')*D;

    fprintf('==============================================\n');
    fprintf('SNR        = %ddB\n', SNRdB);
    for i_it = 1:Nit
%         if i_it == 71 && i_SNR == 2
%             pause
%         end
        if mod(i_it,10) == 0
            fprintf('i_it/Nit   = %3d/%3d\n', i_it, Nit);
        end
        X = (sqrt(2)/2*randn(NumSrs,T)+1j*sqrt(2)/2*randn(NumSrs,T));
        for i_T = 1:T
%             Y(:, i_T) = awgn(H*(alpha.*X(:,i_T)), SNRdBs(i_SNR), 'measured');
            Y(:, i_T) = awgn(H*(alpha.*X(:,i_T)), SNRdB);
        end
%         Y = H;
        
%         theta_MUISC     = MUSIC_MultiPath(Y, NumSrs);  
%         psi = TLS_ESPRIT_Algorithm(Y, NumSrs);
%         theta_ESPRIT = log(psi)/(1j*pi);
%         theta_Corr     = Corr(Y);
%         for i_Srs  = 1:length(theta_MUISC)
%             MseMUSIC(i_SNR, i_Srs)     = MseMUSIC(i_SNR, i_Srs) + abs(asin(theta_MUISC(i_Srs)) - asin(theta(i_Srs)))^2;
%             MseESPRIT(i_SNR, i_Srs)     = MseESPRIT(i_SNR, i_Srs) + abs(asin(theta_ESPRIT(i_Srs)) - asin(theta(i_Srs)))^2;
%         end
%         MseCorr(i_SNR)     = MseCorr(i_SNR) + abs(asin(theta_Corr) - asin(theta))^2;
        sigma2  = 10^(-SNRdB/10);
%         sigma2  = (norm(Y, 'fro')^2-norm(H*X, 'fro')^2)/...
%                 (size(Y, 1)*size(Y, 2));        

        tmp1  = zeros(2,2);
        for i_t = 1:T%该for循环可以优化
            XX = diag([X(1, i_t) X(2, i_t)]);
            tmp2 = XX'*Cst*XX;
            tmp1 = tmp1 + tmp2;
        end
        tmp3 = inv(real(tmp1));
        tmp4 = diag(tmp3);
        tmp5 = tmp4(:).';
        CRLB(i_N, :)     = CRLB(i_N, :) + sigma2/2.*tmp5;
%         if i_SNR == 1
%             lzr(:, i_it) = tmp5.'*1e6;
%         end
    end
end
% MseMUSIC     = MseMUSIC/Nit;
% MseESPRIT     = MseESPRIT/Nit;
% MseCorr     = MseCorr/Nit;
CRLB         = CRLB/Nit;
%% 绘图
set(0,'defaultfigurecolor','w') 
figure; hold on; grid on; box on;
xlabel('Number of Antennas');
ylabel('MSE(rad^2)');
title(['$\sigma^2$ =' num2str(sigma2)], 'Interpreter', 'latex')
set(gca, 'YScale', 'log');

% plot(SNRdBs, MseMUSIC,'g:o', 'LineWidth', 2);
% plot(SNRdBs, MseMUSIC, 'LineWidth', 2);
% plot(SNRdBs, MseESPRIT, 'LineWidth', 2);
% plot(SNRdBs, MseCorr,'b:>', 'LineWidth', 2);
plot(Nall, CRLB(:,1), 'b:>', 'LineWidth', 2);
plot(Nall, CRLB(:,2), 'r:o', 'LineWidth', 2);
legend
l = legend({   
%         'MUSIC~~1'...
%         'MUSIC~~2' ...
%         'ESPRIT~~1'...
%         'ESPRIT~~2' ...
        'CRLB~~~1' ...
        'CRLB~~~2'},...  
        'Interpreter','latex', 'Box','off'); 
l.FontSize = 12; 
%         'Corr'},... 
