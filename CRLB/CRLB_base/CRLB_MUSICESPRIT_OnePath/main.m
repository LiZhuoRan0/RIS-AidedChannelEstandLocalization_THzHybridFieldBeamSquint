%% 
clear
clc
addpath('.\func')
rng(666);
%% parameter setting
N               =           32;                % number of BS antenna
fc              =           100e9;              %100GHz
lambda_c        =           3e8/fc;
d               =           lambda_c/2;
T               =           100;
Nit             =           1e2;
SNRdBs          =           -10:5:20;
% SNRdBs          =           20;
%% 
theta   = 0.5;%pi/6
H   = genSteerVector(theta, N, d, lambda_c);
Y   = zeros(N, T);
A   = eye(N);
a1  = genPartialSteerVector(theta, N, d, lambda_c, 1);
a2  = genPartialSteerVector(theta, N, d, lambda_c, 2);
MseMUSIC    = zeros(length(SNRdBs), 1);
MseCorr    = zeros(length(SNRdBs), 1);
MseESPRIT   = zeros(length(SNRdBs), 1);
MseESPRIT_LS   = zeros(length(SNRdBs), 1);
CRLB        = zeros(length(SNRdBs), 1);
D           = a1;
Cst         = D'*(eye(N)-H*inv(H'*H)*H')*D;
for i_SNR = 1:length(SNRdBs)
    fprintf('==============================================\n');
    fprintf('SNR        = %ddB\n', SNRdBs(i_SNR));
    for i_it = 1:Nit
        
        if mod(i_it,10) == 0
            fprintf('i_it/Nit   = %3d/%3d\n', i_it, Nit);
        end
        X = (sqrt(2)/2*randn(1,T)+1j*sqrt(2)/2*randn(1,T));
        for i_T = 1:T
            Y(:, i_T) = awgn(H*X(i_T), SNRdBs(i_SNR), 'measured');
%             Y(:, i_T) = H*X(i_T);
        end
        
        theta_MUISC     = MUSIC(Y);
        theta_Corr     = Corr(Y);
        psi = TLS_ESPRIT_Algorithm(Y, 1);
        psi_LS = LS_ESPRIT_Algorithm(Y, 1);
        theta_ESPRIT = log(psi)/(1j*pi); 
        theta_ESPRIT_LS = log(psi_LS)/(1j*pi); 
        MseMUSIC(i_SNR)     = MseMUSIC(i_SNR) + abs(asin(theta_MUISC) - asin(theta))^2;
        MseCorr(i_SNR)     = MseCorr(i_SNR) + abs(asin(theta_Corr) - asin(theta))^2;
        MseESPRIT(i_SNR)    = MseESPRIT(i_SNR) + abs(asin(theta_ESPRIT) - asin(theta))^2;
        MseESPRIT_LS(i_SNR)    = MseESPRIT_LS(i_SNR) + abs(asin(theta_ESPRIT_LS) - asin(theta))^2;
%         MseMUSIC(i_SNR)     = MseMUSIC(i_SNR) + abs(pi*(theta_MUISC) - pi*(theta))^2;
%         MseCorr(i_SNR)     = MseCorr(i_SNR) + abs(pi*(theta_Corr) - pi*(theta))^2;
%         MseESPRIT(i_SNR)    = MseESPRIT(i_SNR) + abs(pi*(theta_ESPRIT) - pi*(theta))^2;
%         MseESPRIT_LS(i_SNR)    = MseESPRIT_LS(i_SNR) + abs(pi*(theta_ESPRIT_LS) - pi*(theta))^2;
        
%         sigma2  = 10^(-SNRdBs(i_SNR)/10);
% these two values are asymptotically consistent
        sigma2  = (norm(Y, 'fro')^2-norm(H*X, 'fro')^2)/...
                (size(Y, 1)*size(Y, 2));
%         X_bar   = kron(X.', eye(N));
%         y       = reshape(Y,[],1);
%         CRLB(i_SNR)     = CRLB(i_SNR) + sigma2/2./...
%                             real(-y'*X_bar*a2 + a2'*(X_bar'*X_bar)*H + a1'*(X_bar'*X_bar)*a1);
        CRLB(i_SNR)     = CRLB(i_SNR) + sigma2/2/real((Cst*(X*X')));
    end
end
MseMUSIC     = MseMUSIC/Nit;
MseCorr     = MseCorr/Nit;
MseESPRIT    = MseESPRIT/Nit;
CRLB         = CRLB/Nit;

%% plot
set(0,'defaultfigurecolor','w') 
figure; hold on; grid on; box on;
xlabel('SNR/dB');
ylabel('MSE(rad^2)');
set(gca, 'YScale', 'log');

plot(SNRdBs, MseESPRIT_LS,'c:s', 'LineWidth', 2);
plot(SNRdBs, MseESPRIT,'r:s', 'LineWidth', 2);
plot(SNRdBs, MseMUSIC,'g:o', 'LineWidth', 2);
plot(SNRdBs, MseCorr,'b:>', 'LineWidth', 2);
plot(SNRdBs, CRLB,'m:<', 'LineWidth', 2);

l = legend({  
        'ESPRIT$_{LS}$',...    
        'ESPRIT$_{TLS}$',...    
        'MUSIC',...
        'Corr',...
        'CRLB'},...    
        'Interpreter','latex', 'Box','off'); 
l.FontSize = 12; 
