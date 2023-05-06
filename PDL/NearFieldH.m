function [H] = NearFieldH(K, N, theta, r, r_LastHop, alpha, fc, B, M, flag)
%% ���ɽ����ŵ�
% Inputs:
%       N                           =       ������
%       theta                       =       UE���еĽǶȣ�ά��ΪL��1
%       r                           =       UE���еľ��룬ά��ΪL��1
%       alpha                       =       UE��BS��˥��ϵ��
%       fc                          =       �ز�Ƶ��
%       B                           =       ����
%       M                           =       ���ز���

%  Outputs;
%       H                           =       Ƶ���ŵ�
%%
H = zeros(N, M);
L = length(theta);
% ����
f = fc + (-B/2 : B/M :(B/2 - B/M));
% f = fc + (0 : B/M :(B - B/M));
lambda = 3e8./f;
d = 3e8/fc/2;
% f = fc+((0:M-1)-M/2+1)/M*B;    % frequency of different carriers
% k = 2*pi*3e8./f;
k = 2*pi*f./3e8;
% kc = 2*pi*3e8./fc;
% tmp = randn(N, 1) + 1j*randn(N, 1);
% tmp = 'relative';
if lower(flag) ==  'relative'
    for m = 1:M
        for l = 1: L

%               H(:, m) = H(:, m) + sqrt(N/L)*alpha(l)*genb(theta(l), r(l), N, f(m), d); 
%             H(:, m) = H(:, m) + exp(1/2*10^(K/10)*r(:,l))*lambda(m)*alpha(:,l).*genb(theta(l), r_LastHop(l), N, f(m), d); 
            H(:, m) = H(:, m) + lambda(m)*alpha(:,l).*genb(theta(l), r_LastHop(l), N, f(m), d); 

    %         H(:, m) = H(:, m) + sqrt(N/L)*alpha(l)*exp(-1j*k(m)*r(l))*genb(theta(l), r(l), N, f(m), d);
    %         H(:, m) = H(:, m) + sqrt(N/L)*alpha(l)*genb(theta(l), r(l), N, f(m), d);
    %         H(:, m) = H(:, m) + sqrt(N/L)*alpha(l)*exp(-1j*k(m)*r(l))*(tmp);
    %         H(:, m) = H(:, m) + sqrt(N/L)*alpha(l)*exp(-1j*k(m)*r(l))*genb(theta(l), r(l), N, fc, d);
    %         H(:, m) = H(:, m) + sqrt(N/L)*alpha(l)*exp(-1j*k(m)*r(l))*genSteerVector(theta(l), N, d, lambda(m));
        end    
    end
else%'absolute'
    for m = 1:M
        for l = 1: L
           H(:, m) = H(:, m) + lambda(m)*exp(-1j*k(m)*r(l))*alpha(:,l).*genb(theta(l), r_LastHop(l), N, f(m), d);   
%            H(:, m) = H(:, m) + exp(-1j*k(m)*r(l))*alpha(:,l).*genb(theta(l), r_LastHop(l), N, f(m), d);   

    %         H(:, m) = H(:, m) + sqrt(N/L)*alpha(l)*exp(-1j*k(m)*r(l))*genb(theta(l), r(l), N, f(m), d);
    %         H(:, m) = H(:, m) + sqrt(N/L)*alpha(l)*genb(theta(l), r(l), N, f(m), d);
    %         H(:, m) = H(:, m) + sqrt(N/L)*alpha(l)*exp(-1j*k(m)*r(l))*(tmp);
    %         H(:, m) = H(:, m) + sqrt(N/L)*alpha(l)*exp(-1j*k(m)*r(l))*genb(theta(l), r(l), N, fc, d);
    %         H(:, m) = H(:, m) + sqrt(N/L)*alpha(l)*exp(-1j*k(m)*r(l))*genSteerVector(theta(l), N, d, lambda(m));
        end    
    end
    
end
end