function [H] = NearFieldH(K, N, theta, r, r_LastHop, alpha, fc, B, M, flag)
%% generate the near-field channel
% Inputs:
%       N                           =       number of antennas
%       theta                       =       all the angles of the channel, dimension is L¡Á1
%       r                           =       all the distance of the channel, dimension is L¡Á1
%       alpha                       =       attenuation coefficient
%       fc                          =       center carrier frequency
%       B                           =       bandwidth
%       M                           =       number of carriers

%  Outputs;
%       H                           =       frequency domain channel
%%
H = zeros(N, M);
L = length(theta);
f = fc + (-B/2 : B/M :(B/2 - B/M));
lambda = 3e8./f;
d = 3e8/fc/2;
k = 2*pi*f./3e8;
if lower(flag) ==  'relative'
    for m = 1:M
        for l = 1: L
            H(:, m) = H(:, m) + lambda(m)*alpha(:,l).*genb(theta(l), r_LastHop(l), N, f(m), d); 
        end    
    end
else%'absolute'
    for m = 1:M
        for l = 1: L
           H(:, m) = H(:, m) + lambda(m)*exp(-1j*k(m)*r(l))*alpha(:,l).*genb(theta(l), r_LastHop(l), N, f(m), d);   
        end    
    end
    
end
end