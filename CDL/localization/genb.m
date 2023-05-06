function [b] = genb(theta, r, N, f, d)
%% generate the near-field steering vector

% Inputs£º
%       theta & r=          =       angle and distance
%       N                   =       number of antenna elements
%       fc                  =       carrier frequency
%       d                   =       spacing of the antenna element

% Outputs£º
%       b                   =       near-field steering vector
%%
% b = zeros(N, 1);
kc = 2*pi*f/3e8;
delta_n = ((2*(0:N-1)-N+1)/2).';% central element is the reference point
% delta_n = (0:N-1).';% first element is the reference point
% d = lambda/2;
r_l_n = sqrt(r^2 + delta_n.^2.*d.^2 - 2*r*theta*delta_n*d);

b = exp(-1j*kc*(r_l_n-r))/sqrt(N);
end