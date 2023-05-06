function a = genSteerVector(theta, N, d, lambda)
%%  generate the far-field steering vertor

% inputs
%       theta                   =       virtual angle
%       N                       =       number of antennal elements
%       d                       =       spacing of the antenna elements
%       lambda                  =       wavelength

% outputs
%       a                       =       far-field steering vertor£¨N¡Á1£©
%%
n = ((2*(0:N-1)-N+1)/2).';% central element is the reference point
% n = 1:N;
a = zeros(N, length(theta));
for i = 1:length(theta)
%     a(:, i) = exp(1j*2*pi*1/2*sin(theta(i))*(n - 1)).';% no normalization
%     a(:, i) = exp(1j*2*pi*d/lambda*theta*(n - 1)).'/sqrt(N);% normalization
    a(:, i) = exp(1j*2*pi*d/lambda*theta*n).'/sqrt(N);% normalization
%     a(:, i) = exp(1j*2*pi*1/2*theta*(n - 1)).'/sqrt(N);% normalization
end

end