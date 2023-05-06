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
n = (0:(N-1)).';
a = exp(1j*2*pi*d*theta/lambda.*n);

end