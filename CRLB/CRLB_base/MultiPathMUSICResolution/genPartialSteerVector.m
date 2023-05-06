function a = genPartialSteerVector(theta, N, d, lambda, flag)
%%  generate the partial derivative of steering vector

% inputs
%       theta                   =       virtual angle (after sin)
%       N                       =       number of antenna elements
%       d                       =       elements spacing
%       lambda                  =       wavelength

% outputs
%       a                       =       the partial derivative of steering vector£¨N¡Á1£©
%%
n = (0:(N-1)).';
if flag == 1    
    a = (1j*pi.*n*cos(asin(theta))).*exp(1j*2*pi*d*theta/lambda*n);
%     a = (1j*n).*exp(1j*2*pi*d*theta/lambda*n);
else
    a = -(pi.*n.*cos(asin(theta))).^2.*exp(1j*2*pi*d*theta/lambda*n) -...
        (1j*pi.*n*theta).*exp(1j*2*pi*d*theta/lambda*n);
end
end