function [h] = gen1toArray(xp, yp,...
                           x1, y1,...
                           xN, yN,...
                           N, fc, K)
%generate the channel from a point to an array                       

% Inputs:
%       xp,yp is the coordinate of the point target
%       x1, y1 is the coordinate of the head of the array
%       xN, yN is the coordinate of the tail of the array
%       N is the number of the array elements
% Outputs:
%       h is the channel from a point to an array£¨N¡Á1£©
%%
lambda = 3e8/fc;
kc = 2*pi/lambda;
x = linspace(x1, xN, N).';% generate the x-coordinate of the array
y = linspace(y1, yN, N).';% generate the y-coordinate of the array
h = exp(-1j*kc* sqrt((x-xp).^2 + (y-yp).^2))./sqrt((x-xp).^2 + (y-yp).^2)/sqrt(4*pi)...
    .*(10.^(K.*sqrt((x-xp).^2 + (y-yp).^2)/10));

% [b] = genb(theta, r, N, fc);


end