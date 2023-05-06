function H = genRIS2BS_H(xR1, yR1,...
                        xRN, yRN,...
                        xB1, yB1,...
                        xBN, yBN,...
                        NRIS, NBS, fc, K)
%% generate the channel between the RIS and the BS                
% Inputs:
%       the coordinate of the head and the tail of the RIS (BS) array
%       the number of RIS (BS) elements
%       fc                              =       carrier frequency
% Outputs£º
%       RIS 2 BS channel
%%
H = zeros(NBS, NRIS);
% generate every coordinate of the linear array
xB = linspace(xB1, xBN, NBS);
yB = linspace(yB1, yBN, NBS);
% xR = linspace(xR1, xRN, NRIS);
% yR = linspace(yR1, yRN, NRIS);
for i_NBS = 1:NBS
    H(i_NBS, :) = gen1toArray(xB(i_NBS), yB(i_NBS),...
                   xR1, yR1,...
                   xRN, yRN,...
                   NRIS, fc, K).';
end

end