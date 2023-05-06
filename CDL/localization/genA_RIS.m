function SV = genA_RIS(theta_c, rc, fc, B, N, NRF, i_P, P_RIS)
% design the BS combiner according to the angle, distance from the RIS to
% the BS, center carrier frequency, and the bandwidth
%%
d = 3e8/fc/2;
f = fc - B/2 +B/(NRF*P_RIS)*((i_P-1)*NRF+(1:NRF));
theta = theta_c*fc./f;
SV = zeros(NRF, N);
for i = 1:NRF
%     SV(i,:) = genb(theta(i), rc, N, f(i), d).';
    SV(i,:) = genb(theta_c, rc, N, f(i), d).';
end
end