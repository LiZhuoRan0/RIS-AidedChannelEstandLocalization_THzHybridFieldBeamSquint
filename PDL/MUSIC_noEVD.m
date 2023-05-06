function r_est = MUSIC_noEVD(M, Y, fc, B, flag)
%%
% input
%       M               number of subcarrier
%       Y               received signal
%       fc              central carrier frequency
%       B               bandwidth
%% 
f = fc + (-B/2 : B/M :(B/2 - B/M));
k = 2*pi*f./3e8;
Ryy = Y'*Y;
Ryy_1 = Ryy(1:size(Y, 1),:);
Ryy_2 = Ryy(size(Y, 1)+1:end,:);
D = (Ryy_1*Ryy_1')\(Ryy_1*Ryy_2');
G = [D;-eye(size(Y, 2) - size(Y, 1))];
Pai_G = G/(G'*G)*G';
r_ini = 0;
if B == 10e9
    r_final = 3e8/B*M;
else
    r_final = 100;
end
interval = [0.1 0.01 0.001 0.0001];
for j = 1:length(interval)
    r = r_ini:interval(j):r_final;
    P = zeros(1, length(r));
    for i = 1:length(r)
        a = exp(-1j*k*r(i));
        P(i) = 1/norm(a*Pai_G);
    end
    %%   
    if j==2 && flag == 'B'
        [value, ~] = max(P);
        P_1 = abs(P - 0.3*value);
        [~, index] = sort(P_1);% ascend
        index = (min(index(1:6)) + max(index(1:6)))/2;
    else
        [~, index] = max(P);
    end
    r_est  = r_ini+(index-1)*interval(j);
    if (flag == 'B') && (j == 2)
        break;
    end

    r_ini = r_est - 2*interval(j);
    r_final = r_est + 2*interval(j);
end
end
