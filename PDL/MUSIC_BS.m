function theta_est = MUSIC_BS(Y,A)

[Ev, ~] = eig(Y*Y');
En = Ev(:,1:(end-1));
N = size(A,2);
No = 2;
P = zeros(1, N*No);
for i = 1:N*No
        a_tmp = exp(-1j*pi*i*(0:(N-1))/N/No);
        a = A*a_tmp.';
        P(i) = 1/(a.'*En*En'*conj(a));
end
% figure;plot(abs(P))
[~, index ] = max(abs(P));
theta_est = index/N/No;

%%
range   = 0.01;
thre    = 1e-8;
center  = theta_est;
Nit     = 40;
while range > thre
    theta_1     = linspace(center - range, center + range, Nit);
    P = zeros(1, Nit);
    for i = 1:Nit
        a_tmp   = exp(-1j*pi*theta_1(i)*(0:(N-1)));
        a = A*a_tmp.';
        P(i)    = 1/norm(a.'*En);
    end
    [~, index]  = max(abs(P));
    center  = center - range + 2*range/(Nit-1)*(index-1);
    range   = range/10;
end
theta_est = center;
end