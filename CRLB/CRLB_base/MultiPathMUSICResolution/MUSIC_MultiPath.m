function centerFnl = MUSIC_MultiPath(Y, NumSrs)
[N, ~]  = size(Y);
thre    = 1e-18;
Nit     = 10000;
[Ev, D] = eig(Y*Y');
centerFnl = zeros(1, NumSrs);
%% Obtain initial coarse estimations
center = 0;
range   = 1;
En = Ev(:,1:(end-NumSrs));
% Obtain initial spectrum
theta_1     = linspace(center - range, center + range, Nit);
P = zeros(1, Nit);
for i = 1:Nit
    a   = exp(1j*pi*theta_1(i)*(0:(N-1)));
    P(i)    = 1/norm(conj(a)*En);
end
% find the local peak
Pmin = min(P);
Pmax = max(P);
HyperPara = 12;% hyper-parameter
PP = P - ((Pmax - Pmin)/HyperPara + Pmin);
PP(PP<0) = 0;
DiffP = diff(PP);
locNeg = find(DiffP<0);
locPos = find(DiffP>0);
loc = [];
flag = 0;
for i_loc = 1:length(locNeg)   
    if flag == 0
        loc = [loc locNeg(i_loc)];
    end
    if i_loc < length(locNeg)
        if locNeg(i_loc+1) - locNeg(i_loc) == 1
            flag = 1;
        else 
            flag = 0;
        end
    end
end
for i = 1:length(loc)
    loc(i)  = center - range + 2*range/(Nit-1)*(loc(i)-1);
end
%% Angle estimation for each incoming wave
for iSrs = 1:NumSrs
    center = loc(iSrs);
    range   = 1/Nit;
    while range > thre
        theta_1     = linspace(center - range, center + range, Nit);
        P = zeros(1, Nit);
        for i = 1:Nit
            a   = exp(1j*pi*theta_1(i)*(0:(N-1)));
            P(i)    = 1/norm(conj(a)*En);
        end
        [~, index]  = max(abs(P));
        center  = center - range + 2*range/(Nit-1)*(index-1);
        range   = range/Nit;
    end
    centerFnl(iSrs) = center;
end
end