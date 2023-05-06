function Scatter_UE2BS = GenScatterXY(N, UE, BS, num_scatter_UE2BS)
% Thre = eps;
Thre = pi/N*10;
theta_UE2BS = atan((UE.y-BS.y)/(UE.x-BS.x));
Scatter_UE2BS.x = zeros(1, num_scatter_UE2BS);
Scatter_UE2BS.y = zeros(1, num_scatter_UE2BS);
for i = 1:num_scatter_UE2BS
    Scatter_UE2BS.x(i) = 5+(20*sqrt(2)-5)*rand;
    Scatter_UE2BS.y(i) = -30+(-5-(-30))*rand;
    while abs(atan((Scatter_UE2BS.y(i)-BS.y)/(Scatter_UE2BS.x(i)-BS.x)) - theta_UE2BS) ...
            < Thre
        Scatter_UE2BS.x(i) = 5+(20*sqrt(2)-5)*rand;
        Scatter_UE2BS.y(i) = -30+(-5-(-30))*rand;
    end
end
end