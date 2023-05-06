%% Parameter Generation
N = 512;
fc = 100e9;
lambda_c = 3e8/fc;
M = 3;
B = 10e9;
f_min = 95e9;
f_max = 105e9;
f = linspace(fc-B/2, fc+B/2, M);
lambda = 3e8./f;
% f = fc - B/2 + B/M*(0:M-1);
% f = [f_min, fc, f_max];
d = lambda_c/2;
atheta_c = pi/6;
theta_c = sin(atheta_c);
r_c = 50/3;
% a_c = genSteerVector(theta_c, N, d, lambda_c);
a_c = genb(theta_c, r_c, N, fc, d);
% distance grids and angle grids
theta_interval_1 = 0.02;
theta_interval_2 = 0.002;
theta_1 = -1:theta_interval_1:(theta_c-0.1);
theta_2 = (theta_c-0.1):theta_interval_2:(theta_c+0.1);
theta_3 = (theta_c+0.1):theta_interval_1:1;
theta = [theta_1 theta_2 theta_3];
r = zeros(1, 10);
for i_r = 1:size(r, 2)
    r(i_r) = 50/i_r;
end
r = [50 40 30 25 20 15 10 5];
% r = [20 50 100 200 400 600 800];
color = ['r';'g';'m'];
Cor = zeros(2, size(theta,2), size(r,2));
%%
for i_M = 1:size(f,2)
    for i_theta = 1:size(theta,2)
        for i_r = 1:size(r,2)
%             a = genSteerVector(theta(i_theta), N, d, lambda(i_M))/r(i_r);
            a = genb(theta(i_theta), r(i_r), N, f(i_M), d)/r(i_r);
%             a = genb(theta(i_theta), r(i_r), N, f(i_M), d);
            Cor(i_M, i_theta, i_r) = abs(a_c'*a);
        end
    end
    fprintf('i_M = %d\n', i_M);
end
%%
% normalization
Cor = Cor/max(max(max(Cor)));

x = r'*sin((asin(theta)));
y = r'*cos((asin(theta)));
figure; hold on; box on; grid on; 

for i_M = 1:M
%     subplot(2,1,2);hold on;
    mesh(x, y, squeeze(Cor(i_M, :, :))')    
end
set(0,'defaultfigurecolor','w') 
set(gca,'FontName','Times New Roman','FontSize',15);
xlabel('x/m')
ylabel('y/m')
zlabel('Normalized magnitude')
% xlim([-50 50]); ylim([0 50])
set(gca,'DataAspectRatio',[1 1 0.03])
colormap('hot')
colorbar
% shading interp