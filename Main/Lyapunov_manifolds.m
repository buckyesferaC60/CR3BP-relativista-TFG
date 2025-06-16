


% Obtain the limit orbits

close all
clear
clc


%% SYSTEM PARAMETERS

mu = 0.01215;
Orbit_struct.mu = mu;

kappa = 1e-30;%9.8*1e-9;
Orbit_struct.kappa = kappa;

col = [0 0 0];
f_col = [0 1 1];

% x0 
% LyL1: -0.009 step = -0.0004
% LyL2: 0.005 step = 0.0004
% LyL3: -0.05 step = -0.004

% for i = 1:7
%     figure(i)
%     grid on
%     hold on
% end


%% OBTAIN MONODROMY MATRIX

% Lagrange points
[L,PN_1,~] = PN_1_Lagrange_MovEquations(mu,kappa);

L(1).x0 = -0.009;
L(2).x0 = 0.005;
L(3).x0 = -0.09;

L(1).step = -0.0004;
L(2).step = 0.0004;
L(3).step = -0.002;

l = 1;
%for  l = 1:3
Orbit_struct.family = strcat('Ly_L',num2str(l));

n = 1;
step = 1;

Orbit_struct.ci = L(l).vec(1) + L(l).x0 + step*(n)*L(l).step;

p = Read_Orbit(Orbit_struct);

figure(1)
plot3(p(:,1),p(:,2),p(:,3),'Color',[0 0 1],'LineWidth',0.1)

Mndrmy = Matrixize(p(end,7:end),[6 6]);

Plot_manifolds(p,Mndrmy,PN_1,mu);

[~,eig_val,nu_re,eigval_list] = Monodromy_eigen(Mndrmy);

