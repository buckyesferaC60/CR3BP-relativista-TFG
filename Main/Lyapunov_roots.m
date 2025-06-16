


% Obtain the limit orbits

close all
clear
clc


%% SYSTEM PARAMETERS

mu = 0.01215;
Orbit_struct.mu = mu;

kappa = 0;
Orbit_struct.kappa = kappa;

col = [0 0 0];
f_col = [0 1 1];

% x0 
% LyL1: -0.009 step = -0.0004
% LyL2: 0.005 step = 0.0004
% LyL3: -0.05 step = -0.004

for i = 1:7
    figure(i)
    grid on
    hold on
end


%% OBTAIN MONODROMY MATRIX

% Lagrange points
[L,PN_1,~] = PN_1_Lagrange_MovEquations(mu,kappa);

L(1).x0 = -0.009;
L(2).x0 = 0.005;
L(3).x0 = -0.09;

L(1).step = -0.0004;
L(2).step = 0.0004;
L(3).step = -0.002;

l = 3;
%for  l = 1:3
Orbit_struct.family = strcat('Ly_L',num2str(l));

n_tot = 420;
step = 1;

for i = 1:n_tot/step

    i

Orbit_struct.ci = L(l).vec(1) + L(l).x0 + step*(i-1)*L(l).step;

p = Read_Orbit(Orbit_struct);

figure(1)
plot3(p(:,1),p(:,2),p(:,3),'Color',[0 0 step*i/n_tot],'LineWidth',0.1)

Mndrmy = Matrixize(p(end,7:end),[6 6]);

%Plot_manifolds(p,Mndrmy,PN_1,mu);

[~,eig_val,nu_re,eigval_list] = Monodromy_eigen(Mndrmy);

% eigval_list

figure(2)
plot(Orbit_struct.ci,nu_re,'o','MarkerSize',10,'Color',[1-step*i/n_tot 0 step*i/n_tot]);
figure(3)
plot((real(eig_val)),imag(eig_val),'.','MarkerSize',10,'Color',[1-step*i/n_tot 0 step*i/n_tot]);
figure(4)
for j = 1:length(eig_val)
    if eigval_list(j) == 2
        plot(Orbit_struct.ci,log(eig_val(j)),'.','MarkerSize',10,'Color',[1-step*i/n_tot 0 step*i/n_tot])
    end
end
figure(5)
for j = 1:length(eig_val)
    if eigval_list(j) == 1
        plot(Orbit_struct.ci,real(eig_val(j)),'.','MarkerSize',10,'Color',[1-step*i/n_tot 0 step*i/n_tot])
    end
end
figure(6)
for j = 1:length(eig_val)
    if eigval_list(j) == 1
        plot(i,imag(eig_val(j)),'.','MarkerSize',10,'Color',[1-step*i/n_tot 0 step*i/n_tot])
    end
end
figure(7)
plot3((i)*ones(size(eig_val)),(real(eig_val)),imag(eig_val),'.','MarkerSize',10,'Color',[1-step*i/n_tot 0 step*i/n_tot])
end