




% Obtención de bifurcaciones por el método de la bisección

%% Parámtertos iniciales

close all
clear
clc

mu = 0.01215;
Orbit_struct.mu = mu;

kappa = 9.8*1e-9;
Orbit_struct.kappa = kappa;

l = 3;
Orbit_struct.family = strcat('Ly_L',num2str(l));

% Lagrange points and movement equations
[L,PN_1,~,var_vec] = PN_1_Lagrange_MovEquations(mu,kappa);

L(1).x0 = -0.009;
L(2).x0 = 0.005;
L(3).x0 = -0.09;

L(1).step = -0.0004;
L(2).step = 0.0004;
L(3).step = -0.002;

x_i = L(l).vec(1) + L(l).x0 + 313*L(l).step;
x_d = L(l).vec(1) + L(l).x0 + 316*L(l).step;

figure(1)
grid on
hold on
figure(2)
grid on
hold on

    x_m = (x_i+x_d)/2;
    % format long
    %     x_i
    %     x_m
    %     x_d
    % format short
    figure(1)
    plot(0,x_i,'+','Color',[1 0 0])
    plot(0,x_m,'+','Color',[0 1 0])
    plot(0,x_d,'+','Color',[0 0 1])

    Orbit_struct.ci = x_i;
    Orb_i = Read_Orbit(Orbit_struct);
    M_i = Matrixize(Orb_i(end,7:end),[6 6]);
    eigval_i = eig(M_i);

    Orbit_struct.ci = x_d;
    Orb_d = Read_Orbit(Orbit_struct);
    M_d = Matrixize(Orb_d(end,7:end),[6 6]);
    eigval_d = eig(M_d);

for i = 1:30

    % x_m = (x_i+x_d)/2;
    % % format long
    % %     x_i
    % %     x_m
    % %     x_d
    % % format short
    % figure(1)
    % plot(i,x_i,'+','Color',[1 0 0])
    % plot(i,x_m,'+','Color',[0 1 0])
    % plot(i,x_d,'+','Color',[0 0 1])

    Orbit_struct.ci = x_i;
    M_i = Matrixize(Orb_i(end,7:end),[6 6]);
    eigval_i = eig(M_i);

    Orbit_struct.ci = x_d;
    M_d = Matrixize(Orb_d(end,7:end),[6 6]);
    eigval_d = eig(M_d);

    Orbit_struct.ci = x_m;
    Orb_m = f_Lyapunov_continuation_CR3BP_PN_1(mu,kappa,l,L,PN_1,var_vec,Orb_i,x_m-x_i);
    M_m = Matrixize(Orb_m(end,7:end),[6 6]);
    eigval_m = eig(M_m);

    figure(2)
    plot(x_i,(max(imag(eigval_i))),'o','Color',[1 0 0])
    plot(x_m,(max(imag(eigval_m))),'o','Color',[0 1 0])
    plot(x_d,(max(imag(eigval_d))),'o','Color',[0 0 1])    

    if max(imag(eigval_m)) > 0
        x_i = x_m;
        Orb_i = Orb_m;
    else
        x_d = x_m;
        Orb_d = Orb_m;
    end

    x_m = (x_i+x_d)/2;
    % format long
    %     x_i
    %     x_m
    %     x_d
    % format short
    figure(1)
    plot(i,x_i,'+','Color',[1 0 0])
    plot(i,x_m,'+','Color',[0 1 0])
    plot(i,x_d,'+','Color',[0 0 1])

end

Orbit_struct.ci = 'Ha_L3_bif';

Write_Orbit(Orbit_struct,Orb_m);






















