


%% Cleaning chores

clear 
close all
clc

%% Parameters

syms x y z vx vy vz ax ay az

var_vec = [x y z vx vy vz ax ay az];


G = 6.67*1e-11; % SI
c = 3*1e8; % SI
M1 = 8*1.989*1e30; % SI
M2 = 2.4*1.989*1e30; %(masa menor de las dos)
r12 = (G*(M1+M2)*(4.8*3600)^2/(2*pi)^2)^(1/3);  % SI


mu = M2/(M1+M2);
kappa = G*(M1+M2)/r12/c^2;

fprintf("\n____________________\nSystem parameters:\n_________________\n")
fprintf("\nM1 = %f Kg\n\nM2 = %f Kg\n\nr12 = %f m",M1,M2,r12)
fprintf("\n------------------\n")
fprintf("\nmu = %.20f\n\nvrot = %.20f c\n\nkappa = %.20f\n\n",mu,kappa^0.5,kappa)

%mu = 0.01542; %0.01542 L-T
%kappa = 0.;

L_points_tol = 0.000000000000001;

    figure(1)
    hold on
    grid on

    plot(-mu,0,'o','Color',[0 0 0])
    plot(1-mu,0,'o','Color',[0 0 0])

%% Classical Lagrange points

    [PN_0,~] = Cl_Mov_Eq(mu);

    [L1,L2,L3,L4,L5]=Lpoints(PN_0,mu,var_vec,L_points_tol);

    plot(L1,0,'+','Color',[0 0 0])
    plot(L2,0,'+','Color',[0 0 0])
    plot(L3,0,'+','Color',[0 0 0])
    plot(L4(1),L4(2),'+','Color',[0 0 0])
    plot(L5(1),L5(2),'+','Color',[0 0 0])

%% PN-1 Lagrange points (According to Hamiltonian)

    [PN_1_H,~] = PN_1_Mov_Eq(mu,kappa);

    PN_1_L = PN_1_H;

    for i = 1:3

        PN_1_H(i).val = PN_1_H(i).val + PN_0(i).val;

    end

    [L1_1,L2_2,L3_3,L4_4,L5_5]=Lpoints(PN_1_H,mu,var_vec,L_points_tol);
    
    figure(1)
    hold on

    plot(L1_1,0,'+','Color',[0 0 1])
    plot(L2_2,0,'+','Color',[0 0 1])
    plot(L3_3,0,'+','Color',[0 0 1])
    plot(L4_4(1),L4_4(2),'+','Color',[0 0 1])
    plot(L5_5(1),L5_5(2),'+','Color',[0 0 1])

%% PN-1 Lagrange points (According to Lagrangian)

    for i = 1:3

        PN_1_L(i).val = subs(subs(subs(PN_1_L(i).val,var_vec(7),PN_0(1).val),var_vec(8),PN_0(2).val),var_vec(9),PN_0(3).val) + PN_0(i).val;

    end

    [L1_1,L2_2,L3_3,L4_4,L5_5]=Lpoints(PN_1_L,mu,var_vec,L_points_tol);
    
    figure(1)
    hold on

    plot(L1_1,0,'+','Color',[1 0 0])
    plot(L2_2,0,'+','Color',[1 0 0])
    plot(L3_3,0,'+','Color',[1 0 0])
    plot(L4_4(1),L4_4(2),'+','Color',[1 0 0])
    plot(L5_5(1),L5_5(2),'+','Color',[1 0 0])

    xlim([-2 2])
    ylim([-1.5 1.5])



