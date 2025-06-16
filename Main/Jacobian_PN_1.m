%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @Author: Pablo Dalda Díaz
%
% 4º del grado en física e instrumentación espacial
%
% Date: 21/05/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

%% 

% Vector de variables
syms xS yS zS xpS ypS zpS xppS yppS zppS mu kappa

var_vec = [xS yS zS xpS ypS zpS xppS yppS zppS];

% Ecuaciones de movimiento
var_vec_0 = var_vec*0;
var_vec_1 = var_vec*0;

[PN_0,H_0,var_vec_0] = Cl_Mov_Eq(mu);
[PN_1,H_1,var_vec_1] = PN_1_Mov_Eq(mu,kappa);

for i = 1:9 % Cambio a las mismas variables
    for j = 1:3
        PN_0(j).val = subs(PN_0(j).val,var_vec_0(i),var_vec(i));
        PN_1(j).val = subs(PN_1(j).val,var_vec_1(i),var_vec(i));
    end

    H_0 = subs(H_0,var_vec_0(i),var_vec(i));
    H_1 = subs(H_1,var_vec_1(i),var_vec(i));
end

H = H_0 + H_1;

    for i = 1:3

        PN_0(i).Eq = solve(PN_0(i).val,var_vec(6+i));

    end

    for i = 1:3

        PN_1(i).val = subs(...
            subs(...
            subs(...
            PN_1(i).val,var_vec(7),PN_0(1).Eq...
            ),var_vec(8),PN_0(2).Eq...
            ),var_vec(9),PN_0(3).Eq);
        PN_1(i).val = PN_1(i).val + PN_0(i).Eq;
        PN_1(i).func = matlabFunction(PN_1(i).val);
 
    end

% sigma(1).x = x;
% sigma(2).x = y;
% sigma(3).x = z;
% sigma(4).x = vx;
% sigma(5).x = vy;
% sigma(6).x = vz;

F_sigma(1).x = var_vec(4);
F_sigma(2).x = var_vec(5);
F_sigma(3).x = var_vec(6);
F_sigma(4).x = PN_1(1).val;
F_sigma(5).x = PN_1(2).val;
F_sigma(6).x = PN_1(3).val;

for i=1:6
    for j=1:6
        D_sigma_F(i,j).x = diff(F_sigma(i).x,var_vec(j));
    end
end

%% Evaluación para puntos xy estacionarios (z=xp=yp=zp=0)

D_sigma_xy = D_sigma_F;

for n=3:6
    for i=1:6
        for j=1:6
            D_sigma_xy(i,j).x = (subs(D_sigma_xy(i,j).x,var_vec(n),0));
        end
    end
end

D_sigma_xy_M = sym(zeros(6,6));

for i = 1:6
    for j=1:6
        D_sigma_xy_M(i,j) = D_sigma_xy(i,j).x;
    end
end

% Obtenenemos poliomio característico

%XY_pol = Pol_carac(D_sigma_xy_M);


%% Evaluación para puntos colineales estacionarios. (antes + y=0)

D_sigma_xcol_M = sym(zeros(6,6));

for i=1:6
    for j=1:6
        D_sigma_xcol_M(i,j) = (subs(D_sigma_xy(i,j).x,var_vec(2),0));
    end
end

dx_Fx = matlabFunction(D_sigma_xcol_M(4,1));

dvy_Fx = matlabFunction(D_sigma_xcol_M(4,5)-2);

dvx_Fy = matlabFunction(D_sigma_xcol_M(5,4)+2);

dy_Fy = matlabFunction(D_sigma_xcol_M(5,2));

dz_Fz = matlabFunction(D_sigma_xcol_M(6,3));

figure(1)
hold on
grid on
figure(2)
hold on
grid on
figure(3)
hold on
grid on
figure(4)
hold on
grid on
figure(5)
hold on
grid on

x_Vec = -2:0.01:2;

%i = 0.5
%j = 0.35
for i=1:4
    for j=0:0.01:0.5

        figure(1)
        plot(x_Vec,dx_Fx(10.^(-2*i),j,x_Vec),"Color",[i/4 0 2*j]);

        figure(2)
        plot(x_Vec,dy_Fy(10.^(-2*i),j,x_Vec),"Color",[i/4 0 2*j]);

        figure(3)
        plot(x_Vec,dvy_Fx(10.^(-2*i),j,x_Vec),"Color",[i/4 0 2*j]);

        figure(4)
        plot(x_Vec,dvx_Fy(10.^(-2*i),j,x_Vec),"Color",[i/4 0 2*j]);

        figure(5)
        plot(x_Vec,dz_Fz(10.^(-2*i),j,x_Vec),"Color",[i/4 0 2*j]);

    end
end

figure(1)
ylim([-10,10])
figure(2)
ylim([-10,10])
figure(3)
ylim([-10,10])
figure(4)
ylim([-10,10])
figure(5)
ylim([-10,10])

%% Obtención de las ecuaciones

syms B1 B2 B3 B4 A1 A2 A3 A4 eta chi etaP chiP Uxx Uyy Uxvy G1 G2 G3 G4 l1 l2 l3 l4

sigma_dif = [eta;chi;etaP;chiP];

%B_vec = [B1;B2;B3;B4];

M_Bsigma = [1 1 1 1; G1 G2 G3 G4; l1 l2 l3 l4;G1*l1 G2*l2 G3*l3 G4*l4];

B_sol = M_Bsigma\sigma_dif;





















