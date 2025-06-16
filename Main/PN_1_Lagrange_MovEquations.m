function [L,PN_1,H,var_vec] = PN_1_Lagrange_MovEquations(mu,kappa)
%PN_1_LAGRANGE_MOVEQUATIONS Summary of this function goes here
%   Detailed explanation goes here

tol = 0.00000001; % tolerancia de obtención de los puntos

% Vector de variables
syms xS yS zS xpS ypS zpS xppS yppS zppS

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

% Puntos de Lagrange
[L1,L2,L3,L4,L5] = Lpoints(PN_1,mu,var_vec,tol);

L(1).vec = [L1;0;0;0;0;0];
L(2).vec = [L2;0;0;0;0;0];
L(3).vec = [L3;0;0;0;0;0];
L(4).vec = [L4(1);L4(2);0;0;0;0];
L(5).vec = [L5(1);L5(2);0;0;0;0];

end

