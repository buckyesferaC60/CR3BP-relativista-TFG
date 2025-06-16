

close all
clear
clc

syms muS kappaS

tol = 0.00000001; % tolerancia de obtención de los puntos

% Vector de variables
syms xS yS zS xpS ypS zpS xppS yppS zppS dxS dyS

var_vec = [xS yS zS xpS ypS zpS xppS yppS zppS];

% Ecuaciones de movimiento
var_vec_0 = var_vec*0;
var_vec_1 = var_vec*0;

[PN_0,H_0,var_vec_0] = Cl_Mov_Eq(muS);
[PN_1,H_1,var_vec_1] = PN_1_Mov_Eq(muS,kappaS);

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
        PN_1(i).func = matlabFunction(PN_1(i).val);%Orden: kappa mu x xp y yp z zp
 
    end

% Ecuaciones de movimiento
F(4).x = PN_1(1).val;
F(5).x = PN_1(2).val;
F(6).x = PN_1(3).val;

dx_F4 = diff(F(4).x,var_vec(1));
dy_F4 = diff(F(4).x,var_vec(2));
dx_F5 = diff(F(5).x,var_vec(1));
dy_F5 = diff(F(5).x,var_vec(2));
dvx_F4 = diff(F(4).x,var_vec(4));
dvy_F4 = diff(F(4).x,var_vec(5));
dvx_F5 = diff(F(5).x,var_vec(4));
dvy_F5 = diff(F(5).x,var_vec(5));

% Jacobiano
J = [0 0 1 0;0 0 0 1;dx_F4 dy_F4 dvx_F4 dvy_F4;dx_F5 dy_F5 dvx_F5 dvy_F5];

J_subs = subs(subs(subs(subs(J,var_vec(3),0),var_vec(4),0),var_vec(5),0),var_vec(6),0);

pretty(simplify(J_subs(3,3)+J_subs(4,4)))

J_func = @(x,y,mu,kappa) subs(subs(subs(subs(J_subs,muS,mu),kappaS,kappa),var_vec(1),x),var_vec(2),y);

%eig_val = eig(J_subs);
mu_vec = 0.01:0.01:0.5;
kappa_vec = [1e-5 1e-7 1e-9 1e-11];

L = zeros(5,2,length(mu_vec),length(kappa_vec));
L_dif = L;



for i = 1:length(mu_vec)
    mu = mu_vec(i)

    for k = 1:3
        PN_0_aux(k).val = subs(PN_0(k).val,muS,mu);
    end

    [L1_cl,L2_cl,L3_cl,~]=Lpoints(PN_0_aux,mu,var_vec,tol);
    
    for j = 1:length(kappa_vec)
        kappa = kappa_vec(j)
        
        for k = 1:3
            PN_1_aux(k).val = subs(subs(PN_1(k).val,kappaS,kappa),muS,mu);
        end

        [L1,L2,L3,L4,L5] = Lpoints(PN_1_aux,mu,var_vec,tol);
        L(1,:,i,j)=[L1;0];
        L_dif(1,:,i,j)=[L1-L1_cl;0];
        L(2,:,i,j)=[L2;0];
        L_dif(2,:,i,j)=[L2-L2_cl;0];
        L(3,:,i,j)=[L3;0];
        L_dif(3,:,i,j)=[L3-L3_cl;0];        
        L(4,:,i,j)=L4;
        L_dif(4,:,i,j)=L4-[0.5-mu;3^0.5/2];
        L(5,:,i,j)=L5;
        L_dif(5,:,i,j)=L5-[0.5-mu;-3^0.5/2];
        figure(1) % Diferencia en x
        hold on
        grid on
        plot(L_dif(4,1,i,j),mu_vec(i),'.','MarkerSize',10,'Color',[kappa_vec(j) 0 mu_vec(i)])
        figure(2) % Diferencia en x
        hold on
        grid on
        plot(L_dif(4,2,i,j),mu_vec(i),'.','MarkerSize',10,'Color',[kappa_vec(j) 0 mu_vec(i)])
        figure(3) % Diferencia en x L1 L2 L3
        hold on
        grid on
        plot(L_dif(1,1,i,j),mu_vec(i),'o','MarkerSize',10,'Color',[1 0 0])
        plot(L_dif(2,1,i,j),mu_vec(i),'+','MarkerSize',10,'Color',[0 0 1])
        plot(L_dif(3,1,i,j),mu_vec(i),'*','MarkerSize',10,'Color',[0 0 0])

    end
end


%% Comparación L4 con la posición clásica


% plot(log(L_dif(4,1,:,:)),log(L_dif(4,2,:,:)),'Color',[0 0 1])
% plot(-log(-L_dif(4,1,:,:)),-log(-L_dif(4,2,:,:)),'Color',[1 0 0])
% figure(2) % Diferencia en y
% hold on
% grid on


%% Observación  de autovalores de J

J_end = zeros(length(mu_vec),length(kappa_vec),4,4);
eig_supervec = zeros(length(mu_vec),length(kappa_vec),4);
Tr_supervec = zeros(length(mu_vec),length(kappa_vec));

for i = 1:length(mu_vec)
    for j = 1:length(kappa_vec)
        J_end(i,j,:,:) = double(subs(subs(subs(subs(J_subs,muS,mu_vec(j))),var_vec(1),L(4,1,i,j)),var_vec(2),L(4,2,i,j)));
        
        eig_supervec(i,j,:) = eig(J_end(i,j,:,:));
        figure(4) % autov de J
        hold on
        grid on
        plot3(mu_vec(i)*ones(4),log(kappa_vec(j))*ones(4),eig_supervec(i,j,:),'.','MarkerSize',10,'Color',[0 0 0])
        Tr_supervec(i,j) = sum(diag(J_end));
        figure(5) % traza de J
        hold on
        grid on
        plot3(mu_vec(i),log(kappa_vec(j)),Tr_supervec(i,j),'.','MarkerSize',10,'Color',[0 0 0])
    end
end





