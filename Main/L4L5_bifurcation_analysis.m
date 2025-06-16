

clear 
close all
clc

%%

syms A1 A2 A3 A4 B1 B2 B3 B4

M = [0 0 1 0; 0 0 0 1; A1 A2 A3 A4; B1 B2 B3 B4];

pretty(simplify(Pol_carac(M)))

eig_vec_4 = eig(M);

eig_vec_4_func = matlabFunction(eig_vec_4);

%%

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

% for i = 1:4
%     for j = 1:4
%         J_func(i,j).val = subs(subs(subs(subs(J(i,j),var_vec(3),0),var_vec(4),0),var_vec(5),0),var_vec(6),0);
%     end
% end

pretty(simplify(J_subs(3,3)+J_subs(4,4)))

J_func = matlabFunction(J_subs); % kappa mu x y

Tr_J = @(kappa,mu,x,y) sum(diag(J_func(kappa,mu,x,y)));

%eig_val = eig(J_subs);
mu_vec = 0.03851:0.0000001:0.03853;
aux2 = -12:0.1:-5;
kappa_vec = [0 10.^aux2];

L = zeros(5,2,length(mu_vec),length(kappa_vec));
L_dif = L;

%Funciones de posición de L4/L5:
x_L45 = @(mu,kappa) 0.5-mu+kappa*(-1.205*mu+0.62504);
y_L4 = @(mu,kappa) 3^0.5/2 + kappa*(0.43307*mu*(1-mu)-0.36087);
J_L4 = @(mu,kappa) J_func(kappa,mu,x_L45(mu,kappa),y_L4(mu,kappa));

%% 



figure(1)
hold on
grid on
title('Bifurcación de estabilidad L4')
xlabel('∆μ')
ylabel('κ')

for k = 1:length(kappa_vec)

kappa = kappa_vec(k);

mu_i = 0;
mu_d = 0.5;



for i=1:15

    mu_vec_2 = mu_i:(mu_d-mu_i)/500:mu_d;

    for j = 1:length(mu_vec_2)
        aux = eig(J_L4(mu_vec_2(j),kappa));
        J_eig_vec(j) = max(abs(real(aux)));
        %plot(max(abs(real(aux))),mu_vec_2(j)*ones(size(real(aux))),'.','MarkerSize',20)
    end
    %J_eig_vec
    
    for j = 1:length(mu_vec_2)-1
        Dif_eig_vec(j) = abs(J_eig_vec(j+1)-J_eig_vec(j));
    end
    %Dif_eig_vec

    Max_place = find(Dif_eig_vec==max(Dif_eig_vec));

    if Max_place>=11
        if Max_place<=489
            mu_i = mu_vec_2(Max_place-10);
            mu_d = mu_vec_2(Max_place+11);
        else
            mu_i = mu_vec_2(Max_place-10);
            mu_d = mu_vec_2(end);
        end
    else
        mu_i = mu_vec_2(1);   
        mu_d = mu_vec_2(Max_place+11);
    end

    if kappa == 0
        cl_mu = (mu_i+mu_d)/2;
    end

end

    plot((mu_i+mu_d)/2-cl_mu,kappa_vec(k),'.','MarkerSize',10,'Color',[0 0 1])

end





% for i = 1:9
%     figure(i)
%     hold on
%     grid on
% end
% 
% eig_val = zeros(length(mu_vec),length(kappa_vec)*4);
% eig_val_dif = zeros(length(mu_vec),length(kappa_vec)*4);
% 
% for i = 1:length(mu_vec)
%     mu = mu_vec(i);
% 
%     eig_val_0 = eig(J_L4(mu,0));
%     for j = 1:length(kappa_vec)
%         kappa = kappa_vec(j);
% 
%         M = J_L4(mu,kappa);
%         eig_val(i,4*j-3:4*j) = eig(M);
%         eig_val_dif(i,4*j-3:4*j) = eig(M) - eig_val_0;
%     end
% end
% 
% figure(1)
% for i = 1:length(mu_vec)
%     mu = mu_vec(i);
%     plot(real(eig_val(i,:)),mu*ones(size(eig_val(i,:))),'.','MarkerSize',3,'Color',[0, 0, 0])
% end
% 
% figure(2)
% for i = 1:length(mu_vec)
%     mu = mu_vec(i);
%     plot(imag(eig_val(i,:)),mu*ones(size(eig_val(i,:))),'.','MarkerSize',3,'Color',[0, 0, 0])
% end
% 
% figure(3)
% for i = 1:length(mu_vec)
%     mu = mu_vec(i);
%     plot(real(eig_val(i,:)),imag(eig_val(i,:)),'.','MarkerSize',3,'Color',[0, 0, 0])
% end
% 
% figure(4)
% for i = 1:length(mu_vec)
%     mu = mu_vec(i);
%     plot(log(abs(real(eig_val(i,:)))),mu*ones(size(eig_val(i,:))),'.','MarkerSize',3,'Color',[0, 0, 0])
% end
% 
% figure(5)
% for i = 1:length(mu_vec)
%     mu = mu_vec(i);
%     plot(log(abs(imag(eig_val(i,:)))),mu*ones(size(eig_val(i,:))),'.','MarkerSize',3,'Color',[0, 0, 0])
% end
% 
% figure(6)
% for i = 1:length(mu_vec)
%     mu = mu_vec(i);
%     plot(log(abs(real(eig_val(i,:)))),log(abs(imag(eig_val(i,:)))),'.','MarkerSize',3,'Color',[0, 0, 0])
% end
% 
% figure(7)
% for i = 1:length(mu_vec)
%     mu = mu_vec(i);
%     plot(log(abs(real(eig_val_dif(i,:)))),mu*ones(size(eig_val(i,:))),'.','MarkerSize',3,'Color',[0, 0, 0])
% end
% 
% figure(8)
% for i = 1:length(mu_vec)
%     mu = mu_vec(i);
%     plot(log(abs(imag(eig_val_dif(i,:)))),mu*ones(size(eig_val(i,:))),'.','MarkerSize',3,'Color',[0, 0, 0])
% end
% 
% figure(9)
% for i = 1:length(mu_vec)
%     mu = mu_vec(i);
%     plot(log(abs(real(eig_val_dif(i,:)))),log(abs(imag(eig_val_dif(i,:)))),'.','MarkerSize',3,'Color',[0, 0, 0])
% end