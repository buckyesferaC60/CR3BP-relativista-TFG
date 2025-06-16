%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @Author: Pablo Dalda Díaz
%
% 4º del grado en física e instrumentación espacial
%
% Date: 21/05/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc;

%% Base functions and parameters of the CR3BP-Schw

syms x y z vx vy vz

mu = 0.01215;  % masa reducida: 0.5 - primarias iguales, 1 - una primaria (no usar)  

Orbit_struct.mu = mu;

kappa_vec = 0;%1e-40;%1.16555e-11;%1e-5;%1.16555*1e-11;%[1e-30]; %7.015*1e-6 9.8*1e-9]; % parámetro de escala: kappa = G(M1+M2)/(L*c^2), siendo M las primarias, L su separación y c la velocidad de la luz
color_vec_vec = [0 0 0;1 0 0;0 0 1;0 1 0];%0 0 0; 1 0 0;

for kappa_loop = 1:length(kappa_vec)

kappa = kappa_vec(kappa_loop);

Orbit_struct.kappa = kappa;

color_vec = color_vec_vec(kappa_loop,:);

tol = 0.00000001; % tolerancia de obtención de los puntos

abs_tol = 1e-20;

rel_tol = 1e-10;

t_end = 7;

% Vector de variables
syms xS yS zS xpS ypS zpS xppS yppS zppS

var_vec = [xS yS zS xpS ypS zpS xppS yppS zppS];

% Ecuaciones de movimiento
var_vec_0 = var_vec*0;
var_vec_1 = var_vec*0;

[PN_0,H_0,var_vec_0] = Cl_Mov_Eq(mu);
%[PN_1,H_1,var_vec_1] = PN_1_Mov_Eq(mu,kappa);

for i = 1:9 % Cambio a las mismas variables
    for j = 1:3
        PN_0(j).val = subs(PN_0(j).val,var_vec_0(i),var_vec(i));
        %PN_1(j).val = subs(PN_1(j).val,var_vec_1(i),var_vec(i));
    end

    H_0 = subs(H_0,var_vec_0(i),var_vec(i));
    %H_1 = subs(H_1,var_vec_1(i),var_vec(i));
end

H = H_0;

    for i = 1:3

        PN_0(i).val = solve(PN_0(i).val,var_vec(6+i));
        PN_0(i).func = matlabFunction(PN_0(i).val);%1: x y yp z %2: x xp y z %3: x y z

    end

    % for i = 1:3
    % 
    %     PN_1(i).val = subs(...
    %         subs(...
    %         subs(...
    %         PN_1(i).val,var_vec(7),PN_0(1).Eq...
    %         ),var_vec(8),PN_0(2).Eq...
    %         ),var_vec(9),PN_0(3).Eq);
    %     PN_1(i).val = PN_1(i).val + PN_0(i).Eq;
    %     PN_1(i).func = matlabFunction(PN_1(i).val);
    % 
    % end

% Puntos de Lagrange
[L1,L2,L3,L4,L5] = Lpoints(PN_0,mu,var_vec,tol);

L(1).vec = [L1;0;0;0;0;0];
L(2).vec = [L2;0;0;0;0;0];
L(3).vec = [L3;0;0;0;0;0];
L(4).vec = [L4(1);L4(2);0;0;0;0];
L(5).vec = [L5(1);L5(2);0;0;0;0];

% Ecuaciones de movimiento
F(4).x = PN_0(1).val;
F(5).x = PN_0(2).val;
F(6).x = PN_0(3).val;

dx_F4 = diff(F(4).x,var_vec(1));
dy_F5 = diff(F(5).x,var_vec(2));
dvy_F4 = diff(F(4).x,var_vec(5));
dvx_F5 = diff(F(5).x,var_vec(4));
dz_F6 = diff(F(6).x,var_vec(3));



%for l = 1:3 % For each L_N

    l = 2; % For a L_N

    % Para el correcto almacenamiento
    Orbit_struct.family = strcat('Ly_L',num2str(l));

    %% Caso concreto LN

M_xy = [0 0 1 0; 0 0 0 1; dx_F4 0 0 dvy_F4; 0 dy_F5 dvx_F5 0];
M_z = [0 1;dz_F6 0];

    L(l).M_xy.val = M_xy;
    L(l).M_z.val = M_z;

for i = 1:6
    L(l).M_xy.val = subs(L(l).M_xy.val,var_vec(i),L(l).vec(i));
    L(l).M_z.val = subs(L(l).M_z.val,var_vec(i),L(l).vec(i));
end

L(l).M_xy.val = double(L(l).M_xy.val);
L(l).M_z.val = double(L(l).M_z.val);

[L(l).M_xy.eigvec,L(l).M_xy.eigval] = eig(L(l).M_xy.val);
[L(l).M_z.eigvec,L(l).M_z.eigval] = eig(L(l).M_z.val);

L(l).M_xy.eigval = diag(L(l).M_xy.eigval);
L(l).M_z.eigval = diag(L(l).M_z.eigval);

Gamma = @(n) ( L(l).M_xy.eigval(n)^2 - L(l).M_xy.val(3,1) )/( L(l).M_xy.eigval(n) * L(l).M_xy.val(3,4) );

% Obtenemos los autovalores que son reales
j = 1;
k = 1;
for i = 1:4 
    if isreal(L(l).M_xy.eigval(i))
        L(l).M_xy.eigval_realpos(j) = i;
        j = j+1;
    else
        L(l).M_xy.eigval_nonrpos(k) = i;
        k = k+1;
    end
end

    figure(7)
    hold on
    grid on

    plot3(L(l).vec(1),0,0,"+","Color",[1 0 0])
    plot3(-mu,0,0,"O","Color",[0 0 0.5])
    plot3(1-mu,0,0,"O","Color",[0 0.5 0])

% Jacobiano    
F_sigma(1).x = var_vec(4);
F_sigma(2).x = var_vec(5);
F_sigma(3).x = var_vec(6);
F_sigma(4).x = PN_0(1).val;
F_sigma(5).x = PN_0(2).val;
F_sigma(6).x = PN_0(3).val;

for i=1:6
    for j=1:6
        D_sigma_F(i,j) = diff(F_sigma(i).x,var_vec(j));
    end
end

J = matlabFunction(D_sigma_F);

STM_i = eye(6);

n_orb = 10;
chietazeta0 = [0.005;0;0;0]; % Condición inicial L1 (-) L2(+) L3(-)
step = 0.0004;

% % Bucle para familias
% for f = 1:n_orb
f = 1;

% Condiciones iniciales del movimiento

chi0 = chietazeta0(1);
eta0 = chietazeta0(2);
zeta0 = chietazeta0(3);
zetaP0 = chietazeta0(4);

% Obtención de coeficientes. Disponemos los de las variedades estables a 0
% por definición, y por lo tanto, obtenemos:

    % Coeficientes para x,xP
A_vec = [Gamma(L(l).M_xy.eigval_nonrpos(2)) -1;-Gamma(L(l).M_xy.eigval_nonrpos(1)) 1]...
    /(Gamma(L(l).M_xy.eigval_nonrpos(2)) - Gamma(L(l).M_xy.eigval_nonrpos(1)))...
    * [chi0; eta0];

    % Coeficientes para y,yP
for i = 1:2
    B_vec(i) = A_vec(i) * Gamma(L(l).M_xy.eigval_nonrpos(i));
end

chietaP0_vec = [L(l).M_xy.eigval(L(l).M_xy.eigval_nonrpos(1)) L(l).M_xy.eigval(L(l).M_xy.eigval_nonrpos(2));...
    L(l).M_xy.eigval(L(l).M_xy.eigval_nonrpos(1))*Gamma(L(l).M_xy.eigval_nonrpos(1)) L(l).M_xy.eigval(L(l).M_xy.eigval_nonrpos(2))*Gamma(L(l).M_xy.eigval_nonrpos(2))]...
    *A_vec;

chiP0 = chietaP0_vec(1);
etaP0 = chietaP0_vec(2);

    % Coeficientes para z,zP
C_vec = [L(l).M_z.eigval(2) -1;-L(l).M_z.eigval(1) 1]...
    /(L(l).M_z.eigval(2) - L(l).M_z.eigval(1))...
    * [zeta0; zetaP0];

% Ecuaciones de movimiento

Mov_base(1).t = @(t) real(...
    ...
    A_vec(1).*exp(L(l).M_xy.eigval(L(l).M_xy.eigval_nonrpos(1)).*t) + ...
    A_vec(2).*exp(L(l).M_xy.eigval(L(l).M_xy.eigval_nonrpos(2)).*t) ...
    );

Mov_base(4).t = @(t) real(...
    ...
    A_vec(1).*L(l).M_xy.eigval(L(l).M_xy.eigval_nonrpos(1)).*exp(L(l).M_xy.eigval(L(l).M_xy.eigval_nonrpos(1)).*t) + ...
    A_vec(2).*L(l).M_xy.eigval(L(l).M_xy.eigval_nonrpos(2)).*exp(L(l).M_xy.eigval(L(l).M_xy.eigval_nonrpos(2)).*t) ...
    );

Mov_base(2).t = @(t) real(...
    ...
    B_vec(1).*exp(L(l).M_xy.eigval(L(l).M_xy.eigval_nonrpos(1)).*t) + ...
    B_vec(2).*exp(L(l).M_xy.eigval(L(l).M_xy.eigval_nonrpos(2)).*t) ...
    );

Mov_base(5).t = @(t) real(...
    ...
    B_vec(1).*L(l).M_xy.eigval(L(l).M_xy.eigval_nonrpos(1)).*exp(L(l).M_xy.eigval(L(l).M_xy.eigval_nonrpos(1)).*t) + ...
    B_vec(2).*L(l).M_xy.eigval(L(l).M_xy.eigval_nonrpos(2)).*exp(L(l).M_xy.eigval(L(l).M_xy.eigval_nonrpos(2)).*t) ...
    );

Mov_base(3).t = @(t) real(...
    ...
    C_vec(1).*exp(L(l).M_z.eigval(1).*t) + ...
    C_vec(2).*exp(L(l).M_z.eigval(2).*t) ...
    );

Mov_base(6).t = @(t) real(...
    ...
    C_vec(1).*L(l).M_z.eigval(1).*exp(L(l).M_z.eigval(1).*t) + ...
    C_vec(2).*L(l).M_z.eigval(2).*exp(L(l).M_z.eigval(2).*t) ...
    );

t_vec = 0:0.025:t_end;

% Representación gráfica del movimiento en la aproximación

if f == 1

figure(1)
hold on

plot(t_vec,Mov_base(1).t(t_vec),"Color",[1 0 0]);
plot(t_vec,Mov_base(2).t(t_vec),"Color",[0 0 1]);
plot(t_vec,Mov_base(3).t(t_vec),"Color",[0 1 0]);
grid on

figure(2)
hold on

plot3(L(l).vec,0,0,"+","Color",[1 0 0])
plot3(-mu,0,0,"O","Color",[0 0 0.5])
plot3(1-mu,0,0,"O","Color",[0 0.5 0])

plot3(Mov_base(1).t(t_vec)+L(l).vec(1),Mov_base(2).t(t_vec),Mov_base(3).t(t_vec),"Color",[0 0 0])
grid on

figure(3)
hold on

plot3(0,0,0,"+","Color",[1 0 0])

plot3(Mov_base(1).t(t_vec),Mov_base(2).t(t_vec),Mov_base(3).t(t_vec),"Color",[0 0 0])
grid on

figure(4)
hold on

plot(t_vec,Mov_base(4).t(t_vec),"Color",[1 0 0]);
plot(t_vec,Mov_base(5).t(t_vec),"Color",[0 0 1]);
plot(t_vec,Mov_base(6).t(t_vec),"Color",[0 1 0]);
grid on

% figure(5)
% hold on
% 
% plot(t_vec,Mov_base(1).t(t_vec),"Color",[1 0 0]);
% plot(t_vec,Mov_base(2).t(t_vec),"Color",[0 0 1]);
% plot(t_vec,Mov_base(3).t(t_vec),"Color",[0 1 0]);
% grid on
% 
% figure(6)
% hold on
% 
% plot3(L(l).vec(1),0,0,"+","Color",[1 0 0])
% plot3(Mov_base(1).t(t_vec)+L(l).vec(1),Mov_base(2).t(t_vec),Mov_base(3).t(t_vec),"Color",[0 0 0])
% grid on

end

%% TRAYECTORIA REAL CON MISMAS CI

if f == 1
    sigma_0 = [chi0+L(l).vec(1);eta0;zeta0;chiP0;etaP0;zetaP0];
% else
%     sigma_0 = sigma_0;%[chi0+L(l).vec(1);eta0;zeta0;chiP0*10;etaP0;zetaP0];% + sigma_0)/2;
%     %sigma_0(5) = sigma_0(5) + etaP0;
end
% 
% if f == 1
%     options = odeset('AbsTol',1e-15,'RelTol',1e-10);
%     %Integración
%     [t00,p00] = ode113(@(t,y) CR3BP_PN_1_ODE(t,y,PN_1,mu), [0 t_end], sigma_0,options);%,'Events', @(t,y) cruce(t,y,mu));
% 
% figure(1)
% hold on
% 
% plot(t00,p00(:,1)-L(l).vec(1),"Color",[0.5 0 0]);
% plot(t00,p00(:,2),"Color",[0 0 0.5]);
% plot(t00,p00(:,3),"Color",[0 0.5 0]);
% grid on
% 
% figure(2)
% hold on
% 
% plot3(p00(:,1),p00(:,2),p00(:,3),"Color",[0 0 1])
% grid on
% 
% figure(3)
% hold on
% 
% plot3(p00(:,1)-L(l).vec(1),p00(:,2),p00(:,3),"Color",[0 0 1])
% grid on
% 
% figure(4)
% hold on
% 
% plot(t00,p00(:,4)-L(l).vec(1),"Color",[0.5 0 0]);
% plot(t00,p00(:,5),"Color",[0 0 0.5]);
% plot(t00,p00(:,6),"Color",[0 0.5 0]);
% grid on
% 
% end

%% PROCESO DE BÚSQUEDA DE LYAPUNOV


% Bucle para familias
for f = 1:n_orb

sigma_0 = [sigma_0(1:6);Vectorize(STM_i)];

    options = odeset('AbsTol',1e-20,'RelTol',1e-10,'Events',@(t,y) cruce(t,y));

% Bucle corrector

tol_Ly = 1e-8;

for i = 1:50
    
    % Simulación
    [t00,p00,tc,pc,ic] = ode113(@(t,y) CR3BP_PN_1_STM_ODE(t,y,PN_0,mu,kappa,J), [0 t_end], sigma_0, options);
    figure(5)
    hold on

    plot3(p00(:,1),p00(:,2),p00(:,3),"Color",[0 0 1])
    grid on
    

    STM = Matrixize(pc(2,7:end),[6 6]);

    tc
    pc(2,1:6)

    if abs(pc(2,4)) < tol_Ly
        break;
    end

        %Método sencillo; Más rápido para L1 y L2
      sigma_0(5) = sigma_0(5) - pc(2,4)/STM(4,5);

    % Método mejorado; Muy bueno para L3

     % dydt = CR3BP_PN_1_ODE(0,p00(end,1:6),PN_0,mu,kappa);
     % aux_M = [STM(2,5) pc(1,5);STM(4,5) dydt(5)];  
     % aux_vec_1 = [0;pc(1,4)];
     % aux_vec_2 = (aux_M)\aux_vec_1;
     % sigma_0(5) = sigma_0(5) - aux_vec_2(1);
    %t_end = t_end - aux_vec_2(2);

    %sigma_0(5) = sigma_0(5) - p00(end,4)/(STM(4,5) - dydt(5)/p00(end,5)*STM(2,5));

    %r = dydt(5)/p00(end,5)*STM(2,5)/STM(4,5);
    %sigma_0(5) = sigma_0(5) - p00(end,4)/STM(4,5) * (1 + r + r^2);

    % figure(5)
    % hold on
    % 
    % plot(t00,p00(:,1)-L(l).vec(1),"Color",[0.5 0 0]);
    % plot(t00,p00(:,2),"Color",[0 0 0.5]);
    % plot(t00,p00(:,3),"Color",[0 0.5 0]);
    % grid on
    % 
    % figure(6)
    % hold on
    % 
    % plot3(p00(:,1),p00(:,2),p00(:,3),"Color",[0 0 1])
    % grid on

    %pc(1,4)

end


    % Simulate the last result to see precision

    options = odeset('AbsTol',1e-20,'RelTol',1e-10);
    % %[t00,p00] = ode113(@(t,y) CR3BP_PN_1_ODE(t,y,PN_1,mu), [0 t_end], sigma_0(1:6), options);
    [t00,p00] = ode113(@(t,y) CR3BP_PN_1_STM_ODE(t,y,PN_0,mu,kappa,J), [0 2*tc(1)], sigma_0, options);

    % Correct the complete orbit

    t_full = t00(end);

    for i = 1:25
        [t00,p00] = ode89(@(t,y) CR3BP_PN_1_STM_ODE(t,y,PN_0,mu,kappa,J), [0 t_full], sigma_0, options);

        % p00(end,4)

        if abs(p00(end,4)) < tol_Ly*100
            break;
        end

        STM = Matrixize(p00(end,7:end),[6 6]);

        dydt = CR3BP_PN_1_ODE(0,p00(end,1:6),PN_1,mu,0);
        aux_M = [STM(2,5) p00(end,5);STM(4,5) dydt(5)];  
        aux_vec_1 = [p00(end,2);p00(end,4)];
        aux_vec_2 = (aux_M)\aux_vec_1;
        sigma_0(5) = sigma_0(5) - aux_vec_2(1);
        %t_full = t00(end) - aux_vec_2(2);

        % figure(6)
        % hold on
        % 
        % plot3(p00(:,1),p00(:,2),p00(:,3),"Color",[1 0 1])
        % grid on

        % p00(end,4)
    end

    if abs(p00(end,4)) > tol_Ly*100
        p00(end,4)
    end



    Orbit_struct.ci = sigma_0(1);

    Write_Orbit(Orbit_struct,p00);



    % figure(5)
    % hold on
    % 
    % plot(t00,p00(:,1)-L(l).vec(1),"Color",[0.1 0 0]);
    % plot(t00,p00(:,2),"Color",[0 0 0.1]);
    % plot(t00,p00(:,3),"Color",[0 0.1 0]);
    % grid on
    % 
    % figure(6)
    % hold on
    % 
    % plot3(p00(:,1),p00(:,2),p00(:,3),"Color",[1 0 0])
    % grid on

    % True periodic orbit: forward & backward half orbit

        options = odeset('AbsTol',1e-20,'RelTol',1e-10,'Events',@(t,y) cruce_end(t,y));
    [t_f,p_f,tc_f,pc_f,ic_f] = ode113(@(t,y) CR3BP_PN_1_STM_ODE(t,y,PN_0,mu,kappa,J), [0 t_end*2], sigma_0, options);
    [t_b,p_b,tc_b,pc_b,ic_b] = ode113(@(t,y) CR3BP_PN_1_STM_ODE(t,y,PN_0,mu,kappa,J), [0 -t_end*2], sigma_0, options);

    % Guardar órbita en archivo
    

    


    %Monodromy analysis:--------------------------------------------------------------------------------

    % STM_f = Matrixize(p_f(end,7:end),[6 6]);
    % STM_b = Matrixize(p_b(end,7:end),[6 6]);
    % 
    % Monodromy_STM = STM_f * inv(STM_b);
    % 
    % [Mndr_eig_vec,Mndr_eig_val,nu_re,nu_im,Mndr_re_eig_vec,Mndr_im_eig_vec] = Monodromy_eigen(Monodromy_STM);

    %--------------------------------------------------------------------------------------------------

    % figure(6)
    % hold on
    %
    %plot3(p_f(:,1),p_f(:,2),p_f(:,3),"Color",[1 0 1])
    %plot3(p_b(:,1),p_b(:,2),p_b(:,3),"Color",[1 0 1])

    figure(7)
    hold on
    grid on

    if mod(f,2) == 1
        plot3(p_f(:,1),p_f(:,2),p_f(:,3),"Color",color_vec,"LineWidth",1.5)
        plot3(p_b(:,1),p_b(:,2),p_b(:,3),"Color",color_vec,"LineWidth",1.5)
    % else
    %     plot3(p_f(:,1),p_f(:,2),p_f(:,3),"Color",[0.5 0.5 0],"LineWidth",0.25)
    %     plot3(p_b(:,1),p_b(:,2),p_b(:,3),"Color",[0.5 0.5 0],"LineWidth",0.25)
    end

    fprintf('L_point: L%i _________Orbit: %i ___________ic: (%f %f %f %f %f %f)\n',l,f,sigma_0(1),sigma_0(2),sigma_0(3),sigma_0(4),sigma_0(5),sigma_0(6));

    sigma_0 = [(sigma_0(1)-L(l).vec(1))*(1 + step/chietazeta0(1))+L(l).vec(1);0;0;0;sigma_0(5);0];%*(1 + step/chietazeta0(1)) 

    chietazeta0 = chietazeta0 * (1 + step/chietazeta0(1));
    
 end % bucle for f = 1:N

%end % bucle for l=1:3

    axis equal

end


%% FUNCIONES

function [value,isterminal,direction] = cruce(t,y)

  value = y(2); % The value that we want to be zero 
  isterminal = 0;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction
  
end

function [value,isterminal,direction] = cruce_end(t,y)

  value = y(2); % The value that we want to be zero 
  isterminal = 1;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction
  
end

