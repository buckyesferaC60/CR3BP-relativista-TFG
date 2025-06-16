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

kappa_vec = 1.16555*1e-11;%1e-40;%9.8*1e-9;%1e-7;%1e-5;%1.16555*1e-11;%1e-30;%[9.8*1e-9]; %7.015*1e-6]; % parámetro de escala: kappa = G(M1+M2)/(L*c^2), siendo M las primarias, L su separación y c la velocidad de la luz
color_vec_vec = [0 0 0;1 0 0;0 0 1;0 1 0];%0 0 0; 1 0 0;

for kappa_loop = 1:length(kappa_vec)

kappa = kappa_vec(kappa_loop);

Orbit_struct.kappa = kappa;

color_vec = color_vec_vec(kappa_loop,:);

[L,PN_1,H,var_vec] = PN_1_Lagrange_MovEquations(mu,kappa);

abs_tol = 1e-20;

rel_tol = 1e-10;

t_end = 15;

% Jacobiano    
F_sigma(1).x = var_vec(4);
F_sigma(2).x = var_vec(5);
F_sigma(3).x = var_vec(6);
F_sigma(4).x = PN_1(1).val;
F_sigma(5).x = PN_1(2).val;
F_sigma(6).x = PN_1(3).val;

for i=1:6
    for j=1:6
        D_sigma_F(i,j) = diff(F_sigma(i).x,var_vec(j));
    end
end

 J = matlabFunction(D_sigma_F);
 
 STM_i = eye(6);
 
 
 n_orb_init = 400;
 n_orb_end = 440;

 l = 3;

 Orbit_struct.family = strcat('Ly_L',num2str(l));

 chietazeta0 = -[0.05;0;0;0]; % Condición inicial L1 (-) L2(+) L3(-)
 step = -0.002;

 f = 1;

%% CONDICIÓN INICIAL

Orbit_struct.ci = L(l).vec(1) + chietazeta0(1) + (n_orb_init-1)*step;

p = Read_Orbit(Orbit_struct);
sigma_0 = p(1,1:end).';

sigma_0(1) = sigma_0(1) + step;

sigma_0(5) = sigma_0(5) - step*10 ;%- 25*step;%; %*25 en las

sigma_m1 = sigma_0;

%% PROCESO DE BÚSQUEDA DE LYAPUNOV

t_orb = t_end;

% Bucle para familias
for f = n_orb_init:n_orb_end

sigma_0 = [sigma_0(1:6);Vectorize(STM_i)];

    options = odeset('AbsTol',1e-20,'RelTol',1e-10,'Events',@(t,y) cruce(t,y));

% Bucle corrector

tol_Ly = 1e-8;
if f>100
    tol_Ly = 1e-7;
end
if f>150
    tol_Ly = 1e-6;
end
if f>200
    tol_Ly = 1e-5;
end
if f>300
    tol_Ly = 1e-4;
end

    t_end = 0.75*t_orb;

while 1
    figure(5)
    hold off
    plot(0,0)
    
    % Simulación
    [t00,p00,tc,pc,ic] = ode45(@(t,y) CR3BP_PN_1_STM_ODE(t,y,PN_1,mu,kappa,J), [0 t_end], sigma_0, options);

    if f > 0
        figure(5)
        hold on

        plot3(p00(:,1),p00(:,2),p00(:,3),"Color",[0 0 1])
        grid on
    end

    if abs(tc(1)) < 1e-1
        p_comp = pc(2,4);
        STM = Matrixize(pc(2,7:end),[6 6]);
    else
        p_comp = pc(1,4);
        STM = Matrixize(pc(1,7:end),[6 6]);
    end
    
    p_comp
    % pc(:,4)
    % tc

    if abs(p_comp) < tol_Ly
        break;
    end

        %Método sencillo; Más rápido para L1 y L2
     sigma_0(5) = sigma_0(5) - p_comp/STM(4,5);

    %Método mejorado; Muy bueno para L3

     % dydt = CR3BP_PN_1_ODE(0,p00(end,1:6),PN_1,mu);
     % aux_M = [STM(2,5) pc(1,5);STM(4,5) dydt(5)];  
     % aux_vec_1 = [0;pc(1,4)];
     % aux_vec_2 = (aux_M)\aux_vec_1;
     % sigma_0(5) = sigma_0(5) + aux_vec_2(1);
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

end

    if abs(tc(1)) < 1e-1
        t_comp = 2*tc(2);
    else
        t_comp = 2*tc(1);
    end

    % Simulate the last result to see precision

    options = odeset('AbsTol',1e-20,'RelTol',1e-10);
    % %[t00,p00] = ode113(@(t,y) CR3BP_PN_1_ODE(t,y,PN_1,mu), [0 t_end], sigma_0(1:6), options);
    [t00,p00] = ode45(@(t,y) CR3BP_PN_1_STM_ODE(t,y,PN_1,mu,kappa,J), [0 t_comp], sigma_0, options);

    % Correct the complete orbit

    t_full = t00(end);
format long
    while 1
        [t00,p00] = ode45(@(t,y) CR3BP_PN_1_STM_ODE(t,y,PN_1,mu,kappa,J), [0 t_full], sigma_0, options);

        % p00(end,4)

        figure(6)
        hold on

        plot3(p00(:,1),p00(:,2),p00(:,3),"Color",[1 0 1])
        grid on

        if abs(p00(end,4)) < tol_Ly*100
            break;
        end

        STM = Matrixize(pc(1,7:end),[6 6]);

        dydt = CR3BP_PN_1_ODE(0,p00(end,1:6),PN_1,mu);
        aux_M = [STM(2,5) p00(end,5);STM(4,5) dydt(5)];  
        aux_vec_1 = [0;p00(end,4)];
        aux_vec_2 = (aux_M)\aux_vec_1;
        sigma_0(5) = sigma_0(5) + aux_vec_2(1);
        %t_full = t00(end) - aux_vec_2(2);

        p00(end,4)
    end

    % if abs(p00(end,4)) > 1e-8
    %     p00(end,4)
    % end
format short

    Orbit_struct.ci = sigma_0(1);

    Write_Orbit(Orbit_struct,p00);

    t_orb = t00(end);

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

    %     options = odeset('AbsTol',1e-20,'RelTol',1e-10,'Events',@(t,y) cruce_end(t,y));
    % [t_f,p_f,tc_f,pc_f,ic_f] = ode113(@(t,y) CR3BP_PN_1_STM_ODE(t,y,PN_1,mu,J), [0 t_end*2], sigma_0, options);
    % [t_b,p_b,tc_b,pc_b,ic_b] = ode113(@(t,y) CR3BP_PN_1_STM_ODE(t,y,PN_1,mu,J), [0 -t_end*2], sigma_0, options);

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

    % figure(7)
    % hold on
    % grid on
    % 
    % if mod(f,5) == 1
    %     plot3(p_f(:,1),p_f(:,2),p_f(:,3),"Color",color_vec,"LineWidth",1.5)
    %     plot3(p_b(:,1),p_b(:,2),p_b(:,3),"Color",color_vec,"LineWidth",1.5)
    % % else
    % %     plot3(p_f(:,1),p_f(:,2),p_f(:,3),"Color",[0.5 0.5 0],"LineWidth",0.25)
    % %     plot3(p_b(:,1),p_b(:,2),p_b(:,3),"Color",[0.5 0.5 0],"LineWidth",0.25)
    % end

    fprintf('_________L_point: L%i _________Orbit: %i ___________ic: (%f %f %f %f %f %f)\n',l,f,sigma_0(1),sigma_0(2),sigma_0(3),sigma_0(4),sigma_0(5),sigma_0(6));

    
    sigma_aux = [(sigma_0(1)-L(l).vec(1) + step)+L(l).vec(1);0;0;0;sigma_0(5)+(sigma_0(5)-sigma_m1(5));0];%%%*(1 + step/chietazeta0(1))

    sigma_m1 = sigma_0;

    sigma_0 = sigma_aux;

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

