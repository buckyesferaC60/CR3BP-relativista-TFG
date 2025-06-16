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

kappa_vec = 1e-30;%[9.8*1e-9]; %7.015*1e-6]; % parámetro de escala: kappa = G(M1+M2)/(L*c^2), siendo M las primarias, L su separación y c la velocidad de la luz
color_vec_vec = [0 0 0;1 0 0;0 0 1;0 1 0];%0 0 0; 1 0 0;

for kappa_loop = 1:length(kappa_vec)

kappa = kappa_vec(kappa_loop);

Orbit_struct.kappa = kappa;

color_vec = color_vec_vec(kappa_loop,:);

[L,PN_1,H,var_vec] = PN_1_Lagrange_MovEquations(mu,kappa);

abs_tol = 1e-20;

rel_tol = 1e-10;

t_end = 7;

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
 
 
 n_orb_init = 1;
 n_orb_end = 100;

 l = 2;



 %chietazeta0 = -[0.09;0;0;0]; % Condición inicial L1 (-) L2(+) L3(-)
 %step_x = -0.002;
 step_z = 0.002;

 f = 1;



%% CONDICIÓN INICIAL

if n_orb_init ==1
    Orbit_struct.family = strcat('Ly_L',num2str(l));

    Orbit_struct.ci = strcat('Ha_L',num2str(l),'_bif');

    Orb_bif = Read_Orbit(Orbit_struct);
else
end

Orbit_struct.family = strcat('Ha_L',num2str(l));

%p = Read_Orbit(Orbit_struct);
sigma_0 = Orb_bif(1,1:end).';

sigma_0(3) = sigma_0(3) + step_z;

%sigma_0(5) = sigma_0(5) + step*25; %*25 en las



%% PROCESO DE BÚSQUEDA DE HALO


% Bucle para familias
for f = n_orb_init:n_orb_end

    if f==13
        f
    end

sigma_m1 = sigma_0;

sigma_0 = [sigma_0(1:6);Vectorize(STM_i)];

    options = odeset('AbsTol',1e-20,'RelTol',1e-10,'Events',@(t,y) cruce(t,y));

% Bucle corrector

tol_Ha_x = 1e-7;
tol_Ha_z = 1e-7;

% figure(5)
% 
% hold off
% plot(0,0,'.')

while 1
    
    % Simulación
    [t00,p00,tc,pc,ic] = ode45(@(t,y) CR3BP_PN_1_STM_ODE(t,y,PN_1,mu,J), [0 t_end], sigma_0, options);

    % figure(5)
    % 
    % 
    % plot3(p00(:,1),p00(:,2),p00(:,3),"Color",[0 0 1])
    % hold on
    % grid on

    STM = Matrixize(pc(1,7:end),[6 6]);

    if abs(pc(1,4)) < tol_Ha_x
        if abs(pc(1,6)) < tol_Ha_z
            [pc(1,4),pc(1,6),(pc(1,4)^2+pc(1,6)^2)^0.5]
            break;
        end
    end

        %Método sencillo; Más rápido para L1 y L2
     
     % sigma_0(5) = sigma_0(5) - pc(1,4)/STM(4,5);
     aux_M = [STM(4,1) STM(4,5);STM(6,1) STM(6,5)];  
     aux_vec_1 = [pc(1,4);pc(1,6)];
     aux_vec_2 = (aux_M)\aux_vec_1;
     sigma_0(1) = sigma_0(1) - aux_vec_2(1);
     sigma_0(5) = sigma_0(5) - aux_vec_2(2);


    % Método mejorado; Muy bueno para L3

     % dydt = CR3BP_PN_1_ODE(0,p00(end,1:6),PN_1,mu);
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
    
    %[pc(1,4),pc(1,6),(pc(1,4)^2+pc(1,6)^2)^0.5]

end


    % Simulate the last result to see precision

    options = odeset('AbsTol',1e-20,'RelTol',1e-10);
    % %[t00,p00] = ode113(@(t,y) CR3BP_PN_1_ODE(t,y,PN_1,mu), [0 t_end], sigma_0(1:6), options);
    [t00,p00] = ode113(@(t,y) CR3BP_PN_1_STM_ODE(t,y,PN_1,mu,J), [0 2*tc(1)], sigma_0, options);

    % Correct the complete orbit

    t_full = t00(end);
format long
    while 1
        [t00,p00] = ode45(@(t,y) CR3BP_PN_1_STM_ODE(t,y,PN_1,mu,J), [0 t_full], sigma_0, options);

        % p00(end,4)

        figure(6)
        hold on

        plot3(p00(:,1),p00(:,2),p00(:,3),"Color",[1 0 1])
        grid on

        if abs(p00(end,4)) < tol_Ha_x*100
            if abs(p00(end,6)) < tol_Ha_z*100
                [p00(end,4),p00(end,6),(p00(end,4)^2+p00(end,6)^2)^0.5]
                break;
            end
        end

        STM = Matrixize(p00(end,7:end),[6 6]);

        aux_M = [STM(4,1) STM(4,5);STM(6,1) STM(6,5)];  
        aux_vec_1 = [pc(1,4);pc(1,6)];
        aux_vec_2 = (aux_M)\aux_vec_1;
        sigma_0(1) = sigma_0(1) + aux_vec_2(1);
        sigma_0(5) = sigma_0(5) + aux_vec_2(2);

        %[p00(end,4),p00(end,6),(p00(end,4)^2+p00(end,6)^2)^0.5]
    end

    % if abs(p00(end,4)) > 1e-8
    %     p00(end,4)
    % end
format short

    Orbit_struct.ci = sigma_0(1);
    Orbit_struct.ci2 = sigma_0(3);

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

    
    sigma_aux = [sigma_0(1) + sigma_0(1)-sigma_m1(1);0;sigma_0(3) + step_z;0;sigma_0(5) + (sigma_0(5)-sigma_m1(5));0];%*(1 + step/chietazeta0(1))

    sigma_m1 = sigma_0;

    sigma_0 = sigma_aux;

    % chietazeta0 = chietazeta0 * (1 + step/chietazeta0(1));
    
 end % bucle for f = 1:N

%end % bucle for l=1:3

    %axis equal

end


%% FUNCIONES

function [value,isterminal,direction] = cruce(t,y)

  value = y(2); % The value that we want to be zero 
  isterminal = 0;  % Halt integration 
  direction = 1;   % The zero can be approached from either direction
  
end

function [value,isterminal,direction] = cruce_end(t,y)

  value = y(2); % The value that we want to be zero 
  isterminal = 1;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction
  
end

