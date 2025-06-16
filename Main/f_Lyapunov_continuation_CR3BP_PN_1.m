%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @Author: Pablo Dalda Díaz
%
% 4º del grado en física e instrumentación espacial
%
% Date: 21/05/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Base functions and parameters of the CR3BP-Schw

function Orb = f_Lyapunov_continuation_CR3BP_PN_1(mu,kappa,l,L,PN_1,var_vec,Orb_pre,step) 

Orbit_struct.mu = mu;
Orbit_struct.kappa = kappa;
Orbit_struct.family = strcat('Ly_L',num2str(l));

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

% Condición inicial L1 (-) L2(+) L3(-)
%step = -0.002;


%% CONDICIÓN INICIAL

Orbit_struct.ci = Orb_pre(1,1);

p = Orb_pre;
sigma_0 = p(1,1:end).';

sigma_0(1) = sigma_0(1) + step;

%sigma_0(5) = sigma_0(5) + step*50; %*25 en las

%% PROCESO DE BÚSQUEDA DE LYAPUNOV

sigma_0 = [sigma_0(1:6);Vectorize(STM_i)];

options = odeset('AbsTol',1e-20,'RelTol',1e-10,'Events',@(t,y) cruce(t,y));

% Bucle corrector

tol_Ly = 1e-7;

while 1
    
    % Simulación
    [t00,p00,tc,pc,ic] = ode89(@(t,y) CR3BP_PN_1_STM_ODE(t,y,PN_1,mu,J), [0 t_end], sigma_0, options);

    % figure(6)
    % hold on
    % 
    % plot3(p00(:,1),p00(:,2),p00(:,3),"Color",[0 0 1])
    % grid on

    STM = Matrixize(pc(1,7:end),[6 6]);

    if abs(pc(1,4)) < tol_Ly
        break;
    end

        %Método sencillo; Más rápido para L1 y L2
    sigma_0(5) = sigma_0(5) - pc(1,4)/STM(4,5);

    % Método mejorado; Muy bueno para L3

     % dydt = CR3BP_PN_1_ODE(0,p00(end,1:6),PN_1,mu);
     % aux_M = [STM(2,5) pc(1,5);STM(4,5) dydt(5)];  
     % aux_vec_1 = [0;pc(1,4)];
     % aux_vec_2 = (aux_M)\aux_vec_1;
     % sigma_0(5) = sigma_0(5) - aux_vec_2(1);
    %t_end = t_end - aux_vec_2(2);

    % pc(1,4)

end

    fprintf('---\n')

    % Simulate the last result to see precision

    options = odeset('AbsTol',1e-20,'RelTol',1e-10);
    % %[t00,p00] = ode113(@(t,y) CR3BP_PN_1_ODE(t,y,PN_1,mu), [0 t_end], sigma_0(1:6), options);
    [t00,p00] = ode113(@(t,y) CR3BP_PN_1_STM_ODE(t,y,PN_1,mu,J), [0 2*tc(1)], sigma_0, options);

    % Correct the complete orbit

    t_full = t00(end);

    while 1
        [t00,p00] = ode89(@(t,y) CR3BP_PN_1_STM_ODE(t,y,PN_1,mu,J), [0 t_full], sigma_0, options);

        %p00(end,4)

        figure(6)
        hold on

        plot3(p00(:,1),p00(:,2),p00(:,3),"Color",[1 0 1])
        grid on

        if abs(p00(end,4)) < tol_Ly*1000
            break;
        end

        STM = Matrixize(pc(1,7:end),[6 6]);

        sigma_0(5) = sigma_0(5) - pc(1,4)/STM(4,5);

        % dydt = CR3BP_PN_1_ODE(0,p00(end,1:6),PN_1,mu);
        % aux_M = [STM(2,5) p00(end,5);STM(4,5) dydt(5)];  
        % aux_vec_1 = [p00(end,2);p00(end,4)];
        % aux_vec_2 = (aux_M)\aux_vec_1;
        % sigma_0(5) = sigma_0(5) + aux_vec_2(1);
        %t_full = t00(end) - aux_vec_2(2);

        % p00(end,4)
    end

    Orbit_struct.ci = sigma_0(1);

    Write_Orbit(Orbit_struct,p00);

    Orb = p00;

end


%% FUNCIONES

function [value,isterminal,direction] = cruce(t,y)

  value = y(2); % The value that we want to be zero 
  isterminal = 0;  % Halt integration 
  direction = -1;   % The zero can be approached from either direction
  
end

function [value,isterminal,direction] = cruce_end(t,y)

  value = y(2); % The value that we want to be zero 
  isterminal = 1;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction
  
end

