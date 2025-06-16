

close all
clear
clc

cd('C:\Matlab\CodigosTFG\PN-1')

%% SYSTEM PARAMETERS

mu = 0.01215;
Orbit_struct.mu = mu;

kappa_vec = [1e-5 1e-7 9.8*1e-9 1.16555*1e-11 0];

kappa_vec = 1.16555*1e-11;

tol_Ly = 1e-8;

% Lagrange points
%[L,~,~] = PN_1_Lagrange_MovEquations(mu,0);

% for i = 1:length(kappa_vec)
%     
% end

%for l = 1:3

l = 3;

Orbit_struct.family = strcat('Ly_L',num2str(l));

n_tot(3) = 440;
n_tot(2) = 370;
n_tot(1) = 440;

step = 1;



color_vec = [0 0 0];


%% PROCESO DE LECTURA Y CORRECCIÓN


for k = 1:length(kappa_vec)

    [L,PN_1,~,var_vec] = PN_1_Lagrange_MovEquations(mu,kappa_vec(k)); 

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

    L(1).x0 = -0.009;
    L(2).x0 = 0.005;
    L(3).x0 = -0.09;

    L(1).step = -0.0004;
    L(2).step = 0.0004;
    L(3).step = -0.002;

    L(1).norb_0 = 6;
    L(2).norb_0 = 10;
    L(3).norb_0 = 6;

    L(1).norb_t = 42;
    L(2).norb_t = 37;
    L(3).norb_t = 42;

    L(1).stepstep = 10;
    L(2).stepstep = 10;
    L(3).stepstep = 10;
    
    kappa = kappa_vec(k);

    Orbit_struct.kappa = kappa;

    t_end = 7;

    % Orbit_struct.ci = L(l).vec(1) + L(l).x0 + (L(l).norb_0-1)*L(l).step;
    % 
    % p = Read_Orbit(Orbit_struct);
    % 
    % sigma_0 = p(1,1:end);

    n_start = 1;

    for f = n_start:L(l).norb_t

        Orbit_struct.ci = L(l).vec(1) + L(l).x0 + ((f-1)*L(l).stepstep + L(l).norb_0 - 1)*L(l).step;

        p = Read_Orbit(Orbit_struct);

        sigma_0 = p(1,1:end);

        options = odeset('AbsTol',1e-20,'RelTol',1e-10,'Events',@(t,y) cruce(t,y));

        t_comp = t_end;

        while 1

            t_comp = t_comp*0.75;

            figure(5)
            hold off
            plot(0,0)

            [t00,p00,tc,pc,ic] = ode89(@(t,y) CR3BP_PN_1_STM_ODE(t,y,PN_1,mu,kappa,J), [0 t_comp], sigma_0, options);

            if f > 0
                figure(5)
                hold on
        
                plot3(p00(:,1),p00(:,2),p00(:,3),"Color",[0 0 1])
                grid on
            end            

            if abs(tc(1)) < 1e-1
                p_comp = pc(2,4);
                STM = Matrixize(pc(2,7:end),[6 6]);
                t_comp = 2*tc(2);
            else
                p_comp = pc(1,4);
                STM = Matrixize(pc(1,7:end),[6 6]);
                t_comp = 2*tc(1);
            end

            p_comp

            if abs(p_comp) < tol_Ly
                break;
            end   

            % if ((l == 3) & f<11)
            % 
            %     dydt = CR3BP_PN_1_ODE(0,pc(1,1:6),PN_1,mu,kappa);
            %     aux_M = [STM(2,5) pc(1,5);STM(4,5) dydt(5)];  
            %     aux_vec_1 = [0;pc(1,4)];
            %     aux_vec_2 = (aux_M)\aux_vec_1;
            %     sigma_0(5) = sigma_0(5) + aux_vec_2(1);
            % 
            % else
            % 
                sigma_0(5) = sigma_0(5) - p_comp/STM(4,5);
            % 
            % end          

        end

        [t00,p00] = ode89(@(t,y) CR3BP_PN_1_STM_ODE(t,y,PN_1,mu,kappa,J), [0 t_comp], sigma_0, options);
        
        figure(6)
        hold on

        plot3(p00(:,1),p00(:,2),p00(:,3),"Color",[1 0 0])
        grid on
        
        p00(end,4)

        Write_Orbit_final(Orbit_struct,[t00 p00])

    end

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
