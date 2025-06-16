%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @author: Pablo Dalda Díaz
% Alumno de 4º del grado en FIE
% 18/05/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parámetros del modelo

close all
clear
clc

fprintf('___________________________\nInitializing variables\n___________________________')

dt = 0.1;     %paso temporal
t0 = 0;         %instante inicial (s)
tf = 1000;      %instante final (s)

n_max = 100;
intervalos = 3;

mu = 0.01528;

tol = 0.00000000001;

syms xS yS zS xpS ypS zpS xppS yppS zppS

var_vec = [xS yS zS xpS ypS zpS xppS yppS zppS];

b_values_init = 4;
b_values_end = 4;

for b = b_values_init:b_values_end

 kappa_vec = [0 1e-11 1e-2 1e-6 1e-3 ];
kappa = kappa_vec(b);

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




matlabFunction(H)

    for i = 1:3

        PN_0(i).Eq = solve(PN_0(i).val,var_vec(6+i));

    end

    for i = 1:3

        PN_1(i).val = subs(...
            subs(...
            subs(...
            PN_1(i).val,var_vec(7),PN_0(1).Eq...
            ),var_vec(8),PN_0(2).Eq...
            ),var_vec(9),PN_0(3).Eq...
            ) + PN_0(i).Eq,var_vec(6+i);
        PN_1(i).func = matlabFunction(PN_1(i).val);
 
    end

[L1,L2,L3,L4,L5] = Lpoints(PN_1,mu,var_vec,tol);

L1 = [L1 0 0];
L2 = [L2 0 0];
L3 = [L3 0 0];
L4 = [L4(1) L4(2) 0];
L5 = [L5(1) L5(2) 0];

H_mat = matlabFunction(H);

C_aux.L(1).val = -2*H_mat(L1(1),0,0,0,0,0);
C_aux.L(2).val = -2*H_mat(L2(1),0,0,0,0,0);
C_aux.L(3).val = -2*H_mat(L3(1),0,0,0,0,0);
C_aux.L(4).val = -2*H_mat(L4(1),0,L4(2),0,0,0);
C_aux.L(5).val = -2*H_mat(L5(1),0,L5(2),0,0,0);

%H = -2;
C_vec = C_aux.L(4).val+C_aux.L(3).val; %H=-1.5

for c_aux = 1:length(C_vec)

figure(b+c_aux*(b_values_end-b_values_init+1))
hold on
plot3(L1(1),L1(2),L1(3),'+',"Color",[kappa 0.5 2*mu],'MarkerSize',10);
plot3(L2(1),L2(2),L2(3),'+',"Color",[kappa 0.5 2*mu],'MarkerSize',10);
plot3(L3(1),L3(2),L3(3),'+',"Color",[kappa 0.5 2*mu],'MarkerSize',10);
%plot3(L4(1),L4(2),L4(3),'+',"Color",[kappa 0.5 2*mu],'MarkerSize',10);
%plot3(L5(1),L5(2),L5(3),'+',"Color",[kappa 0.5 2*mu],'MarkerSize',10);
plot3(-mu,0,0,'*',"Color",[0 0 1],'MarkerSize',10);
plot3(1-mu,0,0,'*',"Color",[0 0 1],'MarkerSize',10);

grid on

C = C_vec(c_aux);
Href = -C/2;   

% Y_vec_s = solve(H,var_vec(5));
% 
% Y_vec(1).y = Y_vec_s(1);
% Y_vec(2).y = Y_vec_s(2);

%Y_vec(1).y = matlabFunction(Y_vec(1).y);  %Hs Vx X Y
%Y_vec(2).y = matlabFunction(Y_vec(2).y);

% x(:,1) = linspace(-2,2,n_max);
% x2(:,1) = linspace(-2,2,intervalos);

x(1).init = linspace(-mu-1,-mu-0.1, round(n_max));
x(2).init = linspace(-mu+0.1, 1-mu-0.1, round(n_max));
x(3).init = linspace(1-mu+0.1, 2-mu, round(n_max));

%x(1).init = linspace(-0.3,0,n_max);

y = 0;%linspace(-2,2,n_max);
vx = 0;%linspace(-5,5,n_max);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long

for j = 1:1:intervalos
for i = 1:1:length(x(j).init)
    fprintf('\n\n_______________Condición número: %i\n',i)
for k = 1:4 % Cambiar en función de la solución

    x_cr_v1 = NaN;
    x_cr_v2 = NaN;
    vx_cr_v1 = NaN;
    vx_cr_v2 = NaN;
    vy_cr_v = NaN;

for w = 1:1:length(y)
for g = 1:1:length(vx)


x00 = x(j).init(i);
y00 = y(w);
z00 = 0;

vx00 = vx(g);
vy00 = var_vec(5);
vz00 = 0;

    i0 = [x00; y00; z00; vx00; vy00; vz00];

    H_aux = H;

    for H_subs = 1:6

        H_aux = subs(H_aux,var_vec(H_subs),i0(H_subs));

    end

    Y_vec_s = solve(H_aux-Href,var_vec(5));

    Solution_amount = length(Y_vec_s);

    jump = 0;

    if k > Solution_amount

        jump = 1
        k
        Solution_amount
        break

    end

    vy00 = vpa(Y_vec_s(k));

if isreal(vy00)

    vy_lim = 0; %%% REVISAR POR EL TEMA DE LA NO EXACTITUD %%%;

    i0 = double([x00; y00; z00; vx00; vy00; vz00]);

    options = odeset('AbsTol',1e-10,'RelTol',1e-10,'Events',@(t,y) cruce(t,y));
    %Integración
    [t00,p00,te,pe,ie] = ode113(@(t,y) CR3BP_PN_1_ODE(t,y,PN_1,mu), [t0 tf], i0, options);

    fprintf('Calculated:    %i  -  %i  -  %i  -  %i  -  %i\n',g,w,k,i,j)
    %Sección de Poincaré del eje x
    
    %pe
    pe_length = size(pe);
    pe_length = pe_length(1);

    if ~isempty(pe)
    % 
             x_cr_v1 = [x_cr_v1; pe(:,1)];
    %         x_cr_v2 = [x_cr_v2; pe(:,1)];
             vx_cr_v1 = [vx_cr_v1; pe(:,4)];
    %         vx_cr_v2 = [vx_cr_v2; pe(:,4)];
    %         vy_cr_v = [vy_cr_v; pe(:,5)];
    % 
    %     for i_2 = 1:length(pe_length)
    % 
    %           if vy_cr_v(i_2)+vy_lim < 0
    %               x_cr_v1(end-pe_length+i_2) = NaN;
    %               vx_cr_v1(end-pe_length+i_2) = NaN;                  
    %           else
    %               x_cr_v2(end-pe_length+i_2) = NaN;
    %               vx_cr_v2(end-pe_length+i_2) = NaN;                  
    %           end
    % 
    % 
    % 
    %           % if abs(pe(i_2,1)) < 3
    %           %     if abs(pe(i_2,4)) < 20
    %           %         plot(pe(i_2,1),(pe(i_2,4)),".","Color",[0 0 k-1],'MarkerSize',5);
    %           %     end
    %           % end
    % 
    %     end
    % end
end

end
    if jump == 1

        break;

    end
end
    

    figure(b+c_aux*(b_values_end-b_values_init+1))

    title('C =', num2str(C))
    hold on

    % for g = 1:length(x_cr_v1)
    % 
    %     if abs(x_cr_v1(g))>2 || abs(x_cr_v2(g))>2
    %         x_cr_v1(g)=NaN;
    %         vx_cr_v1(g)=NaN;
    %         x_cr_v2(g)=NaN;
    %         vx_cr_v2(g)=NaN;
    %     elseif abs(vx_cr_v1(g))>20 || abs(vx_cr_v2(g))>20
    %         x_cr_v1(g)=NaN;
    %         vx_cr_v1(g)=NaN;
    %         x_cr_v2(g)=NaN;
    %         vx_cr_v2(g)=NaN;
    %     end
    %     % if H(x_cr_v(g),0,vx_cr_v(g),vy_cr_v(g))-0.1 > -C/2
    %     %     x_cr_v(g)=NaN;
    %     %     vx_cr_v(g)=NaN;
    %     % end
    % 
    % end
    
    %plot(x_cr_v1,vx_cr_v1,".","Color",[0 0 1],'MarkerSize',4)
    %plot(x_cr_v2,vx_cr_v2,".","Color",[0 0 0],'MarkerSize',4)

    for  x_cr_var = 1:length(x_cr_v1)
        if abs(x_cr_v1(x_cr_var))>2.4 || abs(vx_cr_v1(x_cr_var))>5
            x_cr_v1(x_cr_var) = NaN;
            vx_cr_v1(x_cr_var) = NaN;
        end
    end

    plot(x_cr_v1,vx_cr_v1,".","Color",[0 0 0],'MarkerSize',2)

end
end

    %ylim([-6 6])
    %xlim([-mu-1.5 1-mu+1.5])
end
end
end
end

function [value,isterminal,direction] = cruce(t,y)

  value = y(2); % The value that we want to be zero 
  isterminal = 0;  % Halt integration 
  direction = 1;   % The zero can be approached from either direction

    % if ((y(1)+mu)^2 + y(2)^2 + y(3)^2)^0.5 - 0.005 < 0 || ((y(1)-1+mu)^2 + y(2)^2 + y(3)^2)^0.5 - 0.005 < 0
    % 
    %     value = 0;
    %     isterminal = 1;
    %     direction = 2;
    % 
    % end


end

