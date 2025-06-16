%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @author: Pablo Dalda Díaz
% Alumno de 4º del grado en FIE
% 17/05/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Órbitas peculiares    :   mu   k     x         y   z   vx  vy         vy   vz  
%   Lyapunov pequeña L3 :   0.01 1e-11 L3(1)+0.1 0   0   0   -0.2075    0    0
%   
%   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parámetros del modelo

close all
clear
clc

dt = 0.1;     %paso temporal
t0 = 0;         %instante inicial (s)
tf = 0.4;      %instante final (s)

mu = 0.01528; %0.01215;

%k_vec = 
k = 1e-3;

tol = 0.00000000001;

syms xS yS zS xpS ypS zpS xppS yppS zppS

var_vec = [xS yS zS xpS ypS zpS xppS yppS zppS];

var_vec_0 = var_vec*0;
var_vec_1 = var_vec*0;

[PN_0,H_0,var_vec_0] = Cl_Mov_Eq(mu);
[PN_1,H_1,var_vec_1] = PN_1_Mov_Eq(mu,k);

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

        PN_1(i).val = subs(subs(subs(PN_1(i).val,var_vec(7),PN_0(1).Eq),var_vec(8),PN_0(2).Eq),var_vec(9),PN_0(3).Eq) + PN_0(i).Eq,var_vec(6+i);
        PN_1(i).func = matlabFunction(PN_1(i).val);
 
    end

%%

[L1,L2,L3,L4,L5] = Lpoints(PN_1,mu,var_vec,tol);

L1 = [L1 0 0];
L2 = [L2 0 0];
L3 = [L3 0 0];
L4 = [L4(1) L4(2) 0];
L5 = [L5(1) L5(2) 0];

figure(1)
hold on
plot3(L1(1),L1(2),L1(3),'+',"Color",[k 0 2*mu],'MarkerSize',10);
plot3(L2(1),L2(2),L2(3),'+',"Color",[k 0 2*mu],'MarkerSize',10);
plot3(L3(1),L3(2),L3(3),'+',"Color",[k 0 2*mu],'MarkerSize',10);
plot3(L4(1),L4(2),L4(3),'+',"Color",[k 0 2*mu],'MarkerSize',10);
plot3(L5(1),L5(2),L5(3),'+',"Color",[k 0 2*mu],'MarkerSize',10);
plot3(-mu,0,0,'.','Color','k','MarkerSize',10);
plot3(1-mu,0,0,'.','Color','k','MarkerSize',10);

% Y_vec_s = solve(H_sym,ypS);
% 
% Y_vec(1).y = Y_vec_s(1);
% Y_vec(2).y = Y_vec_s(2);
% 
% Y_vec(1).y = matlabFunction(Y_vec(1).y);  %Hs Vx X Y
% Y_vec(2).y = matlabFunction(Y_vec(2).y);  

%C = 3;
%Href =  -C/2;  %H(0.8234,0,0,0.1263); %   
%C = -2*Href


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long

%for w=1:2

x00 = 0.852559195979900;
y00 = 0;
z00 = 0;

vx00 = 0;
vy00 = 0.14738157953651346067301360485155;
vz00 = 0;

H = subs(subs(subs(subs(subs(subs(H,var_vec(1),x00),var_vec(2),y00),var_vec(3),z00),var_vec(4),vx00),var_vec(5),vy00),var_vec(6),vz00);

    i0 = [x00; y00; z00; vx00; vy00; vz00];

    options = odeset('AbsTol',1e-20,'RelTol',1e-10);
    %Integración
    [t00,p00,~] = ode113(@(t,y) CR3BP_PN_1_ODE(t,y,PN_1,mu), [t0 tf], i0);%,'Events', @(t,y) cruce(t,y,mu));

fprintf('\n\n______________________\n H = %d \n_______________________\n\n',double(H));

    figure(1)
    hold on

    plot3(p00(:,1),p00(:,2),p00(:,3),"Color",[0 0 1]);
    %axis equal
    grid on

    xlim([-2 2])
    ylim([-2 2])

    %figure(2)
    hold on

    %plot3(p00(:,4),p00(:,5),p00(:,6),"Color",[0 0 1]);
    %axis equal
    %grid on

    %xlim([-10 10])
    %ylim([-10 10])

    grid on

    axis equal

    % figure(3)
    % plot(H(p00(:,1),p00(:,2),p00(:,4),p00(:,5)),"Color",[0 0 1])
    % 
    % figure(4)
    % plot(real(H(p00(:,1),p00(:,2),p00(:,4),p00(:,5))),"Color",[0 0 1])
    % 
    % grid on

%end
function [value,isterminal,direction] = cruce(t,y,mu)

  value = y(2); % The value that we want to be zero 
  isterminal = 0;  % Halt integration 
  direction = 1;   % The zero can be approached from either direction

    if ((y(1)+mu)^2 + y(2)^2 + y(3)^2)^0.5 - 0.005 < 0 || ((y(1)-1+mu)^2 + y(2)^2 + y(3)^2)^0.5 - 0.005 < 0

        value = 0;
        isterminal = 1;
        direction = 2;

    end


end
%end    