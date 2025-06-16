close all
clear
clc

cd('C:\Matlab\CodigosTFG\PN-1')

%% PARÁMETROS

mu = 0.01215;
Orbit_struct.mu = mu;

kappa = 0;
Orbit_struct.kappa = kappa;

    [L,PN_1,~,var_vec] = PN_1_Lagrange_MovEquations(mu,kappa);

    L(1).x0 = -0.009;
    L(2).x0 = 0.005;
    L(3).x0 = -0.09;

    L(1).step = -0.0004;
    L(2).step = 0.0004;
    L(3).step = -0.002;

    L(1).norb_0 = 6;
    L(2).norb_0 = 10;
    L(3).norb_0 = 6;

    L(1).norb_t = 21;
    L(2).norb_t = 18;
    L(3).norb_t = 21;

    L(1).stepstep = 20;
    L(2).stepstep = 20;
    L(3).stepstep = 20;

    figure(1)
    hold on
    grid on


    color_vec = @(l) [Delta(l-1) Delta(l-2)*0.5 Delta(l)];

%l = 3;

%% PROCESO DE DIBUJADO

n_start = 1;

for f = n_start:21

    for  l = 1:3

        if f>L(l).norb_t
            break;
        end

    Orbit_struct.family = strcat('Ly_L',num2str(l));

    Orbit_struct.ci = L(l).vec(1) + L(l).x0 + ((f-1)*L(l).stepstep + L(l).norb_0 - 1)*L(l).step;
    
    p = Read_Orbit_final(Orbit_struct);
    plot3(p(:,2),p(:,3),p(:,4),'Color',color_vec(l))%,'DisplayName',strcat('Familia de L',num2str(l)))
    
    end
end

plot3(1-mu,0,0,'*','Color',[0 0 0],'MarkerSize',7)
plot3(-mu,0,0,'*','Color',[0 0 0],'MarkerSize',7)
for i = 1:5
   plot3(L(i).vec(1),L(i).vec(2),0,'+','Color',[0 0 0],'MarkerSize',7)
end
legend('Familia de L1','Familia de L2','Familia de L3','')
title('Familias de órbitas de Lyapunov')
xlabel('x')
ylabel('y')

axis equal

%% FUNCIONES 

function Out = Delta(x)
    if x == 1
        Out = 1;
    else
        Out = 0;
    end
end
