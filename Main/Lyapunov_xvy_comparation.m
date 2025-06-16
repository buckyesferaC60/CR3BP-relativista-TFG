% Obtain the limit orbits

close all
clear
clc

cd('C:\Matlab\CodigosTFG\PN-1')

%% SYSTEM PARAMETERS

mu = 0.01215;
Orbit_struct.mu = mu;

kappa_vec = [1e-5 1e-7 9.8*1e-9 1.16555*1e-11];

% for i = 1:1
%     figure(i)
%     grid on
%     hold on
% end

% Lagrange points
[L,~,L(1).H,~] = PN_1_Lagrange_MovEquations(mu,0);
L(1).Hfunc = matlabFunction(L(1).H); % x xP y yP z zP

for i = 1:length(kappa_vec)
    [L_PN(i).L,~,L_PN(i).L(1).H,~] = PN_1_Lagrange_MovEquations(mu,kappa_vec(i));
    L_PN(i).L(1).Hfunc = matlabFunction(L_PN(i).L(1).H);
end

for l = 1:3

%l = 1;

Orbit_struct.family = strcat('Ly_L',num2str(l));

n_tot(3) = 401;
n_tot(2) = 370;
n_tot(1) = 440;

step = 1;

%% Classical Lyapunov parameters

kappa = 0;



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

Orbit_struct.kappa = kappa;

% Ly(l).Mndrmy_0 = zeros(n_tot(l),36);
% Ly(l).eig_val_0 = zeros(n_tot(l),6);
% Ly(l).max_re_0 = zeros(n_tot(l),1);
% Ly(l).nu_re_0 = zeros(n_tot(l),1);
% Ly(l).max_im_0 = zeros(n_tot(l),1);
% Ly(l).eigval_list_0 = zeros(n_tot(l),6);

% for j = 1:length(kappa_vec)
%     Ly(l).K(j).Mndrmy = zeros(n_tot(l),36);
%     Ly(l).K(j).eig_val = zeros(n_tot(l),6);
%     Ly(l).K(j).max_re= zeros(n_tot(l),1);
%     Ly(l).K(j).nu_re = zeros(n_tot(l),1);
%     Ly(l).K(j).max_im = zeros(n_tot(l),1);
%     Ly(l).K(j).eigval_list = zeros(n_tot(l),6);
% end

for f = 1:L(l).norb_t

    f

    Orbit_struct.kappa = 0;

    Orbit_struct.ci = L(l).vec(1) + L(l).x0 + ((f-1)*L(l).stepstep + L(l).norb_0 - 1)*L(l).step;
    p = Read_Orbit_final(Orbit_struct);

    Ly(l).Full_Orb(f).val = p;

    for j = 1:length(kappa_vec)

        kappa = kappa_vec(j);

        Orbit_struct.kappa = kappa;

        Orbit_struct.ci = L_PN(j).L(l).vec(1) + L(l).x0 + ((f-1)*L(l).stepstep + L(l).norb_0 - 1)*L(l).step;
    
        p = Read_Orbit_final(Orbit_struct);

        Ly(l).K(j).Full_Orb(f).val = p;

    end

end

% Condiciones iniciales

    for k = 1:length(kappa_vec)
    
        figure(3*((l-1)*length(kappa_vec)+k)-2)
        grid on
        hold on
        title(strcat('Diferencias en E           L',num2str(l),'  κ=',num2str(kappa_vec(k)),'  x-xL=',num2str(L(l).step*L(l).stepstep)))

        for f = 1:L(l).norb_t
            plot(f,real(log10(- L(1).Hfunc(Ly(l).Full_Orb(f).val(1,2),0,0,Ly(l).Full_Orb(f).val(1,6),0,0) + L_PN(k).L(1).Hfunc(Ly(l).K(k).Full_Orb(f).val(1,2),0,0,Ly(l).K(k).Full_Orb(f).val(1,6),0,0))),'o','MarkerSize',3,'Color',[0 0 0])
            plot(f,- L(1).Hfunc(Ly(l).Full_Orb(f).val(1,2),0,0,Ly(l).Full_Orb(f).val(1,6),0,0) + L_PN(k).L(1).Hfunc(Ly(l).K(k).Full_Orb(f).val(1,2),0,0,Ly(l).K(k).Full_Orb(f).val(1,6),0,0),'+','MarkerSize',3,'Color',[0 0 1])
        end

        xlabel('x-xL')
        ylabel('Diferencia de E')

        legend('log10(E relativista - E clásica)','E relativista - E clásica')
    %-----------------------------------------------------------------------------------
        figure(3*((l-1)*length(kappa_vec)+k)-1)
        grid on
        hold on
        title(strcat('Diferencias en vy0           L',num2str(l),'  κ=',num2str(kappa_vec(k)),'  x-xL=',num2str(L(l).step*L(l).stepstep)))

        for f = 1:L(l).norb_t
            plot(f,real(log10(Ly(l).Full_Orb(f).val(1,6)-Ly(l).K(k).Full_Orb(f).val(1,6))),'o','MarkerSize',3,'Color',[0 0 0])
            plot(f,Ly(l).Full_Orb(f).val(1,6)-Ly(l).K(k).Full_Orb(f).val(1,6),'+','MarkerSize',3,'Color',[0 0 1])
        end

        xlabel('x-xL')
        ylabel('Diferencia de vy0')

        legend('log10(vy0 clásico - vy0 relativista)','vy0 clásico - vy0 relativista')
    %------------------------------------------------------------------------------------
        figure(3*((l-1)*length(kappa_vec)+k))
        grid on
        hold on
        title(strcat('Diferencias en T           L',num2str(l),'  κ=',num2str(kappa_vec(k)),'  x-xL=',num2str(L(l).step*L(l).stepstep)))

        for f = 1:L(l).norb_t
            plot(f,real(log10(Ly(l).Full_Orb(f).val(end,1)-Ly(l).K(k).Full_Orb(f).val(end,1))),'o','MarkerSize',3,'Color',[0 0 0])
            plot(f,Ly(l).Full_Orb(f).val(end,1)-Ly(l).K(k).Full_Orb(f).val(end,1),'+','MarkerSize',3,'Color',[0 0 1])
        end

        xlabel('x-xL')
        ylabel('Diferencia de T')

        legend('log10(T relativista - T clásico)','T relativista - T clásico')

    end

end


