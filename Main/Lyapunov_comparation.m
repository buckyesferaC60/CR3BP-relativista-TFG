% Obtain the limit orbits

close all
clear
clc


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
[L,~,~] = PN_1_Lagrange_MovEquations(mu,0);

for i = 1:length(kappa_vec)
    [L_PN(i).L,~,~] = PN_1_Lagrange_MovEquations(mu,kappa_vec(i));
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



%Orbit_struct.kappa = kappa;

Ly(l).Mndrmy_0 = zeros(n_tot(l),36);
Ly(l).eig_val_0 = zeros(n_tot(l),6);
Ly(l).max_re_0 = zeros(n_tot(l),1);
Ly(l).nu_re_0 = zeros(n_tot(l),1);
Ly(l).max_im_0 = zeros(n_tot(l),1);
Ly(l).eigval_list_0 = zeros(n_tot(l),6);

for j = 1:length(kappa_vec)
    Ly(l).K(j).Mndrmy = zeros(n_tot(l),36);
    Ly(l).K(j).eig_val = zeros(n_tot(l),6);
    Ly(l).K(j).max_re= zeros(n_tot(l),1);
    Ly(l).K(j).nu_re = zeros(n_tot(l),1);
    Ly(l).K(j).max_im = zeros(n_tot(l),1);
    Ly(l).K(j).eigval_list = zeros(n_tot(l),6);
end

for i = 1:n_tot(l)/step

    i

    Orbit_struct.kappa = 0;

    Orbit_struct.ci = L(l).vec(1) + L(l).x0 + step*(i-1)*L(l).step;

    p = Read_Orbit(Orbit_struct);

    Ly(l).Full_Orb(i).val = p;

    Ly(l).Mndrmy_0(i,:) = p(end,7:end);

    [~,eig_val,nu_re,eigval_list] = Monodromy_eigen(Matrixize(p(end,7:end),[6 6]));

    Ly(l).eig_val_0(i,:) = eig_val;
    Ly(l).max_re_0(i) = max(abs(real(eig_val)));
    Ly(l).nu_re_0(i) = nu_re;
    Ly(l).max_im_0(i) = max(abs(imag(eig_val)));
    Ly(l).eigval_list_0(i,:) = eigval_list;

    for j = 1:length(kappa_vec)

        kappa = kappa_vec(j);

        Orbit_struct.kappa = kappa;

        Orbit_struct.ci = L_PN(j).L(l).vec(1) + L(l).x0 + step*(i-1)*L(l).step;

        p = Read_Orbit(Orbit_struct);

        Ly(l).K(j).Full_Orb(i).val = p;

        Ly(l).K(j).Mndrmy(i,:) = p(end,7:end);

        [~,eig_val,nu_re,eigval_list] = Monodromy_eigen(Matrixize(p(end,7:end),[6 6]));

        Ly(l).K(j).eig_val(i,:) = eig_val;
        Ly(l).K(j).max_re(i) = max(abs(real(eig_val)));
        Ly(l).K(j).nu_re(i) = nu_re;
        Ly(l).K(j).max_im(i) = max(abs(imag(eig_val)));
        Ly(l).K(j).eigval_list(i,:) = eigval_list;

    end

end

% Coeficientes de estabilidad

for k = 1:4

figure(16*l-3-(4*(k-1)))
grid on
hold on
%plot(1:(size(Ly(l).nu_re_0)),Ly(l).nu_re_0,'.','MarkerSize',5,'Color',[0 0 0])
plot(1:(size(Ly(l).max_re_0)),real(log10(Ly(l).max_re_0-Ly(l).K(k).max_re)),'.','MarkerSize',5,'Color',[0 0 0])
plot(1:(size(Ly(l).max_re_0)),real(log10(Ly(l).max_im_0-Ly(l).K(k).max_im)),'.','MarkerSize',5,'Color',[0 0 1])


figure(16*l-2-(4*(k-1)))
grid on
hold on
plot(1:(size(Ly(l).nu_re_0)),Ly(l).K(k).max_re,'.','MarkerSize',5,'Color',[0 0 0])
plot(1:(size(Ly(l).nu_re_0)),Ly(l).max_re_0,'.','MarkerSize',5,'Color',[1 0 0])
plot(1:(size(Ly(l).max_im_0)),Ly(l).K(k).max_im,'.','MarkerSize',5,'Color',[0 0 1])
plot(1:(size(Ly(l).max_im_0)),Ly(l).max_im_0,'.','MarkerSize',5,'Color',[1 0 1])


%plot(1:(size(Ly(l).nu_re_0)),real(log10(Ly(l).nu_re_0-Ly(l).K(2).nu_re)),'.','MarkerSize',5,'Color',[0 0 0.5])
%plot(1:(size(Ly(l).nu_re_0)),real(log10(Ly(l).nu_re_0-Ly(l).K(3).nu_re)),'.','MarkerSize',5,'Color',[0 0 0.75])
%plot(1:(size(Ly(l).nu_re_0)),real(log10(Ly(l).nu_re_0-Ly(l).K(4).nu_re)),'.','MarkerSize',5,'Color',[0 0 1])
% plot(1:(size(Ly(l).nu_re_0)),imag(log10(Ly(l).nu_re_0-Ly(l).K(1).nu_re)),'.','MarkerSize',5,'Color',[0.25 0 0])
% plot(1:(size(Ly(l).nu_re_0)),imag(log10(Ly(l).nu_re_0-Ly(l).K(2).nu_re)),'.','MarkerSize',5,'Color',[0.5 0 0])
% plot(1:(size(Ly(l).nu_re_0)),imag(log10(Ly(l).nu_re_0-Ly(l).K(3).nu_re)),'.','MarkerSize',5,'Color',[0.75 0 0])
% plot(1:(size(Ly(l).nu_re_0)),imag(log10(Ly(l).nu_re_0-Ly(l).K(4).nu_re)),'.','MarkerSize',5,'Color',[1 0 0])

% Condiciones iniciales

% figure(16*l-1-(4*(k-1)))
% grid on
% hold on
%     for i = 1:n_tot(l)/step
%         plot(i,real(log10(Ly(l).Full_Orb(i).val(1)-Ly(l).K(k).Full_Orb(i).val(1))),'.','MarkerSize',5,'Color',[0 0 0])
%     end
% figure(16*l-(4*(k-1)))
% grid on
% hold on
%     for i = 1:n_tot(l)/step
%         plot(i,real(log10(Ly(l).Full_Orb(i).val(5)-Ly(l).K(k).Full_Orb(i).val(5))),'.','MarkerSize',5,'Color',[0 0 0])
%     end
% end

end

end
