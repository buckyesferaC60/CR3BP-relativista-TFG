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
[L,~,~] = PN_1_Lagrange_MovEquations(mu,0);

for f = 1:length(kappa_vec)
    [L_PN(f).L,~,~] = PN_1_Lagrange_MovEquations(mu,kappa_vec(f));
end

for l = 1:3

%l = 1;

Orbit_struct.family = strcat('Ly_L',num2str(l));



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
    
    Ly(l).Mndrmy_0 = zeros(L(l).norb_t,36);
    Ly(l).eig_val_0 = zeros(L(l).norb_t,6);
    Ly(l).max_re_0 = zeros(L(l).norb_t,1);
    Ly(l).nu_re_0 = zeros(L(l).norb_t,1);
    Ly(l).max_im_0 = zeros(L(l).norb_t,1);
    Ly(l).eigval_list_0 = zeros(L(l).norb_t,6);
    
    for j = 1:length(kappa_vec)
        Ly(l).K(j).Mndrmy = zeros(L(l).norb_t,36);
        Ly(l).K(j).eig_val = zeros(L(l).norb_t,6);
        Ly(l).K(j).max_re= zeros(L(l).norb_t,1);
        Ly(l).K(j).nu_re = zeros(L(l).norb_t,1);
        Ly(l).K(j).max_im = zeros(L(l).norb_t,1);
        Ly(l).K(j).eigval_list = zeros(L(l).norb_t,6);
    end

    for f = 1:L(l).norb_t
    
        f
    
        Orbit_struct.kappa = 0;
    
        Orbit_struct.ci = L(l).vec(1) + L(l).x0 + ((f-1)*L(l).stepstep + L(l).norb_0 - 1)*L(l).step;
        p = Read_Orbit_final(Orbit_struct);
    
        Ly(l).Full_Orb(f).val = p;
    
        Ly(l).Mndrmy_0(f,:) = p(end,8:end);
    
        [~,eig_val,nu_re,eigval_list] = Monodromy_eigen(Matrixize(p(end,8:end),[6 6]));
    
        Ly(l).eig_val_0(f,:) = eig_val;
        Ly(l).max_re_0(f) = max(abs(real(eig_val)));
        Ly(l).nu_re_0(f) = nu_re;
        Ly(l).max_im_0(f) = max(abs(imag(eig_val)));
        Ly(l).eigval_list_0(f,:) = eigval_list;
    
        for j = 1:length(kappa_vec)
    
            kappa = kappa_vec(j);
    
            Orbit_struct.kappa = kappa;
    
            Orbit_struct.ci = L_PN(j).L(l).vec(1) + L(l).x0 + ((f-1)*L(l).stepstep + L(l).norb_0 - 1)*L(l).step;
    
            p = Read_Orbit_final(Orbit_struct);
    
            Ly(l).K(j).Full_Orb(f).val = p;
    
            Ly(l).K(j).Mndrmy(f,:) = p(end,8:end);
    
            [~,eig_val,nu_re,eigval_list] = Monodromy_eigen(Matrixize(p(end,8:end),[6 6]));
    
            Ly(l).K(j).eig_val(f,:) = eig_val;
            Ly(l).K(j).max_re(f) = max(abs(real(eig_val)));
            Ly(l).K(j).nu_re(f) = nu_re;
            Ly(l).K(j).max_im(f) = max(abs(imag(eig_val)));
            Ly(l).K(j).eigval_list(f,:) = eigval_list;
    
        end
    
    end

% Máximos de parte real e imaginaria (diferencia)

    for k = 1:length(kappa_vec)
    
        figure(2*((l-1)*length(kappa_vec)+k)-1)
        grid on
        hold on
        title(strcat('λ límite    L',num2str(l),'  κ=',num2str(kappa_vec(k)),'  x-xL=',num2str(L(l).step*L(l).stepstep)))
    
        %plot(1:(size(Ly(l).nu_re_0)),Ly(l).nu_re_0,'.','MarkerSize',5,'Color',[0 0 0])
        plot(1:(size(Ly(l).max_re_0)),real(log10(Ly(l).max_re_0-Ly(l).K(k).max_re)),'Color',[0 0 0])%,'.','MarkerSize',5
        plot(1:(size(Ly(l).max_re_0)),imag(log10(Ly(l).max_re_0-Ly(l).K(k).max_re)),'Color',[1 0 0])
        plot(1:(size(Ly(l).max_re_0)),real(log10(Ly(l).max_im_0-Ly(l).K(k).max_im)),'Color',[0 0 1])
        plot(1:(size(Ly(l).max_re_0)),imag(log10(Ly(l).max_im_0-Ly(l).K(k).max_im)),'Color',[1 0 1])
    
        xlabel('x-xL')
        ylabel('λ')

        % legend('Re(log10(max|λ| clásico - max|λ| relativista))',...
        %        'Im(log10(max|λ| clásico - max|λ| relativista))',...
        %        'Re(log10(max(Im(λ)) clásico - max(Im(λ)) relativista))',...
        %        'Im(log10(max(Im(λ)) clásico - max(Im(λ)) relativista))')
    
        figure(2*((l-1)*length(kappa_vec)+k))
        grid on
        hold on
        title(strcat('λ límite    L',num2str(l),'  κ=',num2str(kappa_vec(k)),'  x-xL=',num2str(L(l).step*L(l).stepstep)))
    
        %plot(1:(size(Ly(l).nu_re_0)),Ly(l).nu_re_0,'.','MarkerSize',5,'Color',[0 0 0])
        plot(1:(size(Ly(l).max_re_0)),Ly(l).max_re_0-Ly(l).K(k).max_re,'Color',[0 0 0])%,'.','MarkerSize',5
        plot(1:(size(Ly(l).max_re_0)),Ly(l).max_im_0-Ly(l).K(k).max_im,'Color',[0 0 1])
    
        xlabel('x-xL')
        ylabel('λ')

        % legend('max|λ| clásico - max|λ| relativista',...
        %        'max(Im(λ)) clásico - max(Im(λ)) relativista')
    end

end


