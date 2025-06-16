function Out = Plot_manifolds(Orb,M,PN_1,mu,kappa,L)
%PLOT_MANIFOLDS Summary of this function goes here
%   Detailed explanation goes here

    Out = 0;

    Orb = Orb(:,1:6);

    for i = 1:5
        figure(99+i)
        hold on
        grid on

        plot3(Orb(:,1),Orb(:,2),Orb(:,3),'Color',[0 0 0],'LineWidth',1.5)
        plot3(-mu,0,0,'.','MarkerSize',20,'Color',[0 1 1],'MarkerEdgeColor',[0 0 0],'DisplayName','M1')
        plot3(1-mu,0,0,'.','MarkerSize',20,'Color',[0 1 1],'MarkerEdgeColor',[0 0 0],'DisplayName','M2')
        plot3(L(1).vec(1),0,0,'+','MarkerSize',10,'Color',[0 0 0],'DisplayName','L1')
        plot3(L(2).vec(1),0,0,'+','MarkerSize',10,'Color',[0 0 0],'DisplayName','L2')
        plot3(L(3).vec(1),0,0,'+','MarkerSize',10,'Color',[0 0 0],'DisplayName','L3')
        plot3(L(4).vec(1),L(4).vec(2),0,'+','MarkerSize',10,'Color',[0 0 0],'DisplayName','L4')
        plot3(L(5).vec(1),L(5).vec(2),0,'+','MarkerSize',10,'Color',[0 0 0],'DisplayName','L5')        
    end
           
    % Autovalores y autovectores

    format long
    [eig_vec,eig_val,~,eigval_list] = Monodromy_eigen(M);
    format short

    Aux = find(eigval_list == 2);

    if abs(eig_val(Aux(1))) > abs(eig_val(Aux(2)))
        WU_vec = eig_vec(:,Aux(1));
        WS_vec = eig_vec(:,Aux(2));
    else
        WU_vec = eig_vec(:,Aux(2));
        WS_vec = eig_vec(:,Aux(1));
    end

    % Generación de condiciones iniciales
    n_sim = 50;%round(length(Orb(:,1)));
    t_sim = 9;
    Delta = 0.0001;

    step = round(length(Orb(:,1))/n_sim);

    for i = 1:n_sim
        sigma_0(i).WUM = round(Orb(mod(i*step,length(Orb(:,1)))+1,:) + Delta*WU_vec.',15);
        sigma_0(i).WUm = round(Orb(mod(i*step,length(Orb(:,1)))+1,:) - Delta*WU_vec.',15);
        sigma_0(i).WSM = round(Orb(mod(i*step,length(Orb(:,1)))+1,:) + Delta*WS_vec.',15);
        sigma_0(i).WSm = round(Orb(mod(i*step,length(Orb(:,1)))+1,:) - Delta*WS_vec.',15);       
    end

    

    options = odeset('AbsTol',1e-20,'RelTol',1e-10);

    for i = 1:n_sim

    [tWUM,pWUM] = ode113(@(t,y) CR3BP_PN_1_ODE(t,y,PN_1,mu,kappa), [0 t_sim], sigma_0(i).WUM, options);
    [tWUm,pWUm] = ode113(@(t,y) CR3BP_PN_1_ODE(t,y,PN_1,mu,kappa), [0 t_sim], sigma_0(i).WUm, options);
    [tWSM,pWSM] = ode113(@(t,y) CR3BP_PN_1_ODE(t,y,PN_1,mu,kappa), [0 -t_sim], sigma_0(i).WSM, options);
    [tWSm,pWSm] = ode113(@(t,y) CR3BP_PN_1_ODE(t,y,PN_1,mu,kappa), [0 -t_sim], sigma_0(i).WSm, options);

    figure(100)
    plot3(pWUM(:,1),pWUM(:,2),pWUM(:,3),'Color',[1 0 0],'LineWidth',0.1)
    plot3(pWUm(:,1),pWUm(:,2),pWUm(:,3),'Color',[1 0 0],'LineWidth',0.1)
    plot3(pWSM(:,1),pWSM(:,2),pWSM(:,3),'Color',[0 0 1],'LineWidth',0.1)
    plot3(pWSm(:,1),pWSm(:,2),pWSm(:,3),'Color',[0 0 1],'LineWidth',0.1)

    figure(101)
    plot3(pWUM(:,1),pWUM(:,2),pWUM(:,3),'Color',[1 0 0],'LineWidth',0.1)    
    %plot3(pWUm(:,1),pWUm(:,2),pWUm(:,3),'Color',[1 0 0],'LineWidth',0.1)

    figure(102)
    %plot3(pWSM(:,1),pWSM(:,2),pWSM(:,3),'Color',[0 0 1],'LineWidth',0.1)
    plot3(pWSm(:,1),pWSm(:,2),pWSm(:,3),'Color',[0 0 1],'LineWidth',0.1)

    figure(103)
    %plot3(pWUM(:,1),pWUM(:,2),pWUM(:,3),'Color',[1 0 0],'LineWidth',0.1)
    plot3(pWSM(:,1),pWSM(:,2),pWSM(:,3),'Color',[0 0 1],'LineWidth',0.1)

    figure(104)
    plot3(pWUm(:,1),pWUm(:,2),pWUm(:,3),'Color',[1 0 0],'LineWidth',0.1)
    %plot3(pWSm(:,1),pWSm(:,2),pWSm(:,3),'Color',[0 0 1],'LineWidth',0.1)

    end

    %axis equal
end

