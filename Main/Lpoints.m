function [L1,L2,L3,L4,L5] = Lpoints(Eq,params,var_vec,tol)

    format long
    
    mu = params(1);

    n_max = -3*log(tol);

    L1 = 0.75-mu;
    L2 = 1.01-mu;
    L3 = -0.75-mu;

    L4 = [0.5-mu; 3^0.5/2];
    L5 = [0.5-mu; -3^0.5/2];

    syms xvar yvar

    F_col = subs(Eq(1).val,var_vec(1),xvar);

    for i = 2:length(var_vec)
        F_col = subs(F_col,var_vec(i),0); 
    end

    F_col_dif = diff(F_col,xvar);

    F_col = matlabFunction(F_col);
    F_col_dif = matlabFunction(F_col_dif);

    F_tri_x = subs(subs(Eq(1).val,var_vec(1),xvar),var_vec(2),yvar);
    F_tri_y = subs(subs(Eq(2).val,var_vec(1),xvar),var_vec(2),yvar);

    for i = 3:length(var_vec)
        F_tri_x = subs(F_tri_x,var_vec(i),0); 
        F_tri_y = subs(F_tri_y,var_vec(i),0); 
    end


    F_tri_x_difx = diff(F_tri_x,xvar);
    F_tri_x_dify = diff(F_tri_x,yvar);
    F_tri_y_difx = diff(F_tri_y,xvar);
    F_tri_y_dify = diff(F_tri_y,yvar);

    F_tri_x = matlabFunction(F_tri_x);
    F_tri_y = matlabFunction(F_tri_y);
    F_tri_x_difx = matlabFunction(F_tri_x_difx);
    F_tri_x_dify = matlabFunction(F_tri_x_dify);
    F_tri_y_difx = matlabFunction(F_tri_y_difx);
    F_tri_y_dify = matlabFunction(F_tri_y_dify);

%------------------------------------------------------------------------------------------

    disp('\n L1 ----\n')
    count = 0;
    while tol < abs(F_col(L1)/F_col_dif(L1))

        L1 = L1 - F_col(L1)/F_col_dif(L1);
        count = count + 1;
        fprintf('_______Iteration_:_%d___L1_x_value_:_%.15f\n',count,L1)
    
    end
    % L1 = [L1;0;0];

    disp('\n- L2 ---\n')
    count = 0;
    while tol < abs(F_col(L2)/F_col_dif(L2))

        L2 = L2 - F_col(L2)/F_col_dif(L2); 
        count = count + 1;
        fprintf('_______Iteration_:_%d___L2_x_value_:_%.15f\n',count,L2)

    end
    % L2 = [L2;0;0];

    disp('\n-- L3 --\n')
    count = 0;
    while tol < abs(F_col(L3)/F_col_dif(L3))

        L3 = L3 - F_col(L3)/F_col_dif(L3);        
        count = count + 1;
        fprintf('_______Iteration_:_%d___L3_x_value_:_%.15f\n',count,L3)

    end   
    % L3 = [L3;0;0];

%---------------------------------------------------------------------------------------------------

    Aux_det = @(x,y) F_tri_x_difx(x,y)*F_tri_y_dify(x,y) - F_tri_x_dify(x,y)*F_tri_y_difx(x,y);
    
    disp('\n--- L4 -\n')
    count = 0;
    L4_alt = 10*tol;
    while tol < norm(L4_alt-L4) 
        
        L4_alt = L4;
        %L4 = L4_alt - [(F_tri_x(L4_alt(1),L4_alt(2))*F_tri_y_dify(L4_alt(1),L4_alt(2)) - F_tri_y(L4_alt(1),L4_alt(2))*F_tri_x_dify(L4_alt(1),L4_alt(2))) / Aux_det(L4_alt(1),L4_alt(2)); (F_tri_y(L4_alt(1),L4_alt(2))*F_tri_x_difx(L4_alt(1),L4_alt(2)) - F_tri_x(L4_alt(1),L4_alt(2))*F_tri_y_dify(L4_alt(1),L4_alt(2))) / Aux_det(L4_alt(1),L4_alt(2))] ;
        L4 = L4 - [F_tri_x(L4(1),L4(2)) / ( F_tri_x_difx(L4(1),L4(2)) - F_tri_x_dify(L4(1),L4(2)) / F_tri_y_dify(L4(1),L4(2)) * F_tri_y_difx(L4(1),L4(2)) ) + F_tri_y(L4(1),L4(2)) / ( F_tri_y_difx(L4(1),L4(2)) - F_tri_y_dify(L4(1),L4(2)) / F_tri_x_dify(L4(1),L4(2)) * F_tri_x_difx(L4(1),L4(2)) );...
                   F_tri_y(L4(1),L4(2)) / ( F_tri_y_dify(L4(1),L4(2)) - F_tri_y_difx(L4(1),L4(2)) / F_tri_x_difx(L4(1),L4(2)) * F_tri_x_dify(L4(1),L4(2)) ) + F_tri_x(L4(1),L4(2)) / ( F_tri_x_dify(L4(1),L4(2)) - F_tri_x_difx(L4(1),L4(2)) / F_tri_y_difx(L4(1),L4(2)) * F_tri_y_dify(L4(1),L4(2)) )];

        count = count+1;
        fprintf('_______Iteration_:_%d___L4_x_value_:_%.15f___L4_y_value_:_%.15f\n',count,L4(1),L4(2))
        if count>n_max
            L4 = [NaN,NaN];
            break;
        end
    end
    if abs(L4(1)) > 2
        L4(1) = NaN;
        L4(2) = NaN;
    end
    if abs(L4(2)) > 2
        L4(1) = NaN;
        L4(2) = NaN;
    end

    disp('\n---- L5 \n')
    count = 0;
    L5_alt = 10*tol;
    while tol < norm(L5_alt-L5) 

        L5_alt = L5;
        %L5 = L5_alt - [(F_tri_x(L5_alt(1),L5_alt(2))*F_tri_y_dify(L5_alt(1),L5_alt(2)) - F_tri_y(L5_alt(1),L5_alt(2))*F_tri_x_dify(L5_alt(1),L5_alt(2))) / Aux_det(L5_alt(1),L5_alt(2)); (F_tri_y(L5_alt(1),L5_alt(2))*F_tri_x_difx(L5_alt(1),L5_alt(2)) - F_tri_x(L5_alt(1),L5_alt(2))*F_tri_y_dify(L5_alt(1),L5_alt(2))) / Aux_det(L5_alt(1),L5_alt(2))];
        L5 = L5 - [F_tri_x(L5(1),L5(2)) / ( F_tri_x_difx(L5(1),L5(2)) - F_tri_x_dify(L5(1),L5(2)) / F_tri_y_dify(L5(1),L5(2)) * F_tri_y_difx(L5(1),L5(2)) ) + F_tri_y(L5(1),L5(2)) / ( F_tri_y_difx(L5(1),L5(2)) - F_tri_y_dify(L5(1),L5(2)) / F_tri_x_dify(L5(1),L5(2)) * F_tri_x_difx(L5(1),L5(2)) );
            F_tri_y(L5(1),L5(2)) / ( F_tri_y_dify(L5(1),L5(2)) - F_tri_y_difx(L5(1),L5(2)) / F_tri_x_difx(L5(1),L5(2)) * F_tri_x_dify(L5(1),L5(2)) ) + F_tri_x(L5(1),L5(2)) / ( F_tri_x_dify(L5(1),L5(2)) - F_tri_x_difx(L5(1),L5(2)) / F_tri_y_difx(L5(1),L5(2)) * F_tri_y_dify(L5(1),L5(2)) )];

        count = count+1;
        fprintf('_______Iteration_:_%d___L5_x_value_:_%.15f___L5_y_value_:_%.15f\n',count,L5(1),L5(2))
        if count>n_max
            L4 = [NaN,NaN];
            break;
        end
    end

    %disp('\n')

    if abs(L5(1)) > 3
        L5(1) = NaN;
        L5(2) = NaN;
    end
    if abs(L5(2)) > 3
        L5(1) = NaN;
        L5(2) = NaN;
    end

    format short

end

