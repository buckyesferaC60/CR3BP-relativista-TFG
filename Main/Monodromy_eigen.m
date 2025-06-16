function [eig_vec,eig_val,nu_re,eigval_list] = Monodromy_eigen(M)
%MONODROMY_EIGEN Summary of this function goes here
%   Detailed explanation goes here

    % Autovalores y autovectores

    [eig_vec,eig_val] = eig(M);

    eig_val = diag(eig_val);
    eigval_list = 0*eig_val;



    % Autovalores con la parte imaginaria mayor
    
    Im_eigval = abs(imag(eig_val));

    Aux = find(Im_eigval==max(Im_eigval));

    for i = 1:length(Aux)
        eigval_list(Aux(i)) = 1;
    end

    % Autovalores con módulo notablemente distinto de 1 (Ws y Wu)

    for  i = 1:length(eig_val)
        if abs(eig_val(i)) == max(abs(eig_val))
            eigval_list(i) = 2;
        end
        if abs(eig_val(i)) == min(abs(eig_val))
            eigval_list(i) = 2;
        end
    end

    % Coeficiente de estabilidad

    nu_re = 0;
    for i = 1:length(eig_val)
        if eigval_list(i) == 2
            nu_re = nu_re + 0.5*abs(eig_val(i));
        end
    end

end

