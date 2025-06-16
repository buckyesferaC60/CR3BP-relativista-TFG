function H= Hamiltonian(Lagrangian,var_vec)

    n_dim = length(var_vec)/3;

    H = -Lagrangian;

    for i = 1:n_dim

        P(i).val = simplify(diff(Lagrangian,var_vec(n_dim+i)));
        H = H + P(i).val*var_vec(n_dim+i);

    end

end

