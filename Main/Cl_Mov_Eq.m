function [Eq,H,var_vec] = Cl_Mov_Eq(mu)

    m1 = 1-mu;
    m2 = mu;
    x1 = -mu;
    x2 = 1-mu;

    syms x y z vx vy vz ax ay az

    var_vec = [x y z vx vy vz ax ay az];
    
    T = vx^2 + vy^2 + vz^2 + x^2 + y^2 + 2*(vy*x - vx*y);

    r1 = ((x-x1)^2 + y^2 + z^2)^0.5;
    r2 = ((x-x2)^2 + y^2 + z^2)^0.5;

    L = T/2 + m1/r1 + m2/r2;

    for i = 1:3

        Eq(i).val = simplify(Euler_Lagrange(L,var_vec,i));

    end

    H = Hamiltonian(L,var_vec);

end

