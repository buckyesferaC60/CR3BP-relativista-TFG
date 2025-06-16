function [Eq,H,var_vec] = PN_1_Mov_Eq(mu,kappa)

    m1 = 1-mu;
    m2 = mu;
    x1 = -mu;
    x2 = 1-mu;

    syms x y z vx vy vz ax ay az w G M1 M2 r12 c

    var_vec = [x y z vx vy vz ax ay az];
    
    T = (vx^2+vy^2+vz^2)+2*w*(vy*x-vx*y)+w^2*(x^2+y^2);

    r1 = ((x-x1)^2 + y^2 + z^2)^0.5;
    r2 = ((x-x2)^2 + y^2 + z^2)^0.5;

    w1 = (m1*m2-3)/2;

    L = 1/c^2 *( (T^2/8) + w1*w*(x^2+y^2) + w1*(vy*x-vx*y)...
  - G^2*(M1*M2/r12 * (1/r1+1/r2) + (M1^2/2/r1^2+M2^2/2/r2^2+M1*M2/r1/r2) )...
  + G/2 * (3*(M1/r1+M2/r2)*T + 3*w^2*(M1*x1^2/r1+M2*x2^2/r2) ...
  - 7*w*(vy+w*x)*(M1*x1/r1+M2*x2/r2)...
  - w*y*(vx-w*y)*(M1*x1*(x-x1)/r1^3+M2*x2*(x-x2)/r2^3)...
  - w*y^2*(vy+w*x)*(M1*x1/r1^3+M2*x2/r2^3)...
  - w*y*z*vz*(M1*x1/r1^3+M2*x2/r2^3) ) );

    L = simplify(subs(subs(subs(subs(subs(subs(simplify(L),G,1),M1,m1),M2,m2),c,1/kappa^0.5),w,1),r12,1));

    for i = 1:3

        Eq(i).val = simplify(Euler_Lagrange(L,var_vec,i));

    end

    H = Hamiltonian(L,var_vec);

end

