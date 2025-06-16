function Eq = Euler_Lagrange(Lagrangian,var_vec,var_index)
      
        Eq = diff(Lagrangian,var_vec(var_index)) - d_dt(diff(Lagrangian,var_vec(var_index + length(var_vec)/3)),length(var_vec)/3,var_vec);
  
end

function out = d_dt(f,n_dim,var_vec)
    
    out = 0;
    for i = 1:2*n_dim
        out = out + diff(f,var_vec(i))*var_vec(i+n_dim);
    end

end