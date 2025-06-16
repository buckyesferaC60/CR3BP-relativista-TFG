function pol = Pol_carac(m)

    syms lambda 

    tam = size(m);
    tam = tam(1)

    Ad = eye(tam,tam)*lambda;

    m_2 = m-Ad

    pol = det(m_2);

end

