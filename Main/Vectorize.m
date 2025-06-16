function Vec = Vectorize(M)
%VECTORIZE Transform a axb matrix into a ab column vector
%  Horizontal preference
    n = size(M);

    Vec = zeros(n(1)*n(2),1);

    for i = 1:n(1)
        Vec((i-1)*n(2)+1:i*n(2)) = M(i,:); 
    end

end

