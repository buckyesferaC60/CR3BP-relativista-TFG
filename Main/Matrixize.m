function M = Matrixize(Vec,M_size)
%MATRIXIZE Inverse of Vectorize
    M = zeros(M_size(1),M_size(2));

    for i = 1:M_size(1)
        M(i,:) = Vec((i-1)*M_size(2)+1:i*M_size(2)) ; 
    end
end

