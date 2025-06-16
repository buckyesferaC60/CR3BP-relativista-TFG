function out = Traza(Matriz)


    n = size(Matriz);

    n = n(1);

    out = 0;

    for i = 1:n
        out = out + Matriz(i,i);
    end


end