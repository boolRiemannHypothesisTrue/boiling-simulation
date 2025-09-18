function val = Qb_T_2(fields,params)
    %
    % краевое условие T(y = 0) = Tw
    % невязка будет равна T(i,N) - Tl на правом краю
    %
    T = fields.T;
    
    Tl = params.T_liquid;

    [M,N] = size(T);

    val = zeros(M,N);

    i = 1:M;
    
    val(i,N) = T(i,N) - Tl;

end

