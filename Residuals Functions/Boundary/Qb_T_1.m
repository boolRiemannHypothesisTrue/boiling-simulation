function val = Qb_T_1(fields,params)
    %
    % краевое условие T(y = 0) = Tw
    % невязка будет равна T(i,1) - Tw на левом краю
    %
    T = fields.T;
    
    Tw = params.T_wall;

    [M,N] = size(T);

    val = zeros(M,N);

    i = 1:M;
    
    val(i,1) = T(i,1) - Tw;

end

