function val = Qb_uy_2(fields)
    %
    % краевое условие u_y(y = inf) = 0
    % невязка будет равна u_y на правом краю
    %
    u_y = fields.u_y;
    
    [M,N] = size(u_y);

    val = zeros(M,N);

    i = 1:M;
    
    val(i,N) = u_y(i,N);

end

