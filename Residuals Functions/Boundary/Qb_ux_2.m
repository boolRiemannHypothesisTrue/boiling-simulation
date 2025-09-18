function val = Qb_ux_2(fields)
    %
    % краевое условие u_x(y = inf) = 0
    % невязка будет равна u_x на правом краю
    %
    u_x = fields.u_x;
    
    [M,N] = size(u_x);

    val = zeros(M,N);

    i = 1:M;
    
    val(i,N) = u_x(i,N);

end

