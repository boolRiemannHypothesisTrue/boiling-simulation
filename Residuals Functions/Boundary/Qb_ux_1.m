function val = Qb_ux_1(fields)
    %
    % краевое условие u_x(y = 0) = 0
    % невязка будет равна u_x на левом краю
    %
    u_x = fields.u_x;
    
    [M,N] = size(u_x);

    val = zeros(M,N);

    i = 1:M;
    
    val(i,1) = u_x(i,1);

end

