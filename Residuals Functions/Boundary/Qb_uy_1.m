function val = Qb_uy_1(fields)
    %
    % краевое условие u_y(y = 0) = 0
    % невязка будет равна u_y на левом краю
    %
    u_y = fields.u_y;
    
    [M,N] = size(u_y);

    val = zeros(M,N);

    i = 1:M;
    
    val(i,1) = u_y(i,1);

end

