function val = Qb_p(fields)
    %
    % краевое условие dp/dy ( y = inf ) = 0;
    % невязка будет равна p(i,N) - p(i,N-1) на правом(и предправом) краю
    %
    p = fields.P;

    [M,N] = size(p);

    val = zeros(M,N);

    i = 1:M;
    
    val(i,N) = p(i,N) - p(i,N-1);

end

