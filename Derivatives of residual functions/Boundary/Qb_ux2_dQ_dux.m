function val = Qb_ux2_dQ_dux(Q)



   
    [M, N] = size(Q);

    val = zeros(M, N);
    
    for i = 1:M
     
            val(i,N) = Q(i,N);
     
    end


end
