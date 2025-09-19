function val = Qb_T2_dQ_dT(Q)



   
    [M, N] = size(Q);

    val = zeros(M, N);
    
    for i = 1:M
     
            val(i,N) = Q(i,N);
     
    end


end
