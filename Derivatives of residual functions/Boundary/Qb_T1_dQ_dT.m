function val = Qb_T1_dQ_dT(Q)



   
    [M, N] = size(Q);

    val = zeros(M, N);
    
    for i = 1:M
     
            val(i,1) = Q(i,1);
     
    end


end
