function val = Qb_uy1_dQ_duy(Q)



   
    [M, N] = size(Q);

    val = zeros(M, N);
    
    for i = 1:M
     
            val(i,1) = Q(i,1);
     
    end


end
