function val = Qb_ux1_dQ_dux(Q)



   
    [M, N] = size(Q);

    val = zeros(M, N);
    
    for i = 1:M
     
            val(i,1) = Q(i,1);
     
    end


end
