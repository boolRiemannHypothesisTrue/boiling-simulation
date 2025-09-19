function val = Qb_uy2_dQ_duy(Q)



   
    [M, N] = size(Q);

    val = zeros(M, N);
    
    for i = 1:M
     
            val(i,N) = Q(i,N);
     
    end


end
