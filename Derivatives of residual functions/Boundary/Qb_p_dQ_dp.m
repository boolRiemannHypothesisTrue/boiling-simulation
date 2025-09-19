function val = Qb_p_dQ_dp(Q)



   
    [M, N] = size(Q);

    val = zeros(M, N);
    
    for i = 1:M
     
            val(i,N-1) = Q(i,N-1);
            
            val(i,N) = - Q(i,N);
    end


end
