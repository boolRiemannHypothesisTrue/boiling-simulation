function val = Qu_dQu_duy(params, phase, Q, dy)

    % ∑_i▒∑_j▒〖Q_u^(i1,j1)  (∂Q_u^(i1,j1))/(∂u_y^(i,j) )〗
    % 

    rho = params.rho_liquid * (1 - phase) + params.rho_vapor * phase;
    [M, N] = size(Q);

    val = zeros(M, N);
    
    for i = 2:N-1
        for j = 2:M-1
            val(i,j) = rho(i,j) * (Q(i,j-1) - Q(i,j+1))/ (2*dy);
        end
    end


end
