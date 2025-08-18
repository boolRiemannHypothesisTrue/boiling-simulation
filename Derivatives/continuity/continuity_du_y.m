function Sy = continuity_du_y(params, phase, Q, dy)
    % Центральная разностная схема по (21) и (22) для внутренней области
    % Сетка равномерная, края обнуляются

    rho = params.rho_liquid * (1 - phase) + params.rho_vapor * phase;
    [M, N] = size(Q);

    Sy = zeros(M, N);
    
    for i = 2:N-1
        for j = 2:M-1
            Sy(i,j) = rho(i,j) * (Q(i,j-1) - Q(i,j+1))/ (2*dy);
        end
    end


end
