function Sx = continuity_du_x(params, phase, Q, dx)
    % Центральная разностная схема по (21) и (22) для внутренней области
    % Сетка равномерная, края обнуляются

    rho = params.rho_liquid * (1 - phase) + params.rho_vapor * phase;
    [M, N] = size(Q);

    Sx = zeros(M, N);
    
    for i = 2:N-1
        for j = 2:M-1
            Sx(i,j) = rho(i,j) * (Q(i-1,j) - Q(i+1,j))/ (2*dx);
        end
    end


end
