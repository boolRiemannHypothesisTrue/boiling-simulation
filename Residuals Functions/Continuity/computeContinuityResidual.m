function val = computeContinuityResidual(fields, params, phase, dx, dy)
    % Вычисление невязки уравнения непрерывности div(rho * u) = 0
    % Только во внутренней области, края = 0

    u_x = fields.u_x;
    u_y = fields.u_y;


    rho = params.rho_liquid * (1 - phase) + params.rho_vapor * phase;


    rhou_x = rho .* u_x;
    rhou_y = rho .* u_y;


    [M, N] = size(rho);

    val = zeros(M, N);

    i = 2:M-1;
    j=2:N-1;
    
    val(i,j) = (rhou_x(i+1,j) - rhou_x(i-1,j)) / (2 * dx) + ...
                (rhou_y(i,j+1) - rhou_y(i,j-1)) / (2 * dy);


   % края нулевые
end

