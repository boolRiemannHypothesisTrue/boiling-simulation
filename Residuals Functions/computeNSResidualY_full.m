function val = computeNSResidualY_full(fields, params, phase, dx, dy)
    % Невязка по Y-компоненте уравнения Навье–Стокса (структурированная сетка)

    % Извлекаем поля
    u_x = fields.u_x;
    u_y = fields.u_y;
    P   = fields.P;

    % Свойства
    rho = params.rho_liquid * (1 - phase) + params.rho_vapor * phase;
    mu  = params.mu_liquid  * (1 - phase) + params.mu_vapor  * phase;

    % Размеры
    [M, N] = size(u_y);

    val = zeros(M,N); 

    i = 2:M;
    j = 2:N;


    
end





