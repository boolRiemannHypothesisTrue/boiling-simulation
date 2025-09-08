function val = computeNSResidualX_full(fields, params, phase, dx, dy)
    % невязка Х-проекции ур-й Навье Стокса
    % о векторизации можно пока забыть.

    % Извлекаем поля
    u_x = fields.u_x;
    u_y = fields.u_y;
    p   = fields.P;

    % Размеры
    [M, N] = size(phase);

    % Свойства среды
    rho = params.rho_liquid * (1 - phase) + params.rho_vapor * phase;
    mu  = params.mu_liquid  * (1 - phase) + params.mu_vapor  * phase;
    g   = params.gravity;
    val = zeros(M,N);
    i = 2:M-1; j = 2:N-1;

    val(i,j) = rho(i,j) * ( u_x(i,j) * ( u_x(i+1,j) - u_x(i-1,j) ) / (2*dx) + u_y(i,j) * ( u_x(i,j+1) - u_x(i,j-1) ) / (2*dy) + ...
        + g) + ( p(i+1,j) - p(i-1,j) ) / (2*dx) - ...
        (4/3) * ( mu(i+1,j) * ( u_x(i+1,j) - u_x(i,j) ) / dx^2 + mu(i,j) * ( u_x(i-1,j) - u_x(i,j) ) / dx^2 ) + ...
        (2/3) * ( ( mu(i+1,j) * u_y(i+1,j+1) - mu(i+1,j) * u_y(i+1,j-1) ) / (2*dx*2*dy) - ( mu(i-1,j) * u_y(i-1,j+1) - mu(i-1,j) * u_y(i-1,j-1) ) / (2*dx*2*dy) ) - ...
        ( mu(i,j+1) * ( u_x(i,j+1) - u_x(i,j) ) / dy^2 + mu(i,j) * ( u_x(i,j-1) - u_x(i,j) ) / dy^2 ) - ...
        ( mu(i,j+1) * ( u_y(i+1,j+1) - u_y(i-1,j+1) ) / (2*dx*2*dy) - mu(i,j-1) * ( u_y(i+1,j-1) - u_y(i-1,j-1) ) / (2*dx*2*dy) );

    % 
    % % Предрасчёт первых производных (вся сетка)
    % dux_dx = getScalarDerivativeFull(u_x, dx, dy, 'x');
    % dux_dy = getScalarDerivativeFull(u_x, dx, dy, 'y');
    % duy_dx = getScalarDerivativeFull(u_y, dx, dy, 'x');
    % duy_dy = getScalarDerivativeFull(u_y, dx, dy, 'y');
    % dp_dx  = getScalarDerivativeFull(p,   dx, dy, 'x');
    % 
    % % Вязкие элементы
    % term_x = (4/3) * dux_dx - (2/3) * duy_dy;
    % term_y = dux_dy + duy_dx;
    % 
    % % ∂/∂x [μ * term_x] и ∂/∂y [μ * term_y]
    % dvisc_x_dx = getScalarDerivativeFull(mu .* term_x, dx, dy, 'x');
    % dvisc_y_dy = getScalarDerivativeFull(mu .* term_y, dx, dy, 'y');
    % 
    % % Конвективный член и невязка
    % conv = rho .* (u_x .* dux_dx + u_y .* dux_dy + g);
    % val = conv + dp_dx - dvisc_x_dx - dvisc_y_dy;
end


