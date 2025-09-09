function val = computeEnergyResidual(fields, params, phase, dx, dy)
    % Вычисляет невязку уравнения энергии Q_T
    % fields: структура с полями u_x, u_y, T
    % params: структура с параметрами rho_liquid, rho_vapor, Cp_liquid, Cp_vapor, lambda_liquid, lambda_vapor, mu_liquid, mu_vapor
    % phase: матрица фаз (0 или 1)
    % dx, dy: шаги сетки
    
        [M, N] = size(phase);
        val = zeros(M, N);
    
        % Извлекаем поля
        u_x = fields.u_x;
        u_y = fields.u_y;
        T   = fields.T;
    
        % Параметры по фазам
        rho = params.rho_liquid * (1 - phase) + params.rho_vapor * phase;
        Cp  = params.Cp_liquid  * (1 - phase) + params.Cp_vapor  * phase;
        lambda = params.lambda_liquid * (1 - phase) + params.lambda_vapor * phase;
        mu = params.mu_liquid * (1 - phase) + params.mu_vapor * phase;
        
        i = 2:M-1;
        j=2:N-1;

        val(i,j) = rho(i,j) * Cp(i,j) * ( u_x(i,j) * ( T(i+1,j) - T(i-1,j) ) / (2 * dx) + ...
            u_y(i,j) * ( T(i,j+1) - T(i,j-1) ) / (2 * dy) ) - ...
            ( lambda(i+1,j) * ( T(i+1,j) - T(i,j) ) / dx^2 + lambda(i,j) * ( T(i-1,j) - T(i,j) ) / dx^2 ) - ...
            ( lambda(i,j+1) * ( T(i,j+1) - T(i,j) ) / dy^2 + lambda(i,j) * ( T(i,j-1) - T(i,j) ) / dy^2 ) - ...
            (4/3) * mu(i,j) * ( ( u_x(i+1,j) - u_x(i-1,j) ) / (2 * dx) + (u_y(i,j+1) - u_y(i,j-1) ) / (2 * dy) )^2 - ...
            mu(i,j) * ( (u_y(i+1,j) - u_y(i-1,j) ) / (2 * dx) + ( (u_x(i,j+1) - u_x(i,j-1) ) / (2 * dy) ) )^2 + ...
            (16/3) * mu(i,j) *  ( u_x(i+1,j) - u_x(i-1,j) ) * ( u_y(i,j+1) - u_y(i,j-1) ) / (2 * dx * 2 * dy); 


        
end
