function R = computeEnergyResidual(fields, params, phase, dx, dy)
    % Вычисляет невязку уравнения энергии Q_T
    % fields: структура с полями u_x, u_y, T
    % params: структура с параметрами rho_liquid, rho_vapor, Cp_liquid, Cp_vapor, lambda_liquid, lambda_vapor, mu_liquid, mu_vapor
    % phase: матрица фаз (0 или 1)
    % dx, dy: шаги сетки
    
        [M, N] = size(phase);
        R = zeros(M, N);
    
        % Извлекаем поля
        u_x = fields.u_x;
        u_y = fields.u_y;
        T   = fields.T;
    
        % Параметры по фазам
        rho = params.rho_liquid * (1 - phase) + params.rho_vapor * phase;
        Cp  = params.Cp_liquid  * (1 - phase) + params.Cp_vapor  * phase;
        lambda = params.lambda_liquid * (1 - phase) + params.lambda_vapor * phase;
        mu = params.mu_liquid * (1 - phase) + params.mu_vapor * phase;
        
        for i = 2:M-1
            for j=2:N-1
                R(i,j) = rho(i,j) * Cp(i,j) * ( u_x(i,j) * ( T(i+1,j) - T(i-1,j) ) / (2 * dx) + ...
                    u_y(i,j) * ( T(i,j+1) - T(i,j-1) ) / (2 * dy) ) - ...
                    ( lambda(i+1,j) * ( T(i+1,j) - T(i,j) ) / dx^2 + lambda(i,j) * ( T(i-1,j) - T(i,j) ) / dx^2 ) - ...
                    ( lambda(i,j+1) * ( T(i,j+1) - T(i,j) ) / dy^2 + lambda(i,j) * ( T(i,j-1) - T(i,j) ) / dy^2 ) - ...
                    (4/3) * mu(i,j) * ( ( u_x(i+1,j) - u_x(i-1,j) ) / (2 * dx) + (u_y(i,j+1) - u_y(i,j-1) ) / (2 * dy) )^2 - ...
                    mu(i,j) * ( (u_y(i+1,j) - u_y(i-1,j) ) / (2 * dx) + ( (u_x(i,j+1) - u_x(i,j-1) ) / (2 * dy) ) )^2 + ...
                    4 * mu(i,j) *  ( u_x(i+1,j) - u_x(i-1,j) ) * ( u_y(i,j+1) - u_y(i,j-1) ) / (4 * dx * dy); 

            end
        end



        
end


% function R = computeEnergyResidual(fields, params, phase, dx, dy)
%     % Векторизованная невязка уравнения энергии
% 
%     u_x = fields.u_x;
%     u_y = fields.u_y;
%     T   = fields.T;
% 
%     [M, N] = size(phase);
% 
%     % Параметры по фазам
%     rho = params.rho_liquid * (1 - phase) + params.rho_vapor * phase;
%     Cp = params.Cp_liquid * (1 - phase) + params.Cp_vapor * phase;
%     lambda = params.lambda_liquid * (1 - phase) + params.lambda_vapor * phase;
%     mu = params.mu_liquid * (1 - phase) + params.mu_vapor * phase;
% 
%     % Производные температуры
%     dT_dx = getScalarDerivativeFull(T, dx, dy, "x");
%     dT_dy = getScalarDerivativeFull(T, dx, dy, "y");
% 
%     % Производные скоростей
%     dux_dx = getScalarDerivativeFull(u_x, dx, dy, "x");
%     duy_dx = getScalarDerivativeFull(u_y, dx, dy, "x");
%     dux_dy = getScalarDerivativeFull(u_x, dx, dy, "y");
%     duy_dy = getScalarDerivativeFull(u_y, dx, dy, "y");
% 
%     % Конвекция
%     conv = rho .* Cp .* (u_x .* dT_dx + u_y .* dT_dy);
% 
%     % Диффузионные члены: ∂/∂x(λ ∂T/∂x), ∂/∂y(λ ∂T/∂y)
%     term_x = lambda .* dT_dx;
%     diff_x = getScalarDerivativeFull(term_x, dx, dy, "x");
% 
%     term_y = lambda .* dT_dy;
%     diff_y = getScalarDerivativeFull(term_y, dx, dy, "y");
% 
%     % Вязкие диссипационные члены
%     shear_1 = dux_dx + duy_dy;
%     shear_2 = dux_dy + duy_dx;
%     visc_term = mu .* ( (4/3) * shear_1.^2 + shear_2.^2 - (16/3) * dux_dx .* duy_dy );
% 
%     % Итог
%     R = conv - diff_x - diff_y - visc_term;
% end