function R = computeContinuityResidual(fields, params, phase, dx, dy)
    % Вычисление невязки уравнения непрерывности div(rho * u) = 0
    % Векторизованная версия

    u_x = fields.u_x;
    u_y = fields.u_y;

    % Плотность по фазе
    rho = params.rho_liquid * (1 - phase) + params.rho_vapor * phase;

    % Потоки импульса
    rhou_x = rho .* u_x;
    rhou_y = rho .* u_y;

    % Производные
    drhou_dx = getScalarDerivativeFull(rhou_x, dx, dy, 'x');
    drhou_dy = getScalarDerivativeFull(rhou_y, dx, dy, 'y');

    % Суммарная дивергенция
    R = drhou_dx + drhou_dy;
end

% function R = computeContinuityResidual(fields, params, phase, dx, dy)
%     % Вычисление невязки уравнения непрерывности div(rho * u) = 0
%     % Только во внутренней области, края = 0
% 
%     u_x = fields.u_x;
%     u_y = fields.u_y;
% 
%     % Плотность
%     rho = params.rho_liquid * (1 - phase) + params.rho_vapor * phase;
% 
%     % Потоки импульса
%     rhou_x = rho .* u_x;
%     rhou_y = rho .* u_y;
% 
%     % Размеры
%     [M, N] = size(rho);
% 
%     % Результат
%     R = zeros(M, N);
% 
%     % Центральные разности только по внутренним точкам
%     R(2:M-1, 2:N-1) = ...
%         ( rhou_x(2:M-1, 3:N) - rhou_x(2:M-1, 1:N-2) ) / (2*dx) + ...
%         ( rhou_y(3:M, 2:N-1) - rhou_y(1:M-2, 2:N-1) ) / (2*dy);
% 
%     % края остаются нулями
% end
