function R = computeEnergyResidual(fields, params, phase, dx, dy)
    % Векторизованная невязка уравнения энергии

    u_x = fields.u_x;
    u_y = fields.u_y;
    T   = fields.T;

    [M, N] = size(phase);

    % Параметры по фазам
    rho = params.rho_liquid * (1 - phase) + params.rho_vapor * phase;
    Cp = params.Cp_liquid * (1 - phase) + params.Cp_vapor * phase;
    lambda = params.lambda_liquid * (1 - phase) + params.lambda_vapor * phase;
    mu = params.mu_liquid * (1 - phase) + params.mu_vapor * phase;

    % Производные температуры
    dT_dx = getScalarDerivativeFull(T, dx, dy, "x");
    dT_dy = getScalarDerivativeFull(T, dx, dy, "y");

    % Производные скоростей
    dux_dx = getScalarDerivativeFull(u_x, dx, dy, "x");
    duy_dx = getScalarDerivativeFull(u_y, dx, dy, "x");
    dux_dy = getScalarDerivativeFull(u_x, dx, dy, "y");
    duy_dy = getScalarDerivativeFull(u_y, dx, dy, "y");

    % Конвекция
    conv = rho .* Cp .* (u_x .* dT_dx + u_y .* dT_dy);

    % Диффузионные члены: ∂/∂x(λ ∂T/∂x), ∂/∂y(λ ∂T/∂y)
    term_x = lambda .* dT_dx;
    diff_x = getScalarDerivativeFull(term_x, dx, dy, "x");

    term_y = lambda .* dT_dy;
    diff_y = getScalarDerivativeFull(term_y, dx, dy, "y");

    % Вязкие диссипационные члены
    shear_1 = dux_dx + duy_dy;
    shear_2 = dux_dy + duy_dx;
    visc_term = mu .* ( (4/3) * shear_1.^2 + shear_2.^2 - (16/3) * dux_dx .* duy_dy );

    % Итог
    R = conv - diff_x - diff_y - visc_term;
end


% UNCOMMENT FOR UNVECTORIZED VERSION

% function R = computeEnergyResidual(fields, params, phase, dx, dy)
%     % Вычисляет невязку уравнения энергии Q_T
%     % fields: структура с полями u_x, u_y, T
%     % params: структура с параметрами rho_liquid, rho_vapor, Cp_liquid, Cp_vapor, lambda_liquid, lambda_vapor, mu_liquid, mu_vapor
%     % phase: матрица фаз (0 или 1)
%     % dx, dy: шаги сетки
% 
%     [M, N] = size(phase);
%     R = zeros(M, N);
% 
%     % Извлекаем поля
%     u_x = fields.u_x;
%     u_y = fields.u_y;
%     T   = fields.T;
% 
%     % Параметры по фазам
%     rho = params.rho_liquid * (1 - phase) + params.rho_vapor * phase;
%     Cp  = params.Cp_liquid  * (1 - phase) + params.Cp_vapor  * phase;
%     lambda = params.lambda_liquid * (1 - phase) + params.lambda_vapor * phase;
%     mu = params.mu_liquid * (1 - phase) + params.mu_vapor * phase;
% 
%     for i = 1:M
%         for j = 1:N
%             % Производные T
%             if j == 1
%                 dT_dx = (T(i,j+1) - T(i,j)) / dx;
%             elseif j == N
%                 dT_dx = (T(i,j) - T(i,j-1)) / dx;
%             else
%                 dT_dx = (T(i,j+1) - T(i,j-1)) / (2*dx);
%             end
% 
%             if i == 1
%                 dT_dy = (T(i+1,j) - T(i,j)) / dy;
%             elseif i == M
%                 dT_dy = (T(i,j) - T(i-1,j)) / dy;
%             else
%                 dT_dy = (T(i+1,j) - T(i-1,j)) / (2*dy);
%             end
% 
%             % Производные скоростей (частные производные)
%             % ∂(u_x)/∂x, ∂(u_y)/∂y, ∂(u_x)/∂y, ∂(u_y)/∂x
% 
%             if j == 1
%                 dux_dx = (u_x(i,j+1) - u_x(i,j)) / dx;
%                 duy_dx = (u_y(i,j+1) - u_y(i,j)) / dx;
%             elseif j == N
%                 dux_dx = (u_x(i,j) - u_x(i,j-1)) / dx;
%                 duy_dx = (u_y(i,j) - u_y(i,j-1)) / dx;
%             else
%                 dux_dx = (u_x(i,j+1) - u_x(i,j-1)) / (2*dx);
%                 duy_dx = (u_y(i,j+1) - u_y(i,j-1)) / (2*dx);
%             end
% 
%             if i == 1
%                 dux_dy = (u_x(i+1,j) - u_x(i,j)) / dy;
%                 duy_dy = (u_y(i+1,j) - u_y(i,j)) / dy;
%             elseif i == M
%                 dux_dy = (u_x(i,j) - u_x(i-1,j)) / dy;
%                 duy_dy = (u_y(i,j) - u_y(i-1,j)) / dy;
%             else
%                 dux_dy = (u_x(i+1,j) - u_x(i-1,j)) / (2*dy);
%                 duy_dy = (u_y(i+1,j) - u_y(i-1,j)) / (2*dy);
%             end
% 
%             % Конвекционный член
%             conv = rho(i,j) * Cp(i,j) * (u_x(i,j)*dT_dx + u_y(i,j)*dT_dy);
% 
%             % Вязкие диффузионные члены: - ∂/∂x(λ ∂T/∂x) - ∂/∂y(λ ∂T/∂y)
%             % Считаем по разностной схеме с учетом границ
% 
%             % ∂/∂x( λ ∂T/∂x )
%             if j == 1
%                 dT_dx_plus = (T(i,j+1) - T(i,j)) / dx;
%                 lambda_x_plus = lambda(i,j+1);
%                 lambda_x = lambda(i,j);
%                 diff_x = (lambda_x_plus * dT_dx_plus - lambda_x * dT_dx) / dx;
%             elseif j == N
%                 dT_dx_minus = (T(i,j) - T(i,j-1)) / dx;
%                 lambda_x = lambda(i,j);
%                 lambda_x_minus = lambda(i,j-1);
%                 diff_x = (lambda_x * dT_dx - lambda_x_minus * dT_dx_minus) / dx;
%             else
%                 dT_dx_plus = (T(i,j+1) - T(i,j)) / dx;
%                 dT_dx_minus = (T(i,j) - T(i,j-1)) / dx;
%                 lambda_x_plus = lambda(i,j+1);
%                 lambda_x_minus = lambda(i,j-1);
%                 diff_x = (lambda_x_plus * dT_dx_plus - lambda_x_minus * dT_dx_minus) / (2*dx);
%             end
% 
%             % ∂/∂y( λ ∂T/∂y )
%             if i == 1
%                 dT_dy_plus = (T(i+1,j) - T(i,j)) / dy;
%                 lambda_y_plus = lambda(i+1,j);
%                 lambda_y = lambda(i,j);
%                 diff_y = (lambda_y_plus * dT_dy_plus - lambda_y * dT_dy) / dy;
%             elseif i == M
%                 dT_dy_minus = (T(i,j) - T(i-1,j)) / dy;
%                 lambda_y = lambda(i,j);
%                 lambda_y_minus = lambda(i-1,j);
%                 diff_y = (lambda_y * dT_dy - lambda_y_minus * dT_dy_minus) / dy;
%             else
%                 dT_dy_plus = (T(i+1,j) - T(i,j)) / dy;
%                 dT_dy_minus = (T(i,j) - T(i-1,j)) / dy;
%                 lambda_y_plus = lambda(i+1,j);
%                 lambda_y_minus = lambda(i-1,j);
%                 diff_y = (lambda_y_plus * dT_dy_plus - lambda_y_minus * dT_dy_minus) / (2*dy);
%             end
% 
%             % Вязкие диссипационные члены:
%             shear_1 = (dux_dx + duy_dy);
%             shear_2 = (dux_dy + duy_dx);
% 
%             visc_term = mu(i,j) * ( (4/3)*shear_1^2 + shear_2^2 - (16/3) * dux_dx * duy_dy ); % 16/3 вместо 4
% 
%             % Итоговая невязка
%             R(i,j) = conv - diff_x - diff_y - visc_term;
%         end
%     end
% end