function R = computeNSResidualY_full(fields, params, phase, dx, dy)
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

    % --- Центральные разности ---
    % d/dx и d/dy для u_y и u_x
    du_y_dx = (u_y([2:M, M], :) - u_y([1, 1:M-1], :)) / (2*dx);
    du_y_dy = (u_y(:, [2:N, N]) - u_y(:, [1, 1:N-1])) / (2*dy);
    du_x_dx = (u_x([2:M, M], :) - u_x([1, 1:M-1], :)) / (2*dx);
    du_x_dy = (u_x(:, [2:N, N]) - u_x(:, [1, 1:N-1])) / (2*dy);
    dP_dy   = (P(:, [2:N, N]) - P(:, [1, 1:N-1])) / (2*dy);

    % --- Конвекция ---
    conv = rho .* (u_x .* du_y_dx + u_y .* du_y_dy);

    % --- Вязкие члены ---
    % Δ_x (μ u_y)
    visc_x = (mu([2:M, M], :) .* u_y([2:M, M], :) ...
            - 2*mu .* u_y ...
            + mu([1, 1:M-1], :) .* u_y([1, 1:M-1], :)) / dx^2;

    % Δ_y (μ u_y)
    visc_y = (mu(:, [2:N, N]) .* u_y(:, [2:N, N]) ...
            - 2*mu .* u_y ...
            + mu(:, [1, 1:N-1]) .* u_y(:, [1, 1:N-1])) / dy^2;

    % Смешанные члены (∂/∂y[ μ ∂u_x/∂x ])
    cross_term = (mu(:, [2:N, N]) .* (u_x([2:M, M], [2:N, N]) - u_x([1, 1:M-1], [2:N, N])) ...
                - mu(:, [1, 1:N-1]) .* (u_x([2:M, M], [1, 1:N-1]) - u_x([1, 1:M-1], [1, 1:N-1])) ) ...
                / (4*dx*dy);

    % --- Итог ---
    R = conv + dP_dy - visc_x - (4/3)*visc_y + (2/3)*cross_term;
end

% 
% 
% 
% 






% 
% 
% 
% function R = computeNSResidualY_full(fields, params, phase, dx, dy)
%     % Векторизованная версия расчета невязки по Y в уравнении Навье–Стокса
% 
%     u_x = fields.u_x;
%     u_y = fields.u_y;
%     P   = fields.P;
% 
%     % Свойства
%     rho = params.rho_liquid * (1 - phase) + params.rho_vapor * phase;
%     mu  = params.mu_liquid  * (1 - phase) + params.mu_vapor  * phase;
% 
%     % --- Производные ---
%     du_y_dx = getScalarDerivativeFull(u_y, dx, dy, "x");
%     du_y_dy = getScalarDerivativeFull(u_y, dx, dy, "y");
%     du_x_dx = getScalarDerivativeFull(u_x, dx, dy, "x");
%     du_x_dy = getScalarDerivativeFull(u_x, dx, dy, "y");
%     dP_dy   = getScalarDerivativeFull(P,   dx, dy, "y");
% 
%     % --- Конвекция ---
%     conv = rho .* (u_x .* du_y_dx + u_y .* du_y_dy);
% 
%     % --- Вязкие члены ---
%     % ∂/∂x [ μ (∂u_x/∂y + ∂u_y/∂x) ]
%     term_x = du_x_dy + du_y_dx;
%     dvisc_x_dx = getProductDerivative(mu, term_x, dx, dy, "x");
% 
%     % ∂/∂y [ μ (4/3 ∂u_y/∂y - 2/3 ∂u_x/∂x) ]
%     term_y = (4/3)*du_y_dy - (2/3)*du_x_dx;
%     dvisc_y_dy = getProductDerivative(mu, term_y, dx, dy, "y");
% 
%     % --- Итоговая невязка ---
%     R = conv + dP_dy - dvisc_x_dx - dvisc_y_dy;
% end
% 
% 

%%%%%%%%% UNCOMMENT FOR UNVECTORIZED VERSION
% function R = computeNSResidualY_full(fields, params, phase, dx, dy)
%     % Вычисляет невязку по Y для уравнения Навье–Стокса
%     % R(i,j) — невязка в точке (i,j)
% 
%     [M, N] = size(phase);
%     R = zeros(M, N);
% 
%     % Извлекаем поля
%     u_x = fields.u_x;
%     u_y = fields.u_y;
%     P   = fields.P;
% 
%     % Распределение плотности и вязкости по фазам
%     rho = params.rho_liquid * (1 - phase) + params.rho_vapor * phase;
%     mu  = params.mu_liquid  * (1 - phase) + params.mu_vapor  * phase;
% 
% 
%     for i = 1:M
%         for j = 1:N
%             % --- Односторонние производные на границах ---
%             if j == 1
%                 du_y_dy = (u_y(i,j+1) - u_y(i,j)) / dy;
%                 du_x_dy = (u_x(i,j+1) - u_x(i,j)) / dy;
%             elseif j == N
%                 du_y_dy = (u_y(i,j) - u_y(i,j-1)) / dy;
%                 du_x_dy = (u_x(i,j) - u_x(i,j-1)) / dy;
%             else
%                 du_y_dy = (u_y(i,j+1) - u_y(i,j-1)) / (2*dy);
%                 du_x_dy = (u_x(i,j+1) - u_x(i,j-1)) / (2*dy);
%             end
% 
%             if i == 1
%                 du_y_dx = (u_y(i+1,j) - u_y(i,j)) / dx;
%                 du_x_dx = (u_x(i+1,j) - u_x(i,j)) / dx;
%                 dP_dy   = (P(i+1,j) - P(i,j)) / dy;  % исправлено на dy
%             elseif i == M
%                 du_y_dx = (u_y(i,j) - u_y(i-1,j)) / dx;
%                 du_x_dx = (u_x(i,j) - u_x(i-1,j)) / dx;
%                 dP_dy   = (P(i,j) - P(i-1,j)) / dy;  % исправлено на dy
%             else
%                 du_y_dx = (u_y(i+1,j) - u_y(i-1,j)) / (2*dx);
%                 du_x_dx = (u_x(i+1,j) - u_x(i-1,j)) / (2*dx);
%                 dP_dy   = (P(i+1,j) - P(i-1,j)) / (2*dy);  % исправлено на dy
%             end
% 
%             % --- Конвекция ---
%             conv = rho(i,j) * (u_x(i,j) * du_y_dx + u_y(i,j) * du_y_dy);
% 
%             % --- Вязкие члены ---
%             % ∂/∂x [ μ (∂u_x/∂y + ∂u_y/∂x) ]
%             if i == 1
%                 mu_x_plus = mu(i+1,j);
%                 mu_x_minus = mu(i,j);
%                 term_x_plus = du_x_dy + du_y_dx; % в точке i+1
%                 term_x = du_x_dy + du_y_dx; % в точке i
%                 visc_x = (mu_x_plus*term_x_plus - mu_x_minus*term_x) / dx;
%             elseif i == M
%                 mu_x_plus = mu(i,j);
%                 mu_x_minus = mu(i-1,j);
%                 term_x = du_x_dy + du_y_dx; % в точке i
%                 term_x_minus = du_x_dy + du_y_dx; % в точке i-1
%                 visc_x = (mu_x_plus*term_x - mu_x_minus*term_x_minus) / dx;
%             else
%                 mu_x_plus = mu(i+1,j);
%                 mu_x_minus = mu(i-1,j);
%                 term_x_plus = du_x_dy + du_y_dx; % в точке i+1
%                 term_x_minus = du_x_dy + du_y_dx; % в точке i-1
%                 visc_x = (mu_x_plus*term_x_plus - mu_x_minus*term_x_minus) / (2*dx);
%             end
% 
%             % ∂/∂y [ μ (4/3 ∂u_y/∂y - 2/3 ∂u_x/∂x) ]
%             term1 = (4/3) * du_y_dy - (2/3) * du_x_dx;
%             if j == 1
%                 mu_y_plus = mu(i,j+1);
%                 mu_y = mu(i,j);
%                 visc_y = (mu_y_plus * term1 - mu_y * term1) / dy;
%             elseif j == N
%                 mu_y = mu(i,j);
%                 mu_y_minus = mu(i,j-1);
%                 visc_y = (mu_y * term1 - mu_y_minus * term1) / dy;
%             else
%                 mu_y_plus = mu(i,j+1);
%                 mu_y_minus = mu(i,j-1);
%                 visc_y = (mu_y_plus * term1 - mu_y_minus * term1) / (2*dy);
%             end
% 
%             % --- Итоговая невязка ---
%             R(i,j) = conv + dP_dy - visc_x - visc_y;
%         end
%     end
% end
