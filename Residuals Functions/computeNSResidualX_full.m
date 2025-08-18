% function R = computeNSResidualX_full(fields, params, phase, dx, dy)
%     % Расчет невязки по x-компоненте уравнения Навье–Стокса
%     % для структурированной равномерной сетки (центральные разности)
% 
%     % Извлекаем поля
%     u_x = fields.u_x;
%     u_y = fields.u_y;
%     p   = fields.P;
% 
%     % Размеры
%     [M, N] = size(phase);
% 
%     % Свойства среды
%     rho = params.rho_liquid * (1 - phase) + params.rho_vapor * phase;
%     mu  = params.mu_liquid  * (1 - phase) + params.mu_vapor  * phase;
%     g   = params.gravity;
% 
%     % Центральные разности (с внутренними точками)
%     dux_dx = (u_x(3:M,2:N-1) - u_x(1:M-2,2:N-1)) / (2*dx);
%     dux_dy = (u_x(2:M-1,3:N) - u_x(2:M-1,1:N-2)) / (2*dy);
%     duy_dx = (u_y(3:M,2:N-1) - u_y(1:M-2,2:N-1)) / (2*dx);
%     dp_dx  = (p(3:M,2:N-1)   - p(1:M-2,2:N-1))   / (2*dx);
% 
%     % Ограничиваем поля до внутренних точек
%     rho_c = rho(2:M-1,2:N-1);
%     mu_c  = mu(2:M-1,2:N-1);
%     u_xc  = u_x(2:M-1,2:N-1);
%     u_yc  = u_y(2:M-1,2:N-1);
% 
%     % Конвекция + источник g
%     conv = rho_c .* (u_xc .* dux_dx + u_yc .* dux_dy + g);
% 
%     % Вязкие члены (лапласиан u_x)
%     d2ux_dx2 = (mu(3:M,2:N-1).*u_x(3:M,2:N-1) - 2*mu_c.*u_xc + mu(1:M-2,2:N-1).*u_x(1:M-2,2:N-1)) / dx^2;
%     d2ux_dy2 = (mu(2:M-1,3:N).*u_x(2:M-1,3:N) - 2*mu_c.*u_xc + mu(2:M-1,1:N-2).*u_x(2:M-1,1:N-2)) / dy^2;
% 
%     visc_x = -(4/3) * (d2ux_dx2 + d2ux_dy2);
% 
%     % Смешанные члены (связь с u_y)
%     cross = (2/3) * ( ...
%         (mu(3:M,2:N-1).*(u_y(3:M,3:N) - u_y(3:M,1:N-2)) ...
%        -mu(1:M-2,2:N-1).*(u_y(1:M-2,3:N) - u_y(1:M-2,1:N-2))) / (4*dx*dy) );
% 
%     % Итоговая невязка
%     R = conv + dp_dx + visc_x + cross;
% end







function R = computeNSResidualX_full(fields, params, phase, dx, dy)
    % Векторизованный расчет невязки по x-компоненте уравнения Навье–Стокса

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

    % Предрасчёт первых производных (вся сетка)
    dux_dx = getScalarDerivativeFull(u_x, dx, dy, 'x');
    dux_dy = getScalarDerivativeFull(u_x, dx, dy, 'y');
    duy_dx = getScalarDerivativeFull(u_y, dx, dy, 'x');
    duy_dy = getScalarDerivativeFull(u_y, dx, dy, 'y');
    dp_dx  = getScalarDerivativeFull(p,   dx, dy, 'x');

    % Вязкие элементы
    term_x = (4/3) * dux_dx - (2/3) * duy_dy;
    term_y = dux_dy + duy_dx;

    % ∂/∂x [μ * term_x] и ∂/∂y [μ * term_y]
    dvisc_x_dx = getScalarDerivativeFull(mu .* term_x, dx, dy, 'x');
    dvisc_y_dy = getScalarDerivativeFull(mu .* term_y, dx, dy, 'y');

    % Конвективный член и невязка
    conv = rho .* (u_x .* dux_dx + u_y .* dux_dy + g);
    R = conv + dp_dx - dvisc_x_dx - dvisc_y_dy;
end


%%%%%%%%%%%%% UNCOMMENT FOR UNVECTORIZED VERSION
% function R = computeNSResidualX_full(fields, params, phase, dx, dy)
%     % Расчет невязки по x-компоненте уравнения Навье–Стокса
%     % с учетом переменных свойств и односторонних разностей на границе
% 
%     [M, N] = size(phase);
%     R = zeros(M, N);
% 
%     % Извлекаем поля
%     u_x = fields.u_x;
%     u_y = fields.u_y;
%     p   = fields.P;
% 
%     % Свойства среды с учетом фазы (phase в диапазоне [0,1])
%     rho = params.rho_liquid * (1 - phase) + params.rho_vapor * phase;
%     mu  = params.mu_liquid  * (1 - phase) + params.mu_vapor  * phase;
%     g   = params.gravity;
% 
%     for i = 1:M
%         for j = 1:N
%             % --- Производные по x ---
%             if j == 1
%                 dux_dx = (u_x(i, j+1) - u_x(i, j)) / dx;
%                 dp_dx  = (p(i, j+1) - p(i, j)) / dx;
%                 duy_dx = (u_y(i, j+1) - u_y(i, j)) / dx;
%             elseif j == N
%                 dux_dx = (u_x(i, j) - u_x(i, j-1)) / dx;
%                 dp_dx  = (p(i, j) - p(i, j-1)) / dx;
%                 duy_dx = (u_y(i, j) - u_y(i, j-1)) / dx;
%             else
%                 dux_dx = (u_x(i, j+1) - u_x(i, j-1)) / (2*dx);
%                 dp_dx  = (p(i, j+1) - p(i, j-1)) / (2*dx);
%                 duy_dx = (u_y(i, j+1) - u_y(i, j-1)) / (2*dx);
%             end
% 
%             % --- Производные по y ---
%             if i == 1
%                 dux_dy = (u_x(i+1, j) - u_x(i, j)) / dy;
%                 duy_dy = (u_y(i+1, j) - u_y(i, j)) / dy;
%             elseif i == M
%                 dux_dy = (u_x(i, j) - u_x(i-1, j)) / dy;
%                 duy_dy = (u_y(i, j) - u_y(i-1, j)) / dy;
%             else
%                 dux_dy = (u_x(i+1, j) - u_x(i-1, j)) / (2*dy);
%                 duy_dy = (u_y(i+1, j) - u_y(i-1, j)) / (2*dy);
%             end
% 
%             % --- Вязкие члены ---
%             % ∂/∂x [μ(4/3 ∂u_x/∂x - 2/3 ∂u_y/∂y)]
%             if j == 1
%                 term_x_plus = (4/3)*((u_x(i, j+1) - u_x(i, j))/dx) - (2/3)*((u_y(i, j+1) - u_y(i, j))/dx);
%                 term_x = (4/3)*dux_dx - (2/3)*duy_dy;
%                 dvisc_x_dx = (mu(i, j+1)*term_x_plus - mu(i, j)*term_x) / dx;
%             elseif j == N
%                 term_x = (4/3)*dux_dx - (2/3)*duy_dy;
%                 term_x_minus = (4/3)*((u_x(i, j) - u_x(i, j-1))/dx) - (2/3)*((u_y(i, j) - u_y(i, j-1))/dx);
%                 dvisc_x_dx = (mu(i, j)*term_x - mu(i, j-1)*term_x_minus) / dx;
%             else
%                 term_x_plus = (4/3)*((u_x(i, j+1) - u_x(i, j))/dx) - (2/3)*((u_y(i, j+1) - u_y(i, j))/dx);
%                 term_x_minus = (4/3)*((u_x(i, j) - u_x(i, j-1))/dx) - (2/3)*((u_y(i, j) - u_y(i, j-1))/dx);
%                 dvisc_x_dx = (mu(i, j+1)*term_x_plus - mu(i, j-1)*term_x_minus) / (2*dx);
%             end
% 
%             % ∂/∂y [μ(∂u_x/∂y + ∂u_y/∂x)]
%             if i == 1
%                 term_y_plus = ((u_x(i+1, j) - u_x(i, j))/dy) + ((u_y(i+1, j) - u_y(i, j))/dy);
%                 term_y = dux_dy + duy_dx;
%                 dvisc_y_dy = (mu(i+1, j)*term_y_plus - mu(i, j)*term_y) / dy;
%             elseif i == M
%                 term_y = dux_dy + duy_dx;
%                 term_y_minus = ((u_x(i, j) - u_x(i-1, j))/dy) + ((u_y(i, j) - u_y(i-1, j))/dy);
%                 dvisc_y_dy = (mu(i, j)*term_y - mu(i-1, j)*term_y_minus) / dy;
%             else
%                 term_y_plus = ((u_x(i+1, j) - u_x(i, j))/dy) + ((u_y(i+1, j) - u_y(i, j))/dy);
%                 term_y_minus = ((u_x(i, j) - u_x(i-1, j))/dy) + ((u_y(i, j) - u_y(i-1, j))/dy);
%                 dvisc_y_dy = (mu(i+1, j)*term_y_plus - mu(i-1, j)*term_y_minus) / (2*dy);
%             end
% 
%             % --- Итоговая невязка ---
%             R(i, j) = rho(i, j)*(u_x(i, j)*dux_dx + u_y(i, j)*dux_dy + g) + dp_dx - dvisc_x_dx - dvisc_y_dy;
%         end
%     end
% end