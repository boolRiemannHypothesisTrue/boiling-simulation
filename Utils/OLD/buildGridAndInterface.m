% function [X, Y, phase, y_interface, x_nodes, y_nodes, dx, dy] = buildGridAndInterface(N, M, x_min, x_max, y_min, y_max, c)
%     % Создание двумерной сетки, фазовой матрицы и интерфейса по полиному
% 
%     % 1. Узлы сетки
%     x_nodes = linspace(x_min, x_max, N+1);
%     y_nodes = linspace(y_min, y_max, M+1);
%     dx = x_nodes(2) - x_nodes(1);
%     dy = y_nodes(2) - y_nodes(1);
% 
%     % 2. Центры ячеек
%     x_centers = x_nodes(1:end-1) + dx/2;
%     y_centers = y_nodes(1:end-1) + dy/2;
%     [X, Y] = meshgrid(x_centers, y_centers);
% 
%     % 3. Интерфейс как полином y = P(x)
%     y_interface = zeros(size(x_centers));
%     for i = 0:length(c)-1
%         y_interface = y_interface + c(i+1) * (x_centers.^i)';
%     end
% 
%     % 4. Ограничение интерфейса в пределах области
%     y_interface_clipped = max(min(y_interface, y_max), y_min);
% 
%     % Если была обрезка — предупредим
%     if any(abs(y_interface - y_interface_clipped) > 1e-10)
%         warning('Интерфейс был обрезан до области [%g, %g] по y.', y_min, y_max);
%     end
% 
%     y_interface = y_interface_clipped;
% 
%     % 5. Формируем фазовую матрицу: 0 — жидкость, 1 — пар
%     phase = zeros(M, N);
%     for ix = 1:N
%         for iy = 1:M
%             if Y(iy, ix) > y_interface(ix)
%                 phase(iy, ix) = 1;  % пар
%             else
%                 phase(iy, ix) = 0;  % жидкость
%             end
%         end
%     end
% end


function [X, Y, phase, y_interface,x_poly, x_nodes, y_nodes, dx, dy] = buildGridAndInterface(N, M, x_min, x_max, y_min, y_max, c)
    % Создание двумерной сетки, фазовой матрицы и интерфейса по полиному

    % 1. Узлы сетки
    x_nodes = linspace(x_min, x_max, N+1);
    y_nodes = linspace(y_min, y_max, M+1);
    dx = x_nodes(2) - x_nodes(1);
    dy = y_nodes(2) - y_nodes(1);

    % 2. Центры ячеек
    x_centers = x_nodes(1:end-1) + dx/2;
    y_centers = y_nodes(1:end-1) + dy/2;
    [X, Y] = meshgrid(x_centers, y_centers);

    % 3. Интерфейс как полином y = P(x)
    x_poly = x_centers';  % столбец
    y_interface = polyval(flip(c), x_poly);  % flip, т.к. MATLAB ожидает [a_n, ..., a_0]
    %y_interface = sin(x_poly);
    % 4. Ограничение интерфейса
    y_interface_clipped = max(min(y_interface, y_max), y_min);

    if any(abs(y_interface - y_interface_clipped) > 1e-10)
        warning('Интерфейс был обрезан до области [%g, %g] по y.', y_min, y_max);
    end

    y_interface = y_interface_clipped;

    % 5. Формирование фазовой матрицы: 0 — жидкость, 1 — пар
    % Расширяем y_interface в M строк (вдоль y), чтобы сравнить с Y
    y_interface_matrix = repmat(y_interface', M, 1);
    phase = Y > y_interface_matrix;
   % plot(x_poly,y_interface);
end
