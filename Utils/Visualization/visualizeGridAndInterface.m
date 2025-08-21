function visualizeGridAndInterface(x, y, dx, dy, phase, c)
    % Визуализация фазовой сетки и интерфейса в пределах прямоугольной области

    [M, N] = size(phase);
    figure; hold on;

    % Отрисовка ячеек
    for ix = 1:N
        for iy = 1:M
            x1 = x(ix);
            y1 = y(iy);
            color = [0.8 0.8 1];  % жидкость
            if phase(iy, ix) == 1
                color = [1 0.8 0.8];  % пар
            end
            fill([x1 x1+dx x1+dx x1], [y1 y1 y1+dy y1+dy], color, 'EdgeColor', 'k');
        end
    end

    % Отображение границы
    x_dense = linspace(x(1), x(end), 1000);  % внутри области x
    y_poly = zeros(size(x_dense));
    for i = 0:length(c)-1
        y_poly = y_poly + c(i+1) * (x_dense.^i);
    end

    % Ограничение по y в пределах области y
    y_min = y(1);
    y_max = y(end);
    y_poly = max(min(y_poly, y_max), y_min);

    % Рисуем ограниченный интерфейс
    plot(x_dense, y_poly, 'k', 'LineWidth', 2);

    xlabel('x (м)');
    ylabel('y (м)');
    title('Сетка и интерфейс');
    axis equal tight;
    hold off;
end
% function visualizeGridAndInterface(x, y, dx, dy, phase, c)
%     % Визуализация фазовой сетки и интерфейса в пределах прямоугольной области
% 
%     [M, N] = size(phase);
%     figure; hold on;
% 
%     % Координаты углов всех ячеек
%     [X, Y] = meshgrid(x, y);
%     X = X(:);
%     Y = Y(:);
% 
%     % Каждая ячейка описывается 4 углами
%     xv = [X, X+dx, X+dx, X]';
%     yv = [Y, Y, Y+dy, Y+dy]';
% 
%     % Цвета: сразу вся матрица
%     C = zeros(3, M*N);
%     mask = phase(:)' == 1;
%     C(:,~mask) = repmat([0.8; 0.8; 1], 1, sum(~mask)); % жидкость
%     C(:, mask) = repmat([1; 0.8; 0.8], 1, sum(mask)); % пар
% 
%     % Рисуем все ячейки одним вызовом
%     patch(xv, yv, C', 'EdgeColor', 'k');
% 
%     % Интерфейс (полином)
%     x_dense = linspace(x(1), x(end), 1000);
%     y_poly = polyval(flip(c), x_dense);
% 
%     % Ограничение по y
%     y_poly = max(min(y_poly, y(end)), y(1));
% 
%     % Рисуем интерфейс
%     plot(x_dense, y_poly, 'k', 'LineWidth', 2);
% 
%     xlabel('x (м)');
%     ylabel('y (м)');
%     title('Сетка и интерфейс');
%     axis equal tight;
%     hold off;
% end
