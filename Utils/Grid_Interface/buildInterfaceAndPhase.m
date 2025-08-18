function [y_interface, phase] = buildInterfaceAndPhase(X, Y, x_interface, c, y_min, y_max)
    % Строит интерфейс и фазовую матрицу по полиному

    y_interface = polyval(flip(c), x_interface);

    % Ограничиваем по y
    y_interface = max(min(y_interface, y_max), y_min);

    % Формируем фазовую матрицу
    y_interface_matrix = repmat(y_interface', size(Y, 1), 1);
    phase = Y > y_interface_matrix;
end
