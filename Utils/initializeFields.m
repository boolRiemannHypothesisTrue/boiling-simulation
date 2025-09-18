function fields = initializeFields(X, Y, phase, params)
    % Инициализация полей u_x, u_y, T, P на сетке
    % phase: матрица 0/1 - жидкость/пар
    % params: структура с физическими параметрами и граничными условиями

    [M, N] = size(phase);

    % Инициализация
    u_x = zeros(M, N);
    u_y = zeros(M, N);
    T = zeros(M, N);
    P = zeros(M, N);

    % Для простоты зададим давление как постоянное (можно потом усложнить)
    P_init_liquid = 101325;  % давление жидкости, Па (пример)
    P_init_vapor = 101325;   % давление пара, Па (пример)

    for i = 1:M
        for j = 1:N
            if phase(i,j) == 0
                % Жидкая фаза
                T(i,j) = params.T_liquid + 10;
                P(i,j) = P_init_liquid;
            else
                % Паровая фаза
                T(i,j) = params.T_wall + 10;  % или можно задать другую температуру для пара
                P(i,j) = P_init_vapor;
            end
            % Скорости НЕнулевые но распределение однородное, везде
            u_x(i,j) = i + j;
            u_y(i,j) = j + abs(asin(i));
        end
    end

    % Возвращаем структуру с полями
    fields.u_x = u_x;
    fields.u_y = u_y;
    fields.T = T;
    fields.P = P;
end
