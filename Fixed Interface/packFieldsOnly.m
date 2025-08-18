function x = packFieldsOnly(fields, M, N)
    % Упаковка только полей (без коэффициентов интерфейса)

    total_points = M * N;

    % Векторизация полей
    ux_vec = reshape(fields.u_x, total_points, 1);
    uy_vec = reshape(fields.u_y, total_points, 1);
    T_vec  = reshape(fields.T,  total_points, 1);
    P_vec  = reshape(fields.P,  total_points, 1);

    % Собираем всё в один столбец
    x = [ux_vec; uy_vec; T_vec; P_vec];
end
