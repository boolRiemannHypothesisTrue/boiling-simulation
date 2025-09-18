function R = computeResidual_6(fields, params, phase, dx, dy)
    [M, N] = size(phase);
    rho = params.rho_liquid * (phase == 0) + params.rho_vapor * (phase == 1);
    ux = fields.u_x;
    uy = fields.u_y;

    R = zeros(M, N);

    shifts = {
        [-1, 0],  [0, -1],  [1, 0],  [0, 1]   % вверх, влево, вниз, вправо
    };
    normals = {
        [0, -1],  [-1, 0],  [0, 1],  [1, 0]
    };

    for d = 1:4
        shift = shifts{d};
        normal = normals{d};

        % Создаем сетку индексов
        [I, J] = ndgrid(1:M, 1:N);

        % Сдвинутые индексы с учётом границ
        I_shift = I + shift(1);
        J_shift = J + shift(2);

        % Маска в пределах массива
        valid = I_shift >= 1 & I_shift <= M & J_shift >= 1 & J_shift <= N;

        % Берем только существующие индексы
        I_valid = I(valid);
        J_valid = J(valid);
        I_shift_valid = I_shift(valid);
        J_shift_valid = J_shift(valid);

        % Фазы
        p1 = phase(sub2ind([M, N], I_valid, J_valid));
        p2 = phase(sub2ind([M, N], I_shift_valid, J_shift_valid));

        % Только на границе фаз
        mask_interface = (p1 ~= p2);

        if ~any(mask_interface)
            continue
        end

        % Фильтруем по интерфейсу
        I_int = I_valid(mask_interface);
        J_int = J_valid(mask_interface);
        I_shift_int = I_shift_valid(mask_interface);
        J_shift_int = J_shift_valid(mask_interface);

        p1_int = p1(mask_interface);
        p2_int = p2(mask_interface);

        % Определяем газ/жидкость
        is_gas1 = (p1_int == 1);
        is_gas2 = (p2_int == 1);

        % Плотности
        idx1 = sub2ind([M, N], I_int, J_int);
        idx2 = sub2ind([M, N], I_shift_int, J_shift_int);
        rho1 = rho(idx1);
        rho2 = rho(idx2);

        % Скорости
        ux1 = ux(idx1); uy1 = uy(idx1);
        ux2 = ux(idx2); uy2 = uy(idx2);

        un1 = ux1 * normal(1) + uy1 * normal(2);
        un2 = ux2 * normal(1) + uy2 * normal(2);

        J1 = rho1 .* un1;
        J2 = rho2 .* un2;

        % Вычисляем скачок - всегда газ - жидкость
        jump = zeros(size(J1));
        mask_gas_liq = is_gas1 & ~is_gas2;  % текущая газ, сосед жидкость
        mask_liq_gas = ~is_gas1 & is_gas2;  % текущая жидкость, сосед газ

        jump(mask_gas_liq) = J1(mask_gas_liq) - J2(mask_gas_liq);
        jump(mask_liq_gas) = J2(mask_liq_gas) - J1(mask_liq_gas);

        % Накопление невязки
        ind_linear = sub2ind([M, N], I_int, J_int);
        R(ind_linear) = R(ind_linear) + jump;
    end
end
