function R = computeResidual_7(fields, params, phase, H, dx, dy)
    [M, N] = size(phase);

    % Физические поля
    rho = params.rho_liquid * (phase == 0) + params.rho_vapor * (phase == 1);
    mu  = params.mu_liquid  * (phase == 0) + params.mu_vapor  * (phase == 1);
    P   = fields.P;
    ux  = fields.u_x;
    uy  = fields.u_y;

    % Инициализация невязки
    R = zeros(M, N);

    % Направления и нормали
    shifts  = {[-1, 0], [0, -1], [1, 0], [0, 1]};
    normals = {[0, -1], [-1, 0], [0, 1], [1, 0]};

    % Производные для расчёта τ
    [dudx, dudy] = gradient(ux, dx, dy);
    [dvdx, dvdy] = gradient(uy, dx, dy);
    div_u = dudx + dvdy;

    % Координатная сетка
    [I, J] = ndgrid(1:M, 1:N);

    for d = 1:4
        shift  = shifts{d};
        normal = normals{d};
        ni = normal(1); nj = normal(2);

        I2 = I + shift(1);
        J2 = J + shift(2);

        inside = I2 >= 1 & I2 <= M & J2 >= 1 & J2 <= N;
        if ~any(inside(:)), continue; end

        % Инициализация промежуточных величин
        Jv = zeros(M, N);
        delta_u = zeros(M, N);
        delta_P = zeros(M, N);
        delta_tau = zeros(M, N);

        % Только те элементы, где сосед допустим
        idx = find(inside);

        % Текущая и соседняя фаза
        p1 = phase(idx);
        idx2 = sub2ind([M, N], I2(idx), J2(idx));
        p2 = phase(idx2);

        % Интерфейс и кто газ
        is_interface = p1 ~= p2;
        is_gas1 = p1 == 1;
        is_gas2 = p2 == 1;

        % Маски
        gas_left = is_interface & is_gas1 & ~is_gas2;
        gas_right = is_interface & ~is_gas1 & is_gas2;

        % Индексы текущей и соседней ячеек
        idx1 = idx;
        idx2 = sub2ind([M, N], I2(idx), J2(idx));

        % Векторные компоненты
        ux1 = ux(idx1); uy1 = uy(idx1);
        ux2 = ux(idx2); uy2 = uy(idx2);
        u1n = ux1 * ni + uy1 * nj;
        u2n = ux2 * ni + uy2 * nj;

        rho1 = rho(idx1); rho2 = rho(idx2);
        mu1  = mu(idx1);  mu2  = mu(idx2);
        P1   = P(idx1);   P2   = P(idx2);

        % Градиенты и вязкие компоненты
        div1 = div_u(idx1); div2 = div_u(idx2);
        tau1 = mu1 .* ((4/3) * div1 * ni^2) + 2 * mu1 .* ...
            viscous_proj(dudx(idx1), dudy(idx1), dvdx(idx1), dvdy(idx1), ni, nj);
        tau2 = mu2 .* ((4/3) * div2 * ni^2) + 2 * mu2 .* ...
            viscous_proj(dudx(idx2), dudy(idx2), dvdx(idx2), dvdy(idx2), ni, nj);

        % Запись в массивы
        tmp = zeros(size(idx));
        tmp(gas_left) = rho1(gas_left) .* u1n(gas_left);
        Jv(idx1) = Jv(idx1) + tmp;

        tmp(gas_left) = u1n(gas_left) - u2n(gas_left);
        delta_u(idx1) = delta_u(idx1) + tmp;

        tmp(gas_left) = P1(gas_left) - P2(gas_left);
        delta_P(idx1) = delta_P(idx1) + tmp;

        tmp(gas_left) = tau1(gas_left) - tau2(gas_left);
        delta_tau(idx1) = delta_tau(idx1) + tmp;

        % Газ справа
        tmp(gas_right) = rho2(gas_right) .* u2n(gas_right);
        Jv(idx1) = Jv(idx1) + tmp;

        tmp(gas_right) = u2n(gas_right) - u1n(gas_right);
        delta_u(idx1) = delta_u(idx1) + tmp;

        tmp(gas_right) = P2(gas_right) - P1(gas_right);
        delta_P(idx1) = delta_P(idx1) + tmp;

        tmp(gas_right) = tau2(gas_right) - tau1(gas_right);
        delta_tau(idx1) = delta_tau(idx1) + tmp;

        % Итоговая невязка за это направление
        R = R + (Jv .* delta_u + delta_P - delta_tau);
    end

    % Добавим член поверхностного натяжения
    R = R - 2 * params.sigma * H;
end

function T = viscous_proj(dudx, dudy, dvdx, dvdy, ni, nj)
    T = ni.^2 .* dudx + nj.^2 .* dvdy + ni .* nj .* (dudy + dvdx);
end
