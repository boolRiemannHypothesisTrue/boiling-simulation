function R = computeHeatFluxJump(fields, params, phase, dx, dy)
    [M, N] = size(phase);
    T = fields.T;  % Температурное поле

    % Коэффициенты теплопроводности
    lambda = params.lambda_liquid * (phase == 0) + params.lambda_vapor * (phase == 1);

    % Градиенты температуры
    [dTdx, dTdy] = gradient(T, dx, dy);

    % Инициализация
    R = zeros(M, N);

    % Проходим по 4 соседям (x и y направления)
    shifts = {[-1, 0], [0, -1], [1, 0], [0, 1]};
    normals = {[0, -1], [-1, 0], [0, 1], [1, 0]};

    [I, J] = ndgrid(1:M, 1:N);

    for d = 1:4
        shift = shifts{d};
        normal = normals{d};
        ni = normal(1); nj = normal(2);

        I2 = I + shift(1);
        J2 = J + shift(2);

        inside = I2 >= 1 & I2 <= M & J2 >= 1 & J2 <= N;
        if ~any(inside(:)), continue; end

        idx1 = find(inside);
        idx2 = sub2ind([M N], I2(idx1), J2(idx1));

        p1 = phase(idx1);
        p2 = phase(idx2);
        is_interface = p1 ~= p2;

        gas_left  = is_interface & (p1 == 1); % газ слева
        gas_right = is_interface & (p2 == 1); % газ справа

        % Получаем компоненты ∇T·n
        gradT1_n = dTdx(idx1) * ni + dTdy(idx1) * nj;
        gradT2_n = dTdx(idx2) * ni + dTdy(idx2) * nj;

        k1 = lambda(idx1); k2 = lambda(idx2);

        % Потоки
        q1 = -k1 .* gradT1_n;
        q2 = -k2 .* gradT2_n;

        dQ = zeros(size(idx1));
        dQ(gas_left)  = q1(gas_left) - q2(gas_left);
        dQ(gas_right) = q2(gas_right) - q1(gas_right);

        R(idx1) = R(idx1) + dQ;
    end
end
