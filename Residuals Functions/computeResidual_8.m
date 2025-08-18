function R = computeResidual_8(fields, params, phase, dx, dy)
    [M, N] = size(phase);

    % Физические параметры
    rho = params.rho_liquid * (phase == 0) + params.rho_vapor * (phase == 1);
    mu  = params.mu_liquid  * (phase == 0) + params.mu_vapor  * (phase == 1);

    ux = fields.u_x;
    uy = fields.u_y;

    % Производные для расчета τ
    [dudx, dudy] = gradient(ux, dx, dy);
    [dvdx, dvdy] = gradient(uy, dx, dy);

    R = zeros(M, N);

    % Направления и нормали
    shifts  = {[-1, 0], [0, -1], [1, 0], [0, 1]};
    normals = {[0, -1], [-1, 0], [0, 1], [1, 0]};

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

        idx1 = sub2ind([M, N], I(inside), J(inside));
        idx2 = sub2ind([M, N], I2(inside), J2(inside));

        p1 = phase(idx1);
        p2 = phase(idx2);

        is_interface = p1 ~= p2;
        gas_left  = is_interface & (p1 == 1);
        gas_right = is_interface & (p2 == 1);

        % Скорости
        ux1 = ux(idx1); uy1 = uy(idx1);
        ux2 = ux(idx2); uy2 = uy(idx2);

        % Нормальные компоненты скорости
        u1n = ux1 * ni + uy1 * nj;
        u2n = ux2 * ni + uy2 * nj;

        % Касательные компоненты скорости: τ = [-nj, ni]
        ti = -nj; tj = ni;

        u1t = ux1 * ti + uy1 * tj;
        u2t = ux2 * ti + uy2 * tj;

        rho1 = rho(idx1); rho2 = rho(idx2);
        mu1 = mu(idx1);   mu2 = mu(idx2);

        % Напряжения τ_{nτ}
        tau1 = mu1 .* ( ...
            (dudx(idx1)*ni + dudy(idx1)*nj) * ti + ...
            (dvdx(idx1)*ni + dvdy(idx1)*nj) * tj );

        tau2 = mu2 .* ( ...
            (dudx(idx2)*ni + dudy(idx2)*nj) * ti + ...
            (dvdx(idx2)*ni + dvdy(idx2)*nj) * tj );

        % Локальная невязка
        Rloc = zeros(size(idx1));

        % Газ слева
        Rloc(gas_left) = rho1(gas_left) .* u1n(gas_left) .* ...
                         (u1t(gas_left) - u2t(gas_left)) - ...
                         (tau1(gas_left) - tau2(gas_left));

        % Газ справа
        Rloc(gas_right) = rho2(gas_right) .* u2n(gas_right) .* ...
                          (u2t(gas_right) - u1t(gas_right)) - ...
                          (tau2(gas_right) - tau1(gas_right));

        R(idx1) = R(idx1) + Rloc;
    end
end
