function R = computeResidual_9(fields, params, phase, dx, dy)
    [M, N] = size(phase);

    % Физические параметры по фазе
    rho = params.rho_liquid * (phase == 0) + params.rho_vapor * (phase == 1);
    mu  = params.mu_liquid  * (phase == 0) + params.mu_vapor  * (phase == 1);
    h   = params.h_liquid   * (phase == 0) + params.h_vapor   * (phase == 1);

    % Поля
    ux = fields.u_x;
    uy = fields.u_y;
    [dudx, dudy] = gradient(ux, dx, dy);
    [dvdx, dvdy] = gradient(uy, dx, dy);

    % Инициализация
    R = zeros(M, N);

    % Градиенты температуры для теплового потока
    q_jump = computeHeatFluxJump(fields, params, phase, dx, dy);

    % Соседние направления и нормали
    shifts  = {[-1, 0], [0, -1], [1, 0], [0, 1]};
    normals = {[0, -1], [-1, 0], [0, 1], [1, 0]};

    [I, J] = ndgrid(1:M, 1:N);

    for d = 1:4
        shift  = shifts{d};
        normal = normals{d};
        ni = normal(1); nj = normal(2);

        I2 = I + shift(1);
        J2 = J + shift(2);

        inside = I2 >= 1 & I2 <= M & J2 >= 1 & J2 <= N;
        if ~any(inside(:)), continue; end

        idx1 = find(inside);
        idx2 = sub2ind([M, N], I2(idx1), J2(idx1));

        p1 = phase(idx1);
        p2 = phase(idx2);

        is_interface = p1 ~= p2;
        is_gas1 = p1 == 1;
        is_gas2 = p2 == 1;

        gas_left  = is_interface & is_gas1 & ~is_gas2;
        gas_right = is_interface & ~is_gas1 & is_gas2;

        % Поля
        rho1 = rho(idx1); rho2 = rho(idx2);
        h1   = h(idx1);   h2   = h(idx2);
        ux1  = ux(idx1);  ux2  = ux(idx2);
        uy1  = uy(idx1);  uy2  = uy(idx2);
        mu1  = mu(idx1);  mu2  = mu(idx2);

        dudx1 = dudx(idx1); dudy1 = dudy(idx1);
        dvdx1 = dvdx(idx1); dvdy1 = dvdy(idx1);

        dudx2 = dudx(idx2); dudy2 = dudy(idx2);
        dvdx2 = dvdx(idx2); dvdy2 = dvdy(idx2);

        % Векторные величины
        u1n = ux1 * ni + uy1 * nj;
        u2n = ux2 * ni + uy2 * nj;
        u1sq = ux1.^2 + uy1.^2;
        u2sq = ux2.^2 + uy2.^2;

        % Проекция u_k τ_nk
        tau1 = mu1 .* viscous_proj(dudx1, dudy1, dvdx1, dvdy1, ni, nj);
        tau2 = mu2 .* viscous_proj(dudx2, dudy2, dvdx2, dvdy2, ni, nj);
        utau1 = ux1 .* ni + uy1 .* nj;  % u_k * n_k
        utau2 = ux2 .* ni + uy2 .* nj;

        % Невязка
        Rloc = zeros(size(idx1));

        % Газ слева (ячейка idx1 — пар)
        Rloc(gas_left) = ...
            rho1(gas_left) .* u1n(gas_left) .* ...
            (h1(gas_left) - h2(gas_left) + 0.5 * (u1sq(gas_left) - u2sq(gas_left))) ...
            - utau1(gas_left) .* tau1(gas_left) + utau2(gas_left) .* tau2(gas_left);

        % Газ справа (ячейка idx2 — пар)
        Rloc(gas_right) = ...
            rho2(gas_right) .* u2n(gas_right) .* ...
            (h2(gas_right) - h1(gas_right) + 0.5 * (u2sq(gas_right) - u1sq(gas_right))) ...
            - utau2(gas_right) .* tau2(gas_right) + utau1(gas_right) .* tau1(gas_right);

        % Добавим в итоговую невязку
        R(idx1) = R(idx1) + Rloc;
    end

    % Добавим скачок теплового потока
    R = R + q_jump;
end

function tau = viscous_proj(dudx, dudy, dvdx, dvdy, ni, nj)
    % τ_nk проецированное: u_k * τ_nk = u_i * τ_ni
    tau = ni.^2 .* dudx + nj.^2 .* dvdy + ni .* nj .* (dudy + dvdx);
end
