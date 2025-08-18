function loss = computeTotalResidualFixedInterface(x, params, dx, dy, X, Y, x_interface, y_interface, phase)
    % Распаковка только полей
    fields = unpackFieldsOnly(x, params.M, params.N);

    % Невязки на основе переданной фазовой маски и интерфейса
    residuals = computeAllResiduals(fields, phase, params, dx, dy, y_interface, x_interface);

    % Суммарная невязка
    loss = 0;
    for k = 1:length(residuals)
        loss = loss + norm(residuals{k}(:));
    end

    % Краевые условия
    M = params.M;
    loss = loss + 1e6 * sum(fields.u_x(1,:).^2);                     % u_x(y=0) = 0
    loss = loss + 1e6 * sum(fields.u_y(1,:).^2);                     % u_y(y=0) = 0
    loss = loss + 1e6 * sum((fields.T(1,:) - params.T_wall).^2);     % T(y=0) = T_wall

    loss = loss + 1e6 * sum(fields.u_x(M,:).^2);                     % u_x(y=ymax) = 0
    loss = loss + 1e6 * sum(fields.u_y(M,:).^2);                     % u_y(y=ymax) = 0
    loss = loss + 1e6 * sum((fields.T(M,:) - params.T_liquid).^2);   % T(y=ymax) = T_liquid

    dPdy = (fields.P(M,:) - fields.P(M-1,:)) / dy;                   % ∂P/∂y(y=ymax) ≈ 0
    loss = loss + 1e6 * sum(dPdy.^2);
end
