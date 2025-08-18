function res = computeResidualVector(x, params, dx, dy, X, Y, x_interface)
    [fields, c] = unpackFieldsAndInterface(x, params.M, params.N, length(params.c_init));
    [y_interface, phase] = buildInterfaceAndPhase(X, Y, x_interface, c, params.y_min, params.y_max);

    residuals = computeAllResiduals(fields, phase, params, dx, dy, y_interface, x_interface);
    
    % res = zeros(length(residuals),1);
    res = [];
    for k = 1:length(residuals)
        res = [res; residuals{k}(:)];  % Невязки уравнений
    end

    % Граничные условия
    mult = 1e6; % множитель для штрафа за отклонение от краевых условий
    M = params.M;
    res = [res;
           mult * fields.u_x(1,:)';                         % u_x(y=0) = 0
           mult * fields.u_y(1,:)';                         % u_y(y=0) = 0
           mult * (fields.T(1,:)' - params.T_wall);         % T(y=0) = T_wall
           mult * fields.u_x(M,:)';                         % u_x(y=ymax) = 0
           mult * fields.u_y(M,:)';                         % u_y(y=ymax) = 0
           mult * (fields.T(M,:)' - params.T_liquid);       % T(y=ymax) = T_liquid
           mult * ((fields.P(M,:)' - fields.P(M-1,:)') / dy)];  % ∂P/∂y ≈ 0
end
