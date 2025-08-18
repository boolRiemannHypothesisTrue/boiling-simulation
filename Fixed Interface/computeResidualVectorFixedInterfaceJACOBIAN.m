function [res, J] = computeResidualVectorFixedInterfaceJACOBIAN(x, params, dx, dy, X, Y, x_interface, y_interface, phase)
    % Распаковка только полей
    fields = unpackFieldsOnly(x, params.M, params.N);

    % Основной резидуал
    residuals = computeAllResiduals(fields, phase, params, dx, dy, y_interface, x_interface);
    
    % Склейка невязок
    res = [];
    for k = 1:length(residuals)
        res = [res; residuals{k}(:)];
    end

    % Размерности
    M = params.M;
    N = params.N;
    total_points = M * N;

    % Добавляем граничные условия
    mult = 1e6;
    bc_res = [...
        mult * fields.u_x(1,:)';                         
        mult * fields.u_y(1,:)';                         
        mult * (fields.T(1,:)' - params.T_wall);         
        mult * fields.u_x(M,:)';                         
        mult * fields.u_y(M,:)';                         
        mult * (fields.T(M,:)' - params.T_liquid);       
        mult * ((fields.P(M,:)' - fields.P(M-1,:)') / dy)];
    
    res = [res; bc_res];

    % ====== Численный якобиан (простая версия) ========
    if nargout > 1
        eps_val = 1e-6;
        n = length(x);
        m = length(res);
        J = zeros(m, n);

        for i = 1:n
            x_perturb = x;
            x_perturb(i) = x(i) + eps_val;

            fields_perturb = unpackFieldsOnly(x_perturb, M, N);
            residuals_perturb = computeAllResiduals(fields_perturb, phase, params, dx, dy, y_interface, x_interface);

            res_perturb = [];
            for k = 1:length(residuals_perturb)
                res_perturb = [res_perturb; residuals_perturb{k}(:)];
            end

            % добавим граничные условия
            bc_res_perturb = [...
                mult * fields_perturb.u_x(1,:)';                         
                mult * fields_perturb.u_y(1,:)';                         
                mult * (fields_perturb.T(1,:)' - params.T_wall);         
                mult * fields_perturb.u_x(M,:)';                         
                mult * fields_perturb.u_y(M,:)';                         
                mult * (fields_perturb.T(M,:)' - params.T_liquid);       
                mult * ((fields_perturb.P(M,:)' - fields_perturb.P(M-1,:)') / dy)];

            res_perturb = [res_perturb; bc_res_perturb];

            % дифференцирование по i-й переменной
            J(:, i) = (res_perturb - res) / eps_val;
        end
    end
end
