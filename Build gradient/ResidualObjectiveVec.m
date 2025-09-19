function [J, gradVec] = ResidualObjectiveVec(qVec, fields, params, phase, dx, dy, sz, X, Y)
    % --- 1. Распаковать параметры в поля ---
    [ux, uy, T, p] = unpackParams(qVec, sz);
    fields.u_x = ux;
    fields.u_y = uy;
    fields.T  = T;
    fields.P  = p;

    % --- 2. Невязки по PDE ---
    Qu = computeContinuityResidual(fields, params, phase, dx, dy);
    QT = computeEnergyResidual(fields, params, phase, dx, dy);
    Qx = computeNSResidualX_full(fields, params, phase, dx, dy);
    Qy = computeNSResidualY_full(fields, params, phase, dx, dy);

    % --- 3. Невязки по краевым условиям ---
    Qbux1 = Qb_ux_1(fields);
    Qbux2 = Qb_ux_2(fields);
    Qbuy1 = Qb_uy_1(fields);
    Qbuy2 = Qb_uy_2(fields);
    QbT1  = Qb_T_1(fields,params);
    QbT2  = Qb_T_2(fields,params);
    Qbp   = Qb_p(fields);

    % --- 4. Нормировка невязок ---
    normFactor = max(1, max([max(abs(Qu(:))), max(abs(QT(:))), max(abs(Qx(:))), max(abs(Qy(:))), ...
                             max(abs(Qbux1(:))), max(abs(Qbux2(:))), max(abs(Qbuy1(:))), max(abs(Qbuy2(:))), ...
                             max(abs(QbT1(:))), max(abs(QbT2(:))), max(abs(Qbp(:)))]));
    Qu = Qu / normFactor;
    QT = QT / normFactor;
    Qx = Qx / normFactor;
    Qy = Qy / normFactor;
    Qbux1 = Qbux1 / normFactor;
    Qbux2 = Qbux2 / normFactor;
    Qbuy1 = Qbuy1 / normFactor;
    Qbuy2 = Qbuy2 / normFactor;
    QbT1  = QbT1  / normFactor;
    QbT2  = QbT2  / normFactor;
    Qbp   = Qbp   / normFactor;

    % --- 5. Функция цели ---
    J = 0.5 * ( sum(Qu(:).^2) + sum(QT(:).^2) + sum(Qx(:).^2) + sum(Qy(:).^2) + ...
                sum(Qbux1(:).^2) + sum(Qbux2(:).^2) + sum(Qbuy1(:).^2) + sum(Qbuy2(:).^2) + ...
                sum(QbT1(:).^2) + sum(QbT2(:).^2) + sum(Qbp(:).^2) );

    % --- Проверка на NaN/Inf ---
    if ~isfinite(J)
        error('ResidualObjectiveVec: J is not finite! Check initial conditions and PDE residuals.');
    end

    % --- 6. Градиенты ---
    gradUx = grad_ux(Qu, QT, Qx, Qy, Qbux1, Qbux2, fields, params, phase, dx, dy);
    gradUy = grad_uy(Qu, QT, Qx, Qy, Qbuy1, Qbuy2, fields, params, phase, dx, dy);
    gradT  = grad_T(QT, QbT1, QbT2, fields, params, phase, dx, dy);
    gradP  = grad_p(Qx, Qy, Qbp, dx, dy);

    % --- 7. Упаковать градиенты в вектор ---
    gradVec = packParams(gradUx, gradUy, gradT, gradP);

    % --- 8. Опциональный вывод для отладки ---
    % disp(['J = ', num2str(J), ', ||grad|| = ', num2str(norm(gradVec))]);
end
