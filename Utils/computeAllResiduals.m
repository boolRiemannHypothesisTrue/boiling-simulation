function residuals = computeAllResiduals(fields, phase, params, dx, dy, y_interface, x_interface)
    % Вызов 9 функций невязок
    residuals = cell(8,1);
    H = computeCurvatureFromInterface(x_interface,y_interface); % кривизна под 7
    
    residuals{1} = computeContinuityResidual(fields,params, phase, dx, dy);
    residuals{2} = computeNSResidualX_full(fields ,params,phase, dx, dy);
    residuals{3} = computeNSResidualY_full(fields, params,phase, dx, dy);
    residuals{4} = computeEnergyResidual(fields, params,phase, dx, dy);
    residuals{5} = computeResidual_6(fields, params, phase, dx, dy);
    residuals{6} = computeResidual_7(fields, params, phase, H, dx, dy);
    residuals{7} = computeResidual_8(fields, params, phase, dx, dy);
    residuals{8} = computeResidual_9(fields, params, phase, dx, dy);

    % Каждая computeResidual - функция, возвращающая норму невязок,
    % в computeTotalResidual они будут суммироваться.
end
