function res = computeResidualVectorFixedInterface(x, params, dx, dy, X, Y, x_interface, y_interface, phase)
    % Распаковка только полей (без c)
    fields = unpackFieldsOnly(x, params.M, params.N);

    % Вычисление невязок по фиксированной фазе и интерфейсу
    residuals = computeAllResiduals(fields, phase, params, dx, dy, y_interface, x_interface);
    
    % Объединяем все невязки
    res = [];
    for k = 1:length(residuals)
        res = [res; residuals{k}(:)];
    end

    % Добавляем штраф за граничные условия
    mult = 1e6;
    M = params.M;
    res = [res;
           mult * fields.u_x(1,:)';                         
           mult * fields.u_y(1,:)';                         
           mult * (fields.T(1,:)' - params.T_wall);         
           mult * fields.u_x(M,:)';                         
           mult * fields.u_y(M,:)';                         
           mult * (fields.T(M,:)' - params.T_liquid);       
           mult * ((fields.P(M,:)' - fields.P(M-1,:)') / dy)];
end
