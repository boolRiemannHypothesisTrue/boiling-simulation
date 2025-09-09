function val = QT_dQT_dT(QT, fields,params,phase,dx,dy)
    
    % Реализует: ∑▒∑▒〖Q_T^(i1,j1)  (∂Q_T^(i1,j1))/(∂T^(i,j) )〗
    
    % Извлекаем поля

    u_x = fields.u_x;
    u_y = fields.u_y;
    % Свойства среды
    rho = params.rho_liquid * (1 - phase) + params.rho_vapor * phase;
    lambda = params.lambda_liquid * (1 - phase) + params.lambda_vapor * phase;
    cp = params.Cp_liquid * (1 - phase) + params.Cp_vapor * phase;
    
    [M, N] = size(QT);
    val = zeros(M, N);

    % Внутренние узлы
    i = 2:M-1; j = 2:N-1;
    
    val(i,j) = QT(i,j) * ( (lambda(i,j+1) + lambda(i,j) ) / dy^2 + ( lambda(i+1,j) + lambda(i,j) ) / dx^2 ) + ...
        QT(i-1,j) * ( rho(i-1,j) * cp(i-1,j) * u_x(i-1,j) / (2*dx) - lambda(i,j) / dx^2 ) - QT(i+1,j) * ( rho(i+1,j) * cp(i+1,j) * u_x(i+1,j) / (2*dx) + lambda(i+1,j) / dx^2 ) + ...
        QT(i,j-1) * ( rho(i,j-1) * cp(i,j-1) * u_y(i,j-1) / (2*dy) - lambda(i,j) / dy^2 ) - QT(i,j+1) * ( rho(i,j+1) * cp(i,j+1) * u_y(i,j+1) / (2*dy) + lambda(i,j+1) / dy^2 );

end
