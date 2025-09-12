function val = QT_dQT_dux(QT, fields,params,phase,dx,dy)
    % 09.09.2025, written by Mark 
    % Реализует: ∑∑〖Q_T^(i1,j1)  (∂Q_T^(i1,j1))/(∂u_x^(i,j) )〗
    
    % Извлекаем поля
    
    T = fields.T;
    u_x = fields.u_x;
    u_y = fields.u_y;
    
    % Свойства среды

    rho = params.rho_liquid * (1 - phase) + params.rho_vapor * phase;
    cp = params.Cp_liquid * (1 - phase) + params.Cp_vapor * phase;
    mu = params.mu_liquid * (1 - phase) + params.mu_vapor * phase;
    
    [M, N] = size(QT);
    val = zeros(M, N);

    
    % Внутренние узлы
    i = 3:M-2; j = 3:N-2;

    val(i,j) = QT(i,j) * rho(i,j) * cp(i,j) *  ( T(i+1,j) - T(i-1,j) ) / (2*dx) - ...
        QT(i-1,j) * (8/3) * mu(i-1,j) / (2*dx) * ( (u_x(i,j) - u_x(i-2,j) ) / (2*dx) + ...
        ( u_y(i-1,j+1) - u_y(i-1,j-1) ) / (2*dy) ) + ...
        ...
        QT(i+1,j) * (8/3) * mu(i+1,j) / (2*dx) * ( (u_x(i+2,j) - u_x(i,j) ) / (2*dx) + ...
        ( u_y(i+1,j+1) - u_y(i+1,j-1) ) / (2*dy) ) + ...
        ...
        QT(i-1,j) * (16/3) * mu(i-1,j) * ( u_y(i-1,j+1) - u_y(i-1,j-1) ) / (2*dx*2*dy) - ...
        QT(i+1,j) * (16/3) * mu(i+1,j) * ( u_y(i+1,j+1) - u_y(i+1,j-1) ) / (2*dx*2*dy) - ...
        QT(i,j-1) * 2 * mu(i,j-1) / (2*dy) * ( ( u_y(i+1,j-1) - u_y(i-1,j-1) ) / (2*dx) + ( u_x(i,j) - u_x(i,j-2) ) / (2*dy) ) + ...
        QT(i,j-1) * 2 * mu(i,j+1) / (2*dy) * ( ( u_y(i+1,j+1) - u_y(i-1,j+1) ) / (2*dx) + ( u_x(i,j+2) - u_x(i,j) ) / (2*dy) ) ; 

        % горизонтальные линии
    [i1, j1] = ndgrid([2, M-1], 2:N-1);

    % вертикальные линии
    [i2, j2] = ndgrid(2:M-1, [2, N-1]);

    % объединяем
    i = [i1(:); i2(:)];
    j = [j1(:); j2(:)];

    val(i,j) = QT(i,j) * rho(i,j) * cp(i,j) *  ( T(i+1,j) - T(i-1,j) ) / (2*dx) - ...
        QT(i-1,j) * (16/3) * mu(i-1,j) * ( u_y(i-1,j+1) - u_y(i-1,j-1) ) / (2*dx*2*dy) - ...
        QT(i+1,j) * (16/3) * mu(i+1,j) * ( u_y(i+1,j+1) - u_y(i+1,j-1) ) / (2*dx*2*dy);

end
