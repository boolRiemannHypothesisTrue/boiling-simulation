function val = Qx_dQx_dux(Qx, fields,params,phase,dx,dy)
    
    % Реализует: ∑▒∑▒〖Q_x^(i1,j1)  (∂Q_x^(i1,j1))/(∂u_x^(i,j) )〗
    
    % Извлекаем поля

    u_x = fields.u_x;
    u_y = fields.u_y;
    % Свойства среды
    rho = params.rho_liquid * (1 - phase) + params.rho_vapor * phase;
    mu  = params.mu_liquid  * (1 - phase) + params.mu_vapor  * phase;
   
    [M, N] = size(Qx);
    val = zeros(M, N);

    % Внутренние узлы
    i = 2:M-1; j = 2:N-1;
    
    val(i,j) = Qx(i,j) * rho(i,j) * ( u_x(i+1,j) - u_x(i-1,j) ) / (2*dx) +...
        (4/3) * Qx(i,j) * ( mu(i+1,j) + mu(i,j) ) / (dx^2) + ...
        Qx(i,j) * ( mu(i,j+1) + mu(i,j) ) / (dy^2) + ...
        Qx(i-1,j) * ( u_x(i-1,j) * rho(i-1,j) / (2*dx) - (4/3) * mu(i,j) / dx^2 ) - ...
        Qx(i+1,j) * ( u_x(i+1,j) * rho(i+1,j) / (2*dx) + (4/3) * mu(i+1,j) / dx^2) + ...
        Qx(i,j-1) * ( u_y(i,j-1) / (2*dy) - mu(i,j) / dy^2 ) - ...
        Qx(i,j+1) * ( u_y(i,j+1) / (2*dy) + mu(i,j+1) / dy^2);
end
