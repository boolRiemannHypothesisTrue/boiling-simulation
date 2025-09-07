function val = Qy_dQy_dux(Qy, fields,params,phase,dx,dy)
    
    % Реализует: ∑▒∑▒〖Q_y^(i1,j1)  (∂Q_y^(i1,j1))/(∂u_x^(i,j) )〗
    
    % Извлекаем поля

    u_x = fields.u_x;
    u_y = fields.u_y;

    % Свойства среды
    rho = params.rho_liquid * (1 - phase) + params.rho_vapor * phase;
    mu  = params.mu_liquid  * (1 - phase) + params.mu_vapor  * phase;
   
    [M, N] = size(Qy);
    val = zeros(M, N);

    % Внутренние узлы
    i = 2:M-1; j = 2:N-1;
    
    val(i,j) = Qy(i,j) * rho(i,j) * ( u_y(i+1,j) - u_x(i-1,j) ) / (2*dx) + ...
        Qy(i-1,j) / (2*dx) * ( ( mu(i,j-1) - mu(i,j+1) ) / (2*dy) - ...
        (2/3) * ( mu(i-1,j+1) - mu(i-1,j-1) / (2*dy) ) ) - ...
        Qy(i+1,j) / (2*dx) * ( (mu(i,j-1) - mu(i,j+1) ) / (2*dy) - ...
        (2/3) * ( mu(i+1,j+1) - mu(i+1,j-1) ) / (2*dy));
end
