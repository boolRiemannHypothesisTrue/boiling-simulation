function val = Qx_dQx_duy(Qx, fields,params,phase,dx,dy)
    
    % Реализует: ∑▒∑▒〖Q_x^(i1,j1)  (∂Q_x^(i1,j1))/(∂u_y^(i,j) )〗
    
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
    
    val(i,j) = Qx(i,j) * rho(i,j) * ( u_x(i,j+1) - u_x(i,j-1) ) / (2*dy) + ...
        Qx(i-1,j) / (2*dx) * ( (2/3) * ( mu(i,j-1) - mu(i,j+1) ) / (2*dy) + ...
        ( mu(i-1,j-1) - mu(i-1,j+1) / (2*dy) ) ) - ...
        Qx(i+1,j) / (2*dx) * ( (2/3) * (mu(i,j-1) - mu(i,j+1) ) / (2*dy) + ...
        ( mu(i+1,j+1) - mu(i+1,j-1) ) / (2*dy));
end
