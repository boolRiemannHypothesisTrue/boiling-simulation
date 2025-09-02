function val = dQxdu_x(Qx,params,phase,dx)
    % Частные производная невязки Х-проекции уравнений Навье Стокса Qx по u_x  
    [M, N] = size(Qx);
    val = zeros(M, N);
    rho = params.rho_liquid * (1 - phase) + params.rho_vapor * phase;
    mu = params.mu_liquid * (1 - phase) + params.mu_vapor * phase;
    % Внутренние узлы
    i = 2:M-1; j = 2:N-1;
    val(i,j) = Qx(i,j) * rho(i,j) * ( u_x(i+1,j) - u_x(i-1,j) ) / ( 2* dx) +...
        (4/3) * Qx(i,j) *  ( mu(i+1,j) + mu(i,j) ) / dx^2 + ...
        Qx(i,j) * ( mu(i,j+1) + mu(i,j) ) / dy^2 + ...
        Qx(i-1,j) *
    % TODO: FINISH
end
