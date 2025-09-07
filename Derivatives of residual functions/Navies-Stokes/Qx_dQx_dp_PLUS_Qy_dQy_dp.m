function val = Qx_dQx_dp_PLUS_Qy_dQy_dp(Qx, Qy, dx, dy)
    %
    % Реализует: ∑∑(Q_x^(i1,j1)  (∂Q_x^(i1,j1))/(∂P^(i,j) )+Q_y^(i1,j1)  (∂Q_y^(i1,j1))/(∂P^(i,j) )) 
    % и да, я вставил это из ворда напрямую,лень перепечатывать 

    [M, N] = size(Qx);
    val = zeros(M, N);

    % Внутренние узлы
    i = 2:M-1; j = 2:N-1;
    val(i,j) = ( Qx(i-1,j) - Qx(i+1,j) ) / (2*dx) ...
              + ( Qy(i,  j-1) - Qy(i,  j+1) ) / (2*dy);
end
