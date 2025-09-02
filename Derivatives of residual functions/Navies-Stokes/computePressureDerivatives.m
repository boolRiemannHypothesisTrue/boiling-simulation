function dRdP = computePressureDerivatives(Qx, Qy, dx, dy)
    % Частные производные невязок по давлению P (структурированная сетка, w=1)
    % Реализует:
    % dR/dP(i,j) = (Qx(i-1,j)-Qx(i+1,j))/(2*dx) + (Qy(i,j-1)-Qy(i,j+1))/(2*dy)

    [M, N] = size(Qx);
    dRdP = zeros(M, N);

    % Внутренние узлы
    i = 2:M-1; j = 2:N-1;
    dRdP(i,j) = ( Qx(i-1,j) - Qx(i+1,j) ) / (2 * dx) ...
              + ( Qy(i,  j-1) - Qy(i,  j+1) ) /(2 * dy);
end
