function D2 = getSecondScalarDerivativeFull(F, dx, dy, direction)
    % Вычисляет вторую производную по x или y для всей сетки с
    % векторизацией
    % direction = 'x' или 'y'

    [M, N] = size(F);
    D2 = zeros(M, N);

    if direction == "x"
        % Центральная разность для внутренних точек
        D2(:,2:N-1) = (F(:,3:N) - 2*F(:,2:N-1) + F(:,1:N-2)) / dx^2;

        % Односторонняя (вперед/назад) разность на краях
        D2(:,1) = (F(:,3) - 2*F(:,2) + F(:,1)) / dx^2;
        D2(:,N) = (F(:,N) - 2*F(:,N-1) + F(:,N-2)) / dx^2;

    elseif direction == "y"
        D2(2:M-1,:) = (F(3:M,:) - 2*F(2:M-1,:) + F(1:M-2,:)) / dy^2;

        D2(1,:) = (F(3,:) - 2*F(2,:) + F(1,:)) / dy^2;
        D2(M,:) = (F(M,:) - 2*F(M-1,:) + F(M-2,:)) / dy^2;

    else
        error("Direction must be 'x' or 'y'");
    end
end
