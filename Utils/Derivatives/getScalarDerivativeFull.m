function D = getScalarDerivativeFull(F, dx, dy, direction)
    % ПРОИЗВОДНАЯ ПО ВСЕЙ СЕТКЕ. ВЕКТОРИЗИРОВАНО, БЕЗ ЦИКЛОВ.
    % Края зануленыибо они не нужны.
    [M, N] = size(F);
    D = zeros(M, N);

    if direction == "x"
        D(:,2:N-1) = (F(:,3:N) - F(:,1:N-2)) / (2*dx);
        % D(:,1)     =(F(:,2)   - F(:,1))     / dx;
        % D(:,N)     =(F(:,N)   - F(:,N-1))   / dx;

    elseif direction == "y"
        D(2:M-1,:) = (F(3:M,:) - F(1:M-2,:)) / (2*dy);
        % D(1,:)     = (F(2,:)   - F(1,:))     / dy;
        % D(M,:)     = (F(M,:)   - F(M-1,:))   / dy;

    else
        error("Direction must be 'x' or 'y'");
    end
end
