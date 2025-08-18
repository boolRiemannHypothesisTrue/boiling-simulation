function [R, J] = computeContinuityResidual_JACOBIAN(fields, params, phase, dx, dy)
    % Вычисление невязки уравнения непрерывности div(rho * u) = 0
    % и аналитического якобиана по u_x и u_y
    %
    % Вход:
    %   fields - структура с полями u_x, u_y размером MxN
    %   params - параметры с rho_liquid и rho_vapor
    %   phase  - матрица фазы (M x N)
    %   dx, dy - шаги сетки
    %
    % Выход:
    %   R - вектор невязки (M*N x 1)
    %   J - якобиан (M*N x 2*M*N), производные R по u_x и u_y
        % Вычисление невязки уравнения непрерывности div(rho * u) = 0
    % Векторизованная версия

    u_x = fields.u_x;
    u_y = fields.u_y;

    % Плотность по фазе
    rho = params.rho_liquid * (1 - phase) + params.rho_vapor * phase;

    % Потоки импульса
    rhou_x = rho .* u_x;
    rhou_y = rho .* u_y;

    % Производные
    drhou_dx = getScalarDerivativeFull(rhou_x, dx, dy, 'x');
    drhou_dy = getScalarDerivativeFull(rhou_y, dx, dy, 'y');

    % Суммарная дивергенция
    R = drhou_dx + drhou_dy;

  
    [M, N] = size(u_x);


    % Векторизация
    rho_vec = rho(:);         % (M*N x 1)


    % Построение матриц производных
    [Dx, Dy] = buildCentralDifferenceMatrices(M, N, dx, dy);



    % Якобиан по u_x:
    % dR/du_x = Dx * diag(rho)
    J_ux = Dx * spdiags(rho_vec, 0, M*N, M*N);

    % Якобиан по u_y:
    % dR/du_y = Dy * diag(rho)
    J_uy = Dy * spdiags(rho_vec, 0, M*N, M*N);

    % Итоговый якобиан: столбцы - u_x и u_y соответственно
    % Размер: (M*N) x (2*M*N)
    J = [J_ux, J_uy];
end


% Вспомогательная функция для построения матриц производных (если ещё не сделана):
function [Dx, Dy] = buildCentralDifferenceMatrices(M, N, dx, dy)
    eN = ones(N,1);
    Dx1D = spdiags([-eN, eN], [-1,1], N, N)/(2*dx);
    Dx1D(1,1:2) = [-1,1]/dx;
    Dx1D(N,N-1:N) = [-1,1]/dx;

    eM = ones(M,1);
    Dy1D = spdiags([-eM, eM], [-1,1], M, M)/(2*dy);
    Dy1D(1,1:2) = [-1,1]/dy;
    Dy1D(M,M-1:M) = [-1,1]/dy;

    Ix = speye(M);
    Iy = speye(N);

    Dx = kron(Ix, Dx1D);
    Dy = kron(Dy1D, Iy);
end
