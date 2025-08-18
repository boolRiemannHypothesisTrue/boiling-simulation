function D = getProductDerivative(mu, f, dx, dy, direction)
    % Вычисляет ∂/∂x(mu*f) или ∂/∂y(mu*f), с учетом односторонних разностей на границе

    prod = mu .* f;
    D = getScalarDerivativeFull(prod, dx, dy, direction);
end
