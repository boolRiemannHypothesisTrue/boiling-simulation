function [X, Y, x_interface, x_nodes, y_nodes, dx, dy] = buildGrid(N, M, x_min, x_max, y_min, y_max)
    % Строит сетку (без интерфейса)

    x_nodes = linspace(x_min, x_max, N+1);
    y_nodes = linspace(y_min, y_max, M+1);

    dx = x_nodes(2) - x_nodes(1);
    dy = y_nodes(2) - y_nodes(1);

    x_centers = x_nodes(1:end-1) + dx/2;
    y_centers = y_nodes(1:end-1) + dy/2;

    [X, Y] = meshgrid(x_centers, y_centers);

    x_interface = x_centers';  % это x, по которому считается y_interface
end
