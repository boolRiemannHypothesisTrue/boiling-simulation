clear; clc; close all;
cd("C:\Users\MSI\Desktop\sim") % project folder
addpath(genpath(cd)); % path to folder
% Загрузить параметры

params = config();

% Построение сетки и матрицы фаз
tic

% Постоянная часть: сетка
[X, Y, x_interface, x_nodes, y_nodes, dx, dy] = ...
    buildGrid(params.N, params.M, ...
              params.x_min, params.x_max, ...
              params.y_min, params.y_max);

% Перестраивается на каждой итерации (при новом c)
[y_interface, phase] = buildInterfaceAndPhase(X, Y, x_interface, params.c_init, ...
                                              params.y_min, params.y_max);

toc


% Инициализация полей
fields = initializeFields(X, Y, phase, params);

% Подготовка границ для GA с учётом фаз
total_points = params.M * params.N;
phase_vec = reshape(phase, total_points, 1);

% Ограничения скорости
u_min_liquid = -1;
u_max_liquid = 1;
u_min_vapor  = -10;
u_max_vapor  = 10;

% Ограничения температуры
T_wall = params.T_wall;
T_sat  = params.T_sat;       % Добавь T_sat в params
T_liquid = params.T_liquid;

% Ограничения давления
P0 = 1e5;
P_min = 0.9 * P0;
P_max = 1.1 * P0;

lb_ux = zeros(total_points,1);
ub_ux = zeros(total_points,1);
lb_uy = zeros(total_points,1);
ub_uy = zeros(total_points,1);
lb_T  = zeros(total_points,1);
ub_T  = zeros(total_points,1);
lb_P  = P_min * ones(total_points,1);
ub_P  = P_max * ones(total_points,1);

for i = 1:total_points
    if phase_vec(i) == 1  % жидкость
        lb_ux(i) = u_min_liquid;
        ub_ux(i) = u_max_liquid;
        lb_uy(i) = u_min_liquid;
        ub_uy(i) = u_max_liquid;

        lb_T(i) = params.T_liquid;  % 300 К — нижняя граница для жидкости
        ub_T(i) = params.T_sat;     % 373 К — верхняя граница для жидкости
    else  % пар
        lb_ux(i) = u_min_vapor;
        ub_ux(i) = u_max_vapor;
        lb_uy(i) = u_min_vapor;
        ub_uy(i) = u_max_vapor;

        lb_T(i) = params.T_sat;     % 373 К — нижняя граница для пара
        ub_T(i) = params.T_wall;    % 400 К — верхняя граница для пара
    end
end


% Распаковка и установка границ в общий вектор для GA
if exist('x_opt', 'var')
    x0 = x_opt;
else
    %x0 = packFieldsAndInterface(fields, params.M, params.N, length(params.c_init), params.c_init);
    x0 = packFieldsOnly(fields, params.M, params.N);
end

nvars = length(x0);

% Инициализация границ
lb = -Inf(nvars, 1);
ub = Inf(nvars, 1);

% Присваиваем границы по переменным
lb(1:total_points)                   = lb_ux;
ub(1:total_points)                   = ub_ux;

lb(total_points+1:2*total_points)   = lb_uy;
ub(total_points+1:2*total_points)   = ub_uy;

lb(2*total_points+1:3*total_points) = lb_T;
ub(2*total_points+1:3*total_points) = ub_T;

lb(3*total_points+1:4*total_points) = lb_P;
ub(3*total_points+1:4*total_points) = ub_P;


% Целевая функция
fun = @(x) computeTotalResidualFixedInterface(x,params,dx,dy,X,Y,x_interface,y_interface,phase);

% Параметры оптимизации
target_loss = 1e5;
max_restarts = 20;
restart_count = 0;
fval = Inf;

% Главный цикл
while fval > target_loss && restart_count < max_restarts
    fprintf('\n🌀 GA Рестарт #%d\n', restart_count + 1);

    num_individuals = 30;
    if exist('x_opt', 'var')
        init_pop = [x_opt(:)'; rand(num_individuals - 1, nvars) .* (ub' - lb') + lb'];
    else
        init_pop = rand(num_individuals, nvars) .* (ub' - lb') + lb';
    end

    options = optimoptions('ga', ...
        'Display', 'iter', ...
        'PopulationSize', num_individuals, ...
        'MaxGenerations', 1000, ...
        'FunctionTolerance', 1e-6, ...
        'UseParallel', true, ...
        'InitialPopulationMatrix', init_pop);

    [x_opt, fval, exitflag, output] = ga(fun, nvars, [], [], [], [], lb, ub, [], options);

    fprintf('✅ GA завершил рестарт #%d: fval = %.4e\n', restart_count + 1, fval);
    save(sprintf('opt_restart_%02d.mat', restart_count + 1), 'x_opt', 'fval', 'exitflag', 'output');

    restart_count = restart_count + 1;
end

%

visualizeGridAndInterface(x_nodes, y_nodes, dx, dy, phase, params.c_init);

field_opt = unpackFieldsOnly(x_opt,params.M,params.N);

plotInitialFields(X,Y,fields_opt)
