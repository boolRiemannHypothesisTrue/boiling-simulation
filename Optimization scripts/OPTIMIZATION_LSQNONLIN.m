clear; clc; close all;
cd("C:\Users\MSI\Desktop\sim") % project folder
addpath(genpath(cd)); % path to folder

% Загрузить параметры
params = config();

% Построение сетки и матрицы фаз

[X, Y, x_interface, x_nodes, y_nodes, dx, dy] = ...
    buildGrid(params.N, params.M, ...
              params.x_min, params.x_max, ...
              params.y_min, params.y_max);

[y_interface, phase] = buildInterfaceAndPhase(X, Y, x_interface, params.c_init, ...
                                              params.y_min, params.y_max);

% Инициализация полей
fields = initializeFields(X, Y, phase, params);

% Подготовка границ для lsqnonlin с учётом фаз
total_points = params.M * params.N;
phase_vec = reshape(phase, total_points, 1);

% Ограничения скорости
u_min_liquid = -1;
u_max_liquid = 1;
u_min_vapor  = -10;
u_max_vapor  = 10;

% Ограничения температуры
T_wall = params.T_wall;
T_sat  = params.T_sat;
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

        lb_T(i) = T_liquid;  % 300 К — нижняя граница для жидкости
        ub_T(i) = T_sat;     % 373 К — верхняя граница для жидкости
    else  % пар
        lb_ux(i) = u_min_vapor;
        ub_ux(i) = u_max_vapor;
        lb_uy(i) = u_min_vapor;
        ub_uy(i) = u_max_vapor;

        lb_T(i) = T_sat;     % 373 К — нижняя граница для пара
        ub_T(i) = T_wall;    % 400 К — верхняя граница для пара
    end
end

% Начальное приближение
if exist('x_opt', 'var')
    x0 = x_opt;
else
    x0 = packFieldsOnly(fields, params.M, params.N); % или packFieldsAndInterface
end

nvars = length(x0);

% Инициализация границ для всех переменных
lb = -Inf(nvars, 1);
ub = Inf(nvars, 1);

% Присваиваем границы по переменным: ux, uy, T, P — по блокам
lb(1:total_points)                   = lb_ux;
ub(1:total_points)                   = ub_ux;

lb(total_points+1:2*total_points)   = lb_uy;
ub(total_points+1:2*total_points)   = ub_uy;

lb(2*total_points+1:3*total_points) = lb_T;
ub(2*total_points+1:3*total_points) = ub_T;

lb(3*total_points+1:4*total_points) = lb_P;
ub(3*total_points+1:4*total_points) = ub_P;

% Целевая функция — вектор невязки
fun = @(x) computeResidualVectorFixedInterfaceJACOBIAN(x, params, dx, dy, X, Y, x_interface,y_interface,phase);

% Параметры останова
target_loss = 1e5;
max_restarts = 20;
restart_count = 0;
fval = Inf;

typicalX = ones(nvars, 1);

% Примерный масштаб по блокам:
total_points = params.M * params.N;
typicalX(1:total_points)                     = 1;      % u_x
typicalX(total_points+1:2*total_points)     = 1;      % u_y
typicalX(2*total_points+1:3*total_points)   = 300;    % T — порядка 300-400 K
typicalX(3*total_points+1:4*total_points)   = 1e5;    % P — порядка 1e5


while fval > target_loss && restart_count < max_restarts
    fprintf('\n🔁 LSQNONLIN Рестарт #%d\n', restart_count + 1);

    options = optimoptions('lsqnonlin', ...
        'Display', 'iter', ...
        'UseParallel', true, ...
        'TypicalX', typicalX, ...
        'FunctionTolerance', 1e-6, ...
        'StepTolerance', 1e-10, ...
        'Algorithm', 'trust-region-reflective',...
        'TolX',1e-3,...
        'FiniteDifferenceType', 'central',...
           'ScaleProblem', 'none',...
           'Jacobian','on');

    [x_opt, resvec, ~, exitflag, output] = lsqnonlin(fun, x0, lb, ub, options);
    fval = sum(resvec.^2);

    % 🔄 Сохраняем в один и тот же файл
    save('opt_lsqnonlin_latest.mat');

    fprintf('✅ Рестарт #%d завершён: fval = %.4e\n', restart_count + 1, fval);

    x0 = x_opt;
    restart_count = restart_count + 1;
end

if fval <= target_loss
    fprintf('\n🎯 Сходимость достигнута: fval = %.4e\n', fval);
else
    fprintf('\n⚠️ Сходимость НЕ достигнута после %d рестартов. Итоговая fval = %.4e\n', restart_count, fval);
end

%%
% Визуализация результатов
visualizeGridAndInterface(x_nodes, y_nodes, dx, dy, phase, params.c_init);

fields_opt = unpackFieldsOnly(x_opt, params.M, params.N);
plotInitialFields(X, Y, fields_opt);
