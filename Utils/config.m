function params = config()
    % Параметры сетки
    params.N = 500;          % число ячеек по x
    params.M = 500;          % число ячеек по y

    params.x_min = 0;
    params.x_max = 5e-2;
    params.y_min = 0;
    params.y_max = 5e-2;
    
    % Параметры межфазной границы (полином 15-й степени)
    params.poly_degree = 15;   
    params.c_init = zeros(params.poly_degree+1,1);
    params.c_init(1) = 0.03;   % начальный интерфейс - const y = 0.03
    params.c_init(2) = 0.3;
   
    % Физические параметры (пример)
    params.rho_liquid = 1000;     % плотность жидкости, кг/м^3
    params.rho_vapor = 0.6;       % плотность пара, кг/м^3
    params.mu_liquid = 1e-3;      % вязкость жидкости, Па·с
    params.mu_vapor = 1e-5;       % вязкость пара, Па·с
    params.lambda_liquid = 0.6;   % теплопроводность жидкости, Вт/(м·К)
    params.lambda_vapor = 0.02;   % теплопроводность пара, Вт/(м·К)
    params.Cp_liquid = 4200;      % теплоёмкость жидкости, Дж/(кг·К)
    params.Cp_vapor = 2000;       % теплоёмкость пара, Дж/(кг·К)
    params.sigma = 0.072;  % Поверхностное натяжение (Н/м)
    params.h_liquid = 4.19e5;   % [Дж/кг] — насыщенная вода при 100 °C
    params.h_vapor  = 2.676e6;  % [Дж/кг] — насыщенный пар при 100 °C

    % Граничные условия
    params.T_wall = 400;    % Температура стенки, К
    params.T_liquid = 300;  % Температура жидкости, К
    params.T_sat = 373; % Насыщение
    % Параметры оптимизации
    params.optim_maxIter = 500;   % макс. итераций оптимизации
    params.optim_tol = 1e-6;      % критерий сходимости

    % Другие настройки
    params.gravity = 9.81;        % ускорение свободного падения, м/с^2
end
