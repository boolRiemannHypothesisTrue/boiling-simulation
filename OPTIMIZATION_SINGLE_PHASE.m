clear; clc; close all;
%cd("C:\Users\MSI\Desktop\sim") % project folder
cd("C:\Users\MSI\Documents\GitHub\boiling-simulation");
addpath(genpath(cd)); % path to folder
% Загрузить параметры

params = config();

% Построение сетки и матрицы фаз
tic
% 
% [X, Y, phase, y_interface,x_interface, x_nodes, y_nodes, dx, dy] = ...
%     buildGridAndInterface(params.N, params.M, ...
%                           params.x_min, params.x_max, ...
%                           params.y_min, params.y_max, ...
%                           params.c_init);


% Постоянная часть: сетка
[X, Y, x_interface, x_nodes, y_nodes, dx, dy] = ...
    buildGrid(params.N, params.M, ...
              params.x_min, params.x_max, ...
              params.y_min, params.y_max);


[y_interface, phase] = buildInterfaceAndPhase(X, Y, x_interface, params.c_init, ...
                                              params.y_min, params.y_max);
% one phase
[M,N] = size(phase);
phase = ones(M,N);
toc


% Инициализация полей

fields = initializeFields(X, Y, phase, params);
% размеры сетки
[Nx, Ny] = size(fields.u_x);

% упаковать начальные параметры
q0 = packParams(fields.u_x, fields.u_y, fields.T, fields.P);

% параметры Adam
stepSize = 1e-4;  % можно менять по необходимости
beta1 = 0.9;
beta2 = 0.999;
epsilon = 1e-8;
nEpochSize = 50;
options = optimset('GradObj','on','Display','iter');

% функция цели
fun = @(q) ResidualObjectiveVec(q, fields, params, phase, dx, dy, [Nx, Ny]);

% запуск оптимизации
[qOptVec, fval, exitflag, output] = fmin_adam(fun, q0, stepSize, beta1, beta2, epsilon, nEpochSize, options);

% распаковка оптимальных параметров
[uxOpt, uyOpt, TOpt, pOpt] = unpackParams(qOptVec, [Nx, Ny]);
