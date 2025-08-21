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

% Перестраивается на каждой итерации (при новом c)
[y_interface, phase] = buildInterfaceAndPhase(X, Y, x_interface, params.c_init, ...
                                              params.y_min, params.y_max);

toc


% Инициализация полей

fields = initializeFields(X, Y, phase, params);

%% Визуализация сетки и интерфейса

visualizeGridAndInterface(x_nodes, y_nodes, dx, dy, phase, params.c_init);

%% Визуализация начальных полей

plotInitialFields(X, Y, fields);

%% Уравнение непрерывности (3)
tic

R_continuity = computeContinuityResidual(fields, params, phase, dx, dy);

RR = norm(R_continuity);   % L2-норма (корень суммы квадратов)
fprintf("L2 norm of R = %e\n",RR)

toc

%% Визуализация невязки уравнения непрерывности

figure;

pcolor(X, Y, R_continuity);

shading interp;
colorbar;
colormap jet;
title('dR/du_y');
xlabel('x (м)');
ylabel('y (м)');
axis equal tight;

%% Навье стокс OX (4.1)

tic

R_ns_x = computeNSResidualX_full(fields, params, phase, dx, dy);

norm( R_ns_x ) % L2-норма (корень суммы квадратов)

toc 

%% Визуализация невязки уравнения Навье-Стокса (ОХ)

figure;
pcolor(X, Y, R_ns_x);
shading interp;
colorbar;
colormap jet;
title('Невязка уравнения Навье-Стокса (ОХ)');
xlabel('x (м)');
ylabel('y (м)');
axis equal tight

%% Навье стокс OY (4.2)
tic

R_ns_y = computeNSResidualY_full(fields, params, phase, dx, dy);

norm( R_ns_y )  % L2-норма (корень суммы квадратов)

toc

%% Визуализация невязки уравнения Навье-Стокса (ОУ)

figure;
pcolor(X, Y, R_ns_y);
shading interp;
colorbar;
colormap jet;
title('Невязка уравнения Навье-Стокса (ОY)');
xlabel('x (м)');
ylabel('y (м)');
axis equal tight

%% Уравнение энергии (5)
tic 

R_energy = computeEnergyResidual(fields,params,phase,dx,dy);

norm(R_energy) % L2-норма (корень суммы квадратов)

toc

%% Визуализация невязки уравнения энергии

figure;
pcolor(X, Y, R_energy);
shading interp;
colorbar;
colormap jet;
title('Невязка уравнения Энергии');
xlabel('x (м)');
ylabel('y (м)');
axis equal tight

%% Уравнение совместности 6 

tic

R = computeResidual_6(fields, params, phase, dx, dy);

toc

norm(R)

%% Уравнение 6.Визуализация

figure
pcolor(X, Y, R);
shading interp;
colorbar;
colormap jet;
title('Residual 6');
xlabel('x (м)');
ylabel('y (м)');
axis equal tight

%% Уравнение совместности 7

tic

H = computeCurvatureFromInterface(x_interface,y_interface);

R = computeResidual_7(fields, params, phase, H, dx, dy);

norm(R)

toc

%% Кривизна H по координате X границы раздела фаз (к уравнению совместности 7)

plot(x_interface,H)

%% Уравнение 7.Визуализация

figure
pcolor(X, Y, R);
shading interp;
colorbar;
colormap jet;
title('Residual 7');
xlabel('x (м)');
ylabel('y (м)');
axis equal tight

%% Уравнение совместности 8

R = computeResidual_8(fields, params, phase, dx,dy);

norm(R)

%% Уравнение 8.Визуализация

figure
pcolor(X, Y, R);
shading interp;
colorbar;
colormap jet;
title('Residual 8');
xlabel('x (м)');
ylabel('y (м)');
axis equal tight

%% q_v - q_l на межфазной границе (к уравнению совместности 9)

dq = computeHeatFluxJump(fields, params, phase, dx, dy);

figure
pcolor(X, Y,dq);
shading interp;
colorbar;
colormap jet;
xlabel('x (м)');
ylabel('y (м)');
axis equal tight
title('q_v - q_l');

%% Уравнение совместности 9

R = computeResidual_9(fields, params, phase, dx, dy);

norm(R)

%% Уравнение 9.Визуализация
figure

pcolor(X, Y, R);
shading interp;
colorbar;
colormap jet;

xlabel('x (м)');
ylabel('y (м)');
axis equal tight
title('Residual 9');
%% Сбор всех норм невязок в один массив

residuals = computeAllResiduals(fields, phase, params, dx, dy, y_interface, x_interface);

total_loss = 0;

for k = 1:length(residuals)
    total_loss = total_loss + norm(residuals{k}(:));
end

fprintf('Общая невязка (сумма L2 норм):  %e\n',total_loss)

%% Добавление наказаний за отклонение от граничных условий (делается через computeTotalResidual) 

% === Упаковка всех переменных в x ===
x0 = packFieldsAndInterface(fields, params.M, params.N, length(params.c_init), params.c_init);

% === Вызов функции потерь ===
loss = computeTotalResidual(x0, params, dx, dy, X, Y, x_interface);

fprintf('Значение функции потерь:  %e\n',loss)


