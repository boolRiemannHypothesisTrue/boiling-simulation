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

%% УРАВНЕНИЕ НЕПРЕРЫВНОСТИ (ПРОИЗВОДНЫЕ ПО U_X,U_Y)

R_continuity = computeContinuityResidual(fields, params, phase, dx, dy); % невязка уравнения непрерывности

d1 = continuity_du_x(params,phase,R_continuity,dx);
d2 = continuity_du_y(params,phase,R_continuity,dy);

figure;

subplot(1,2,1); % первый подграфик
pcolor(X, Y, d1);
shading interp;
colorbar;
colormap jet;
title('dR/du_x');
xlabel('x (м)');
ylabel('y (м)');
axis equal tight;

subplot(1,2,2); % второй подграфик
pcolor(X, Y, d2);
shading interp;
colorbar;
colormap jet;
title('dR/du_y');
xlabel('x (м)');
ylabel('y (м)');
axis equal tight;

%% УРАВНЕНИЯ НАВЬЕ-СТОКСА/БАЛАНС МОМЕНТА ИМПУЛЬСА (ПРОИЗВОДНЫЕ ПО P,U_X,U_Y)

% d/dP

R_ns_x = computeNSResidualX_full(fields, params, phase, dx, dy); 
R_ns_y = computeNSResidualY_full(fields, params, phase, dx, dy);

dRdP = computePressureDerivatives(R_ns_x,R_ns_y,dx,dy);
 
figure

pcolor(X, Y, dRdP);
shading interp;
colorbar;
colormap jet;

xlabel('x (м)');
ylabel('y (м)');
axis equal tight
title('dR/dP');

% d/du_x для x-ой проекции уравнений Навье-Стокса

