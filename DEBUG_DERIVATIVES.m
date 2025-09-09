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

%  ∑∑(Q_x^(i1,j1)  (∂Q_x^(i1,j1))/(∂P^(i,j) )+Q_y^(i1,j1)  (∂Q_y^(i1,j1))/(∂P^(i,j) )) 

Qx = computeNSResidualX_full(fields, params, phase, dx, dy); 
Qy = computeNSResidualY_full(fields, params, phase, dx, dy);

a  = Qx_dQx_dp_PLUS_Qy_dQy_dp(Qx,Qy,dx,dy);
 
figure

pcolor(X, Y, a);
shading interp;
colorbar;
colormap jet;

xlabel('x (м)');
ylabel('y (м)');
axis equal tight
title('Qx*dQxdp + Qy*dQydp');

norm(a)

%% ∑∑〖Q_x^(i1,j1)  (∂Q_x^(i1,j1))/(∂u_x^(i,j) )〗



Qx = computeNSResidualX_full(fields, params, phase, dx, dy); 
a = Qx_dQx_dux(Qx,fields,params,phase,dx,dy);

 
figure

pcolor(X, Y, a);
shading interp;
colorbar;
colormap jet;

xlabel('x (м)');
ylabel('y (м)');
axis equal tight
title('Qx * dQxdux');

norm(a)

%% ∑∑〖Q_x^(i1,j1)  (∂Q_x^(i1,j1))/(∂u_y^(i,j) )〗



Qx = computeNSResidualX_full(fields, params, phase, dx, dy); 
a = Qx_dQx_duy(Qx,fields,params,phase,dx,dy);

 
figure

pcolor(X, Y, a);
shading interp;
colorbar;
colormap jet;

xlabel('x (м)');
ylabel('y (м)');
axis equal tight
title('Qx * dQxduy');

norm(a)

%% ∑∑〖Q_y^(i1,j1)  (∂Q_y^(i1,j1))/(∂u_x^(i,j) )〗



Qy = computeNSResidualY_full(fields, params, phase, dx, dy); 
a = Qy_dQy_dux(Qy,fields,params,phase,dx,dy);

 
figure

pcolor(X, Y, a);
shading interp;
colorbar;
colormap jet;

xlabel('x (м)');
ylabel('y (м)');
axis equal tight
title('Qx * dQxduy');

norm(a)


%% ∑∑〖Q_y^(i1,j1)  (∂Q_y^(i1,j1))/(∂u_x^(i,j) )〗



Qy = computeNSResidualY_full(fields, params, phase, dx, dy); 
a = Qy_dQy_duy(Qy, fields,params,phase,dx,dy);

 
figure

pcolor(X, Y, a);
shading interp;
colorbar;
colormap jet;

xlabel('x (м)');
ylabel('y (м)');
axis equal tight
title('Qx * dQxduy');

norm(a)


%% ∑∑〖Q_T^(i1,j1)  (∂Q_T^(i1,j1))/(∂T^(i,j) )〗



QT = computeEnergyResidual(fields, params, phase, dx, dy); 
a = QT_dQT_dT(QT, fields,params,phase,dx,dy);

 
figure

pcolor(X, Y, a);
shading interp;
colorbar;
colormap jet;

xlabel('x (м)');
ylabel('y (м)');
axis equal tight
title('Qx * dQxduy');

norm(a)


%% ∑∑〖Q_T^(i1,j1)  (∂Q_T^(i1,j1))/(∂u_x^(i,j) )〗

QT = computeEnergyResidual(fields, params, phase, dx, dy); 
a = QT_dQT_dux(QT, fields,params,phase,dx,dy);

 
figure

pcolor(X, Y, a);
shading interp;
colorbar;
colormap jet;

xlabel('x (м)');
ylabel('y (м)');
axis equal tight
title('Qx * dQxduy');

norm(a)