%% полукруг 
theta = linspace(0, pi, 200);
R = 50;
x = R * cos(theta);
y = R * sin(theta);

H = computeCurvatureFromInterface(x, y);

% Ожидаемая кривизна -1/R (для выпуклого полукруга сверху вниз)
expected = -1 / R;

fprintf('Средняя кривизна: %.4e, ожидаемая: %.4e\n', mean(H), expected);

plot(theta, H, 'b-', 'LineWidth', 2);
hold on;
yline(expected, 'r--', 'LineWidth', 2);
xlabel('\theta'); ylabel('Кривизна');
title('Кривизна полукруга');
grid on;
legend('Вычисленная', 'Ожидаемая');


%% эллипс


a = 60; % большая ось
b = 40; % малая ось
theta = linspace(0, pi, 200);
x = a * cos(theta);
y = b * sin(theta);

H = -computeCurvatureFromInterface(x, y);

% Теоретическая кривизна эллипса на полукруге (параметрическая формула):
% k(θ) = (a*b) / ((b*cos(θ))^2 + (a*sin(θ))^2)^(3/2)
expected = (a * b) ./ ((b*cos(theta)).^2 + (a*sin(theta)).^2).^(3/2);

fprintf('Средняя вычисленная кривизна: %.4e\n', mean(H));
fprintf('Средняя теоретическая кривизна: %.4e\n', mean(expected));

figure;
plot(theta, H, 'b-', 'LineWidth', 2);
hold on;
plot(theta, expected, 'r--', 'LineWidth', 2);
xlabel('\theta');
ylabel('Кривизна');
title('Кривизна эллипса');
grid on;
legend('Вычисленная', 'Теоретическая');


%% прямая 

 x = linspace(1,5,500);
 y = 2*x + 3;
 H = computeCurvatureFromInterface(x, y);
expected = zeros(size(x));  % вектор нулей той же длины, что и x


fprintf('Средняя вычисленная кривизна: %.4e\n', mean(H));
fprintf('Теоретическая кривизна: %.4e\n', mean(expected));

figure;
plot(x, H, 'b-', 'LineWidth', 2);
hold on;
plot(x, expected, 'r--', 'LineWidth', 2);
xlabel('x');
ylabel('Кривизна');
title('Кривизна прямой линии');
grid on;
legend('Вычисленная', 'Теоретическая');