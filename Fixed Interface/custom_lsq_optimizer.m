function [x_opt, history] = custom_lsq_optimizer(fun, x0, lb, ub, options)
% Кастомный аналог lsqnonlin с контролем шага и адаптивным line search
%
% Входные:
%   fun    — вектор-функция невязки
%   x0     — начальное приближение
%   lb, ub — нижние и верхние границы
%   options.step_init — начальный шаг (напр. 1.0)
%   options.step_min  — минимальный шаг (напр. 1e-6)
%   options.max_iter  — максимум итераций
%   options.grad_eps  — шаг для численного градиента

    % Параметры по умолчанию
    if ~isfield(options, 'step_init'), options.step_init = 10.0; end
    if ~isfield(options, 'step_min'),  options.step_min  = 1e-1; end
    if ~isfield(options, 'max_iter'),  options.max_iter  = 100; end
    if ~isfield(options, 'grad_eps'),  options.grad_eps  = 1e-6; end

    step_size = options.step_init;
    step_min = options.step_min;
    max_iter = options.max_iter;
    grad_eps = options.grad_eps;

    x = x0;
    history = struct('fval', [], 'grad_norm', [], 'step', []);
    n = length(x);

    for iter = 1:max_iter
        r = fun(x);
        fval = sum(r.^2);

        % Численный градиент
        grad = zeros(n, 1);
        for i = 1:n
            x_perturb = x;
            x_perturb(i) = x_perturb(i) + grad_eps;
            r_perturb = fun(x_perturb);
            grad(i) = (sum(r_perturb.^2) - fval) / grad_eps;
        end

        grad_norm = norm(grad);

        % Line search (жёсткий нижний порог)
alpha = step_size;
rho = 0.5;
c = 1e-4;
max_ls = 20;

success = false;

for ls_iter = 1:max_ls
    if alpha < step_min
        alpha = step_min;
        fprintf('⚠️ Достигнут минимальный шаг (%.1e), принудительно применяем\n', alpha);
    end

    x_trial = x - alpha * grad;
    x_trial = max(min(x_trial, ub), lb);
    f_trial = sum(fun(x_trial).^2);

    if f_trial < fval - c * alpha * (grad' * grad)
        success = true;
        break;
    end

    if alpha == step_min
        break; % больше уменьшать нельзя — используем как есть
    end

    alpha = max(alpha * rho, step_min);
end

if ~success
    fprintf('⚠️ Line search не нашёл улучшения — идём с минимальным шагом\n');
end
        

        % Обновление
        x = x - alpha * grad;
        x = max(min(x, ub), lb);

        % Логгирование
        fprintf('🔁 Итерация %d: fval = %.4e, ||∇f|| = %.2e, step = %.2e\n', ...
            iter, fval, grad_norm, alpha);

        history.fval(end+1) = fval;
        history.grad_norm(end+1) = grad_norm;
        history.step(end+1) = alpha;

        % Сохраняем прогресс
        save('checkpoint_opt.mat', 'x', 'fval', 'iter', 'grad', 'alpha');

        % Условие остановки
        if grad_norm < 1e-4 || alpha < step_min
            fprintf('✅ Сходимость достигнута на итерации %d\n', iter);
            break;
        end
    end

    x_opt = x;
end
