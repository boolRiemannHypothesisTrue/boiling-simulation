options = struct();
options.step_init = 10.0;
options.step_min = 1e-1;
options.max_iter = 100;
options.grad_eps = 1e-6;

[x_opt, history] = custom_lsq_optimizer(fun, x0, lb, ub, options);
