%% LoadModelAnalysis.m
% MATLAB analogue of NN4MG/LoadModelAnalysis.ipynb using meshgen solvers.
% This script runs multiple solver configurations (dc/dual/el) and compares
% grid quality metrics (misfit, skewness, cell size ratio).

if isempty(which('meshgen.formulations.dual_solve'))
    startup;
end

plot_ideal_grid = true;
misfit_type = 'standard'; % standard | target1 | target2 | target3
lam_mode = 'fixed';       % fixed | eq9

% ---- Add or remove runs to compare ----
runs_spec = {
    struct('label', 'dc prob5 harmonic-grid', 'solver', 'dc', 'params', struct( ...
        'problemId', 5, 'Nx1', 41, 'Nx2', 41, 'Nt1', 41, 'Nt2', 41, ...
        'sampleNs', 51, 'max_iter', 5, 'doPlot', 0, 'plotInitial', 0, ...
        'showHarmonicPlot', false, 'showHarmonicRes', false, ...
        'with_harmonic_init', true, 'hyperbolic', 1)),
    struct('label', 'dc prob5 hyper-grid', 'solver', 'dc', 'params', struct( ...
        'problemId', 5, 'Nx1', 41, 'Nx2', 41, 'Nt1', 41, 'Nt2', 41, ...
        'sampleNs', 51, 'max_iter', 5, 'doPlot', 0, 'plotInitial', 0, ...
        'showHarmonicPlot', false, 'showHarmonicRes', false, ...
        'with_harmonic_init', false, 'hyperbolic', 1)),
    struct('label', 'dc prob5 og-grid', 'solver', 'dc', 'params', struct( ...
        'problemId', 5, 'Nx1', 41, 'Nx2', 41, 'Nt1', 41, 'Nt2', 41, ...
        'sampleNs', 51, 'max_iter', 5, 'doPlot', 0, 'plotInitial', 0, ...
        'showHarmonicPlot', false, 'showHarmonicRes', false, ...
        'with_harmonic_init', false, 'hyperbolic', 0)),
};

runs = cell(numel(runs_spec), 1);

for k = 1:numel(runs_spec)
    spec = runs_spec{k};
    solver = lower(spec.solver);

    switch solver
        case 'dc'
            params = meshgen.defaults_dual(spec.params);
            [x1, x2, info] = meshgen.formulations.dc_solve(params);
            formulation = 'dual';
            problem_id = params.problemId;
            metric_params = params;
            extra = struct('info', info);
        case 'dual'
            params = meshgen.defaults_dual(spec.params);
            [x1, x2, hist, info] = meshgen.formulations.dual_solve(params);
            formulation = 'dual';
            problem_id = params.problemId;
            metric_params = params;
            extra = struct('hist', hist, 'info', info);
        case 'el'
            cf = meshgen.defaults_el(spec.params);
            [x1, x2, res_list, cost_list, info] = meshgen.formulations.el_solve(cf);
            formulation = 'primal';
            problem_id = cf.problem;
            metric_params = cf;
            extra = struct('res_list', res_list, 'cost_list', cost_list, 'info', info);
        otherwise
            error('Unknown solver: %s', solver);
    end

    [M11, M12, M22] = metric_components(problem_id, x1, x2, metric_params);

    [misfit_field, misfit_avg, misfit_max] = compute_misfit_field(...
        x1, x2, M11, M12, M22, misfit_type, lam_mode);

    skew_field = compute_skewness_field(x1, x2);
    skew_avg = mean(skew_field(:));
    skew_max = max(skew_field(:));

    grad_field = compute_cell_ratio_field(x1, x2);
    grad_avg = mean(grad_field(:));
    grad_max = max(grad_field(:));

    run = struct();
    run.label = spec.label;
    run.solver = solver;
    run.formulation = formulation;
    run.problem_id = problem_id;
    run.X = x1;
    run.Y = x2;
    run.M11 = M11;
    run.M12 = M12;
    run.M22 = M22;
    run.misfit_field = misfit_field;
    run.misfit_avg = misfit_avg;
    run.misfit_max = misfit_max;
    run.skew_field = skew_field;
    run.skew_avg = skew_avg;
    run.skew_max = skew_max;
    run.grad_field = grad_field;
    run.grad_avg = grad_avg;
    run.grad_max = grad_max;
    run.extra = extra;

    if plot_ideal_grid && strcmp(formulation, 'primal')
        [ideal_nx, ideal_ny, X_ideal, Y_ideal] = compute_ideal_grid(x1, x2, M11, M12, M22);
        run.ideal_cells = [ideal_nx, ideal_ny];
        run.X_ideal = X_ideal;
        run.Y_ideal = Y_ideal;
    else
        run.ideal_cells = [];
        run.X_ideal = [];
        run.Y_ideal = [];
    end

    runs{k} = run;

    if strcmp(formulation, 'primal')
        fprintf('[%s] solver=%s problem=%d ideal=%dx%d\n', ...
            run.label, run.solver, run.problem_id, run.ideal_cells(1), run.ideal_cells(2));
    else
        fprintf('[%s] solver=%s problem=%d (dual formulation)\n', ...
            run.label, run.solver, run.problem_id);
    end
end

% ---- Deformed grid ----
num_models = numel(runs);
figure('Name', 'Deformed grid');
for i = 1:num_models
    ax = subplot(1, num_models, i);
    plot_grid(ax, runs{i}.X, runs{i}.Y, runs{i}.label);
end
sgtitle('Deformed grid');

% ---- Ideal grid (primal only) ----
if plot_ideal_grid
    has_primal = any(cellfun(@(r) strcmp(r.formulation, 'primal'), runs));
    if has_primal
        figure('Name', 'Ideal grid');
        for i = 1:num_models
            ax = subplot(1, num_models, i);
            if strcmp(runs{i}.formulation, 'primal')
                label = sprintf('%s (%dx%d)', runs{i}.label, runs{i}.ideal_cells(1), runs{i}.ideal_cells(2));
                plot_grid(ax, runs{i}.X_ideal, runs{i}.Y_ideal, label);
            else
                axis(ax, 'off');
                text(ax, 0.5, 0.5, 'dual formulation', 'HorizontalAlignment', 'center');
            end
        end
        sgtitle('Ideal grid');
    else
        fprintf('Skipping ideal grid: all runs are dual.\n');
    end
end

% ---- Misfit integrand ----
misfit_min = min(cellfun(@(r) min(r.misfit_field(:)), runs));
misfit_max = max(cellfun(@(r) max(r.misfit_field(:)), runs));
figure('Name', 'Misfit integrand');
for i = 1:num_models
    ax = subplot(1, num_models, i);
    plot_scalar_field(ax, runs{i}.X, runs{i}.Y, runs{i}.misfit_field, runs{i}.label, 'viridis', misfit_min, misfit_max);
    add_stats_text(ax, runs{i}.misfit_avg, runs{i}.misfit_max);
end
sgtitle('Misfit integrand');

% ---- Skewness ----
skew_min = min(cellfun(@(r) min(r.skew_field(:)), runs));
skew_max = max(cellfun(@(r) max(r.skew_field(:)), runs));
figure('Name', 'Skewness');
for i = 1:num_models
    ax = subplot(1, num_models, i);
    plot_scalar_field(ax, runs{i}.X, runs{i}.Y, runs{i}.skew_field, runs{i}.label, 'magma', skew_min, skew_max);
    add_stats_text(ax, runs{i}.skew_avg, runs{i}.skew_max, '%.2f');
end
sgtitle('Skewness |90-theta| (deg)');

% ---- Cell size ratio ----
grad_min = min(cellfun(@(r) min(r.grad_field(:)), runs));
grad_max = max(cellfun(@(r) max(r.grad_field(:)), runs));
figure('Name', 'Cell size ratio');
for i = 1:num_models
    ax = subplot(1, num_models, i);
    plot_scalar_field(ax, runs{i}.X, runs{i}.Y, runs{i}.grad_field, runs{i}.label, 'inferno', grad_min, grad_max);
    add_stats_text(ax, runs{i}.grad_avg, runs{i}.grad_max);
end
sgtitle('Cell size ratio');

%% -------- Local helper functions --------
function [M11, M12, M22] = metric_components(problem_id, X, Y, metric_params)
    if nargin >= 4 && ~isempty(metric_params) && problem_id == 5
        M = problems.Prob5Metric(X, Y, metric_params);
    else
        Mfun = meshgen.metrics.for_problem(problem_id);
        M = Mfun(X, Y);
    end
    M11 = M.M11;
    M12 = M.M12;
    M22 = M.M22;
end

function [misfit, avg_val, max_val] = compute_misfit_field(X, Y, M11, M12, M22, misfit_type, lam_mode)
    if nargin < 6
        misfit_type = 'standard';
    end
    if nargin < 7
        lam_mode = 'fixed';
    end
    [n1, n2] = size(X);
    sigma1 = 1 / n1;
    sigma2 = 1 / n2;
    [dx1ds1, dx2ds2] = DCentral(X, Y, sigma1, sigma2);
    [dx2ds1, dx1ds2] = DCentral(Y, X, sigma1, sigma2);

    q_xi = M11 .* dx1ds1.^2 + M22 .* dx2ds1.^2 + 2 .* M12 .* dx1ds1 .* dx2ds1;
    q_eta = M11 .* dx1ds2.^2 + M22 .* dx2ds2.^2 + 2 .* M12 .* dx1ds2 .* dx2ds2;

    if strcmp(lam_mode, 'eq9')
        lam_xi = sigma_from_q(q_xi);
        lam_eta = sigma_from_q(q_eta);
    else
        lam_xi = 1.0;
        lam_eta = 1.0;
    end

    misfit = misfit_terms(q_xi, q_eta, lam_xi, lam_eta, misfit_type);
    avg_val = mean(misfit(:));
    max_val = max(misfit(:));
end

function sigma = sigma_from_q(q)
    num = mean(abs(q(:)));
    den = mean(q(:).^2);
    sigma = sqrt(num / max(den, 1e-12));
end

function integrand = misfit_terms(q_xi, q_eta, lam_xi, lam_eta, misfit_type)
    switch lower(misfit_type)
        case 'standard'
            integrand = (lam_xi.^2) .* q_xi + (lam_eta.^2) .* q_eta;
        case 'target1'
            integrand = (lam_xi.^2 .* q_xi - 1).^2 + (lam_eta.^2 .* q_eta - 1).^2;
        case 'target2'
            integrand = abs((q_xi - mean(q_xi(:))) ./ mean(q_xi(:))) + ...
                        abs((q_eta - mean(q_eta(:))) ./ mean(q_eta(:)));
        case 'target3'
            aux_xi = (q_xi - mean(q_xi(:))) ./ mean(q_xi(:));
            aux_eta = (q_eta - mean(q_eta(:))) ./ mean(q_eta(:));
            integrand = softplus(aux_xi) + softplus(-aux_xi) + ...
                        softplus(aux_eta) + softplus(-aux_eta) - 4 * log(2);
        otherwise
            error('Unknown misfit_type: %s', misfit_type);
    end
end

function y = softplus(x)
    y = log1p(exp(-abs(x))) + max(x, 0);
end

function skew = compute_skewness_field(X, Y)
    eps_val = 1e-12;
    dX_dxi = zeros(size(X));
    dY_dxi = zeros(size(Y));
    dX_deta = zeros(size(X));
    dY_deta = zeros(size(Y));

    dX_dxi(:, 2:end-1) = 0.5 * (X(:, 3:end) - X(:, 1:end-2));
    dY_dxi(:, 2:end-1) = 0.5 * (Y(:, 3:end) - Y(:, 1:end-2));
    dX_dxi(:, 1) = X(:, 2) - X(:, 1);
    dY_dxi(:, 1) = Y(:, 2) - Y(:, 1);
    dX_dxi(:, end) = X(:, end) - X(:, end-1);
    dY_dxi(:, end) = Y(:, end) - Y(:, end-1);

    dX_deta(2:end-1, :) = 0.5 * (X(3:end, :) - X(1:end-2, :));
    dY_deta(2:end-1, :) = 0.5 * (Y(3:end, :) - Y(1:end-2, :));
    dX_deta(1, :) = X(2, :) - X(1, :);
    dY_deta(1, :) = Y(2, :) - Y(1, :);
    dX_deta(end, :) = X(end, :) - X(end-1, :);
    dY_deta(end, :) = Y(end, :) - Y(end-1, :);

    v1x = dX_dxi; v1y = dY_dxi;
    v2x = dX_deta; v2y = dY_deta;

    dotp = v1x .* v2x + v1y .* v2y;
    n1 = sqrt(v1x .* v1x + v1y .* v1y);
    n2 = sqrt(v2x .* v2x + v2y .* v2y);
    cosang = dotp ./ (n1 .* n2 + eps_val);
    cosang = min(max(cosang, -1.0), 1.0);
    angle = acosd(cosang);
    skew = abs(90 - angle);
end

function grad = compute_cell_ratio_field(X, Y)
    eps_val = 1e-12;
    x00 = X(1:end-1, 1:end-1);
    y00 = Y(1:end-1, 1:end-1);
    x10 = X(1:end-1, 2:end);
    y10 = Y(1:end-1, 2:end);
    x11 = X(2:end, 2:end);
    y11 = Y(2:end, 2:end);
    x01 = X(2:end, 1:end-1);
    y01 = Y(2:end, 1:end-1);

    area = 0.5 * abs(...
        x00 .* y10 - y00 .* x10 + ...
        x10 .* y11 - y10 .* x11 + ...
        x11 .* y01 - y11 .* x01 + ...
        x01 .* y00 - y01 .* x00);
    area = max(area, eps_val);

    ratio = ones(size(area));
    if size(area, 1) > 1
        r = max(area(2:end, :) ./ area(1:end-1, :), area(1:end-1, :) ./ area(2:end, :));
        ratio(1:end-1, :) = max(ratio(1:end-1, :), r);
        ratio(2:end, :) = max(ratio(2:end, :), r);
    end
    if size(area, 2) > 1
        r = max(area(:, 2:end) ./ area(:, 1:end-1), area(:, 1:end-1) ./ area(:, 2:end));
        ratio(:, 1:end-1) = max(ratio(:, 1:end-1), r);
        ratio(:, 2:end) = max(ratio(:, 2:end), r);
    end

    accum = zeros(size(X));
    count = zeros(size(X));
    accum(1:end-1, 1:end-1) = accum(1:end-1, 1:end-1) + ratio;
    accum(1:end-1, 2:end) = accum(1:end-1, 2:end) + ratio;
    accum(2:end, 1:end-1) = accum(2:end, 1:end-1) + ratio;
    accum(2:end, 2:end) = accum(2:end, 2:end) + ratio;

    count(1:end-1, 1:end-1) = count(1:end-1, 1:end-1) + 1;
    count(1:end-1, 2:end) = count(1:end-1, 2:end) + 1;
    count(2:end, 1:end-1) = count(2:end, 1:end-1) + 1;
    count(2:end, 2:end) = count(2:end, 2:end) + 1;

    grad = accum ./ max(count, 1);
end

function [ideal_nx, ideal_ny, X_ideal, Y_ideal] = compute_ideal_grid(X, Y, M11, M12, M22)
    [n1, n2] = size(X);
    sigma1 = 1 / n1;
    sigma2 = 1 / n2;
    [dx1ds1, dx2ds2] = DCentral(X, Y, sigma1, sigma2);
    [dx2ds1, dx1ds2] = DCentral(Y, X, sigma1, sigma2);

    p1 = M11 .* dx1ds1.^2 + M22 .* dx2ds1.^2 + 2 .* M12 .* dx1ds1 .* dx2ds1;
    p2 = M11 .* dx1ds2.^2 + M22 .* dx2ds2.^2 + 2 .* M12 .* dx1ds2 .* dx2ds2;

    sig1 = sqrt(sum(abs(p1(:))) / max(sum(p1(:).^2), 1e-12));
    sig2 = sqrt(sum(abs(p2(:))) / max(sum(p2(:).^2), 1e-12));

    ideal_nx = max(1, round(1 / sig1));
    ideal_ny = max(1, round(1 / sig2));

    s1 = linspace(0, 1, ideal_nx + 1);
    s2 = linspace(0, 1, ideal_ny + 1);
    [s1_grid, s2_grid] = meshgrid(s1, s2);

    s1_src = linspace(0, 1, size(X, 2));
    s2_src = linspace(0, 1, size(X, 1));
    [s1_src_grid, s2_src_grid] = meshgrid(s1_src, s2_src);

    X_ideal = interp2(s1_src_grid, s2_src_grid, X, s1_grid, s2_grid, 'linear');
    Y_ideal = interp2(s1_src_grid, s2_src_grid, Y, s1_grid, s2_grid, 'linear');
end

function plot_grid(ax, X, Y, title_str)
    plot(ax, X, Y, 'k', 'LineWidth', 0.1);
    hold(ax, 'on');
    plot(ax, X', Y', 'k', 'LineWidth', 0.1);
    axis(ax, 'equal');
    if nargin > 3 && ~isempty(title_str)
        title(ax, title_str);
    end
    hold(ax, 'off');
end

function plot_scalar_field(ax, X, Y, field, title_str, cmap, vmin, vmax)
    pcolor(ax, X, Y, field);
    shading(ax, 'flat');
    axis(ax, 'equal');
    if nargin > 5 && ~isempty(cmap)
        try
            colormap(ax, cmap);
        catch
            colormap(ax, 'parula');
        end
    end
    if nargin > 6 && ~isempty(vmin) && ~isempty(vmax)
        caxis(ax, [vmin vmax]);
    end
    if nargin > 4 && ~isempty(title_str)
        title(ax, title_str);
    end
    colorbar(ax);
end

function add_stats_text(ax, avg_val, max_val, fmt)
    if nargin < 4 || isempty(fmt)
        fmt = '%.3e';
    end
    txt = sprintf(['avg=' fmt '\nmax=' fmt], avg_val, max_val);
    text(ax, 0.02, 0.98, txt, 'Units', 'normalized', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'white', 'Margin', 2);
end
