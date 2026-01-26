clear
% Export harmonic-transformed metric tensor to text files.
if ~exist('output_dir', 'var') || isempty(output_dir)
    output_dir = fullfile(pwd, 'data/hyperbolic_transform_output');
end

if ~exist('params', 'var') || isempty(params)
    params = struct();
end
params = meshgen.defaults_dual(params);

[x1, x2, M11_fun, M12_fun, M22_fun, boundary_points] = meshgen.problems.init_dual(params);
%[x1, x2] = meshgen.formulations.dual_harmonic_map(x1, x2, boundary_points, params);
metric_data = meshgen.formulations.metric_transform_harmonic(x1, x2, M11_fun, M12_fun, M22_fun, params);

if ~isfolder(output_dir)
    mkdir(output_dir);
end

M11 = metric_data.M11;
M12 = metric_data.M12;
M22 = metric_data.M22;
Mp11 = metric_data.Mp11;
Mp12 = metric_data.Mp12;
Mp22 = metric_data.Mp22;
save(fullfile(output_dir, 'x1.txt'), 'x1', '-ascii');
save(fullfile(output_dir, 'x2.txt'), 'x2', '-ascii');
save(fullfile(output_dir, 'M11.txt'), 'M11', '-ascii');
save(fullfile(output_dir, 'M12.txt'), 'M12', '-ascii');
save(fullfile(output_dir, 'M22.txt'), 'M22', '-ascii');
save(fullfile(output_dir, 'Mp11.txt'), 'Mp11', '-ascii');
save(fullfile(output_dir, 'Mp12.txt'), 'Mp12', '-ascii');
save(fullfile(output_dir, 'Mp22.txt'), 'Mp22', '-ascii');

disp(['Saved transformed metric tensors to: ', output_dir]);
%%
% Plot metric components and eigenvectors scaled by sqrt(eigenvalues).
if isvector(x1) && isvector(x2)
    [X1, X2] = meshgrid(x1, x2);
else
    X1 = x1;
    X2 = x2;
end

figure('Name', 'Metric Components', 'Color', 'w');
subplot(1, 3, 1);
pcolor(X1, X2, M11);
shading interp;
axis equal tight;
colorbar;
title('M11');
xlabel('x'); ylabel('y');

subplot(1, 3, 2);
pcolor(X1, X2, M12);
shading interp;
axis equal tight;
colorbar;
title('M12');
xlabel('x'); ylabel('y');

subplot(1, 3, 3);
pcolor(X1, X2, M22);
shading interp;
axis equal tight;
colorbar;
title('M22');
xlabel('x'); ylabel('y');

% Also plot components on a unit-square grid for pipeline checks.
[n_rows, n_cols] = size(Mp11);
xu = linspace(0, 1, n_cols);
yu = linspace(0, 1, n_rows);
[Xu, Yu] = meshgrid(xu, yu);

figure('Name', 'Metric Components (Unit Square)', 'Color', 'w');
subplot(1, 3, 1);
pcolor(Xu, Yu, Mp11);
shading interp;
axis equal tight;
colorbar;
title('Mp11 (unit square)');
xlabel('x'); ylabel('y');

subplot(1, 3, 2);
pcolor(Xu, Yu, Mp12);
shading interp;
axis equal tight;
colorbar;
title('Mp12 (unit square)');
xlabel('x'); ylabel('y');

subplot(1, 3, 3);
pcolor(Xu, Yu, Mp22);
shading interp;
axis equal tight;
colorbar;
title('Mp22 (unit square)');
xlabel('x'); ylabel('y');

% Eigenvectors scaled by sqrt(eigenvalues) for symmetric 2x2 metric.
a = M11;
b = M12;
c = M22;
trace_term = 0.5 * (a + c);
delta = sqrt((0.5 * (a - c)).^2 + b.^2);
lam1 = trace_term + delta;
lam2 = trace_term - delta;

eps_tol = 1e-12;
v1x = b;
v1y = lam1 - a;
degenerate = abs(v1x) < eps_tol & abs(v1y) < eps_tol;
pick_x = degenerate & (a >= c);
pick_y = degenerate & (a < c);
v1x(pick_x) = 1;
v1y(pick_x) = 0;
v1x(pick_y) = 0;
v1y(pick_y) = 1;

v1norm = sqrt(v1x.^2 + v1y.^2);
v1x = v1x ./ (v1norm + eps_tol);
v1y = v1y ./ (v1norm + eps_tol);
v2x = -v1y;
v2y = v1x;

lam1 = max(lam1, 0);
lam2 = max(lam2, 0);
s1 = sqrt(lam1);
s2 = sqrt(lam2);
e1x = v1x .* s1;
e1y = v1y .* s1;
e2x = v2x .* s2;
e2y = v2y .* s2;

step = max(floor(min(size(Mp11)) / 25), 1);
Xs = X1(1:step:end, 1:step:end);
Ys = X2(1:step:end, 1:step:end);
e1x_s = e1x(1:step:end, 1:step:end);
e1y_s = e1y(1:step:end, 1:step:end);
e2x_s = e2x(1:step:end, 1:step:end);
e2y_s = e2y(1:step:end, 1:step:end);

figure('Name', 'Metric Eigenvectors', 'Color', 'w');
quiver(Xs, Ys, e1x_s, e1y_s, 0, 'Color', [0.0 0.45 0.74]);
hold on;
quiver(Xs, Ys, e2x_s, e2y_s, 0, 'Color', [0.85 0.33 0.10]);
hold off;
axis equal tight;
title('Eigenvectors scaled by sqrt(eigenvalues)');
xlabel('x'); ylabel('y');
