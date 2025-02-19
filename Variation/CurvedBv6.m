clear
addpath('./','./utils')
animation = 1;
pause_time = 0.2;
make_gif = 0;
problem = 4;
method = 1;
gif_name = 'example13.gif';
title_name = 'Curved boundary, Curved metric, Alternative Cost function';
save_output = 0;

C = 0.1;
Nx1 = 80;
Nx2 = 40;
max_iter = 100;
tolerance = 1e-6;
epsilon = 0.2; % Under-relaxation factor

% Grid size
N = Nx1*Nx2;

% Equispaced computational coordinates
s1 = linspace(0, 1, Nx1);
s2 = linspace(0, 1, Nx2);

% Initial physical coordinates (equispaced)
if ismember(problem, [1,2])
    [x1, x2] = meshgrid(s1, s2);
elseif problem == 3
    [x1, x2] = InitProb3(Nx1, Nx2, C);
elseif problem == 4
    [x1, x2] = InitProb4(Nx1, Nx2);
    [x1_exact, x2_exact] = InitProb4(50, 50);
end

% Plot initial grid
plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k');
title(title_name); xlim([-0.1,1.1]); ylim([-0.1,1.1]);
axis equal; hold off
if make_gif; exportgraphics(gcf, gif_name); end

res_list = [];
%%
for iter = 1:max_iter
    [A, b, res] = AssembleLinearSystem(x1, x2, problem, method);
    x_star = A \ b;

    x1_star = x_star(1:N);
    x2_star = x_star(N+1:end);

    dx1 = reshape(x1_star, Nx2, Nx1) - x1;
    dx2 = reshape(x2_star, Nx2, Nx1) - x2;

    res_list = [res_list, res];
    x1 = x1 + epsilon * dx1;
    x2 = x2 + epsilon * dx2;

    if problem == 4
        [x1, x2] = UpdateCorrection(x1, x2, x1_exact, x2_exact);
    end
    
    if animation == 1
        plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k');
        title(title_name); xlim([-0.1,1.1]); ylim([-0.1,1.1]);
        axis equal; hold off;
        pause(pause_time);
    end
    if make_gif; exportgraphics(gcf, gif_name, Append=true); end
end

%%
figure()
plot(1:iter, res_list);
legend('Cost', 'Residual');
grid on
