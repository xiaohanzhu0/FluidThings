function [x1, x2, res_list, cost_list, info] = el_solve(cf)
%EL_SOLVE Euler-Lagrange solver in physical coordinates.
    if nargin < 1
        cf = struct();
    end

    [x1, x2, boundary_points, cf] = meshgen.problems.init_physical(cf);

    if cf.make_gif
        plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k');
        title(cf.title_name);
        xlim([min(x1,[],'all'),max(x1,[],'all')]);
        ylim([min(x2,[],'all'),max(x2,[],'all')]);
        axis equal; hold off
        exportgraphics(gcf, cf.gif_name);
    end

    res_list = [];
    err_list = [];
    cost_list = [];

    cf.forDx = 0;
    Mfun = meshgen.metrics.for_problem(cf.problem);
    [A, b] = AssembleLinearSystem(x1, x2, Mfun, cf);
    res = norm(A*[x1(:);x2(:)] - b);

    err = [x1(:);x2(:)];
    res_list = [res_list, norm(res)];

    disp('Started iterations');
    iter = 0;
    for iter = 1:cf.max_iter
        [dx1, dx2, res, err, cf] = nonlinear_step(x1, x2, Mfun, cf, err);

        dx1(1,1)=0; dx1(1,end)=0; dx1(end,1)=0; dx1(end,end)=0;
        dx2(1,1)=0; dx2(1,end)=0; dx2(end,1)=0; dx2(end,end)=0;

        x1 = x1 + cf.omega*dx1;
        x2 = x2 + cf.omega*dx2;

        res_list = [res_list, norm(res)];
        err_list = [err_list, norm(norm(dx1,'fro'),norm(dx2,'fro'))];

        if should_apply_boundary(cf.problem)
            [x1, x2] = UpdateCorrection(x1, x2, boundary_points);
        end

        if cf.animation == 1 && ~mod(iter, cf.iter_per_frame)
            plot_iteration(x1, x2, cf);
        end

        if cf.plot_res == 1 && ~mod(iter, cf.iter_per_frame)
            figure(6)
            semilogy(res_list); grid on
            title([cf.title_name, ' residual']); hold off
        end

        if cf.tolerance > res_list(end)
            break;
        end

        if ~mod(iter, cf.iter_per_frame)
            [isValid, ~] = checkGridOverlap(x1, x2);
            if ~isValid; disp('Warning: cell overlap expected'); end
        end
    end

    if cf.save_output
        save_all(x1, x2, res_list, cost_list);
    end

    if nargout > 4
        info.iterations = iter;
        info.err_list = err_list;
        info.boundary_points = boundary_points;
    end
end

function tf = should_apply_boundary(problemId)
    tf = ismember(problemId, [0, 3, 4, 5, 8]);
end

function [dx1, dx2, res, err, cf] = nonlinear_step(x1, x2, Mfun, cf, err)
    [Nx2, Nx1] = size(x1);
    N = Nx1 * Nx2;

    if cf.nonlinear == 4
        cf.forDx = 1;
        [A, b] = AssembleLinearSystem(x1, x2, Mfun, cf);
        res = A*err - b;
        err = A \ b;
        res = abs(res(1:N)) + abs(res(N+1:end));
        res = reshape(res, cf.Nx2, cf.Nx1);
        res = res(2:end-1, 2:end-1);

        dx1 = reshape(err(1:N), cf.Nx2, cf.Nx1);
        dx2 = reshape(err(N+1:end), cf.Nx2, cf.Nx1);
    elseif cf.nonlinear == 5
        [A, b] = AssembleLinearSystemApprox(x1, x2, Mfun, cf);
        res = A*err - b;
        err = A \ b;
        res = abs(res(1:N)) + abs(res(N+1:end));
        res = reshape(res, cf.Nx2, cf.Nx1);
        res = res(2:end-1, 2:end-1);

        dx1 = reshape(err(1:N), cf.Nx2, cf.Nx1) - x1;
        dx2 = reshape(err(N+1:end), cf.Nx2, cf.Nx1) - x2;
    elseif cf.nonlinear == 6
        [dx, res] = AssembleLinearSystemGMRE(x1, x2, Mfun, cf);
        dx1 = reshape(dx(1:N), cf.Nx2, cf.Nx1);
        dx2 = reshape(dx(N+1:end), cf.Nx2, cf.Nx1);
    elseif cf.nonlinear == 1 || cf.nonlinear == 2
        [A, b, res] = AssembleLinearSystem(x1, x2, Mfun, cf);
        x_new = A \ b;
        x1_new = reshape(x_new(1:N), cf.Nx2, cf.Nx1);
        x2_new = reshape(x_new(N+1:end), cf.Nx2, cf.Nx1);

        dx1 = x1_new - x1;
        dx2 = x2_new - x2;

        res = abs(res(1:N)) + abs(res(N+1:end));
        res = reshape(res, cf.Nx2, cf.Nx1);
        res = res(2:end-1, 2:end-1);
    elseif cf.nonlinear == 7
        [A, b, res] = AssembleLinearSystemConserve(x1, x2, Mfun, cf);
        x_new = A \ b;
        x1_new = reshape(x_new(1:N), cf.Nx2, cf.Nx1);
        x2_new = reshape(x_new(N+1:end), cf.Nx2, cf.Nx1);

        dx1 = x1_new - x1;
        dx2 = x2_new - x2;

        res = abs(res(1:N)) + abs(res(N+1:end));
        res = reshape(res, cf.Nx2, cf.Nx1);
        res = res(2:end-1, 2:end-1);
    else
        error('Unsupported nonlinear solver id: %d', cf.nonlinear);
    end
end

function plot_iteration(x1, x2, cf)
    figure(2)
    plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k');
    xlim([min(x1,[],'all'),max(x1,[],'all')]);
    ylim([min(x2,[],'all'),max(x2,[],'all')]);
    title(cf.title_name, ' full mesh'); hold off;
    if cf.make_gif; exportgraphics(gcf, cf.gif_name, Append=true); end

    figure(3)
    plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k');
    xlim([-0.5,1.5]);
    ylim([-1,1]);
    title(cf.title_name, ' close look 1'); hold off;

    figure(4)
    plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k');
    xlim([-0.02,0.02]);
    ylim([-0.02,0.02]);
    title(cf.title_name, ' close look 2'); hold off;

    figure(5)
    plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k');
    xlim([0.9,1.1]);
    ylim([-0.1,0.1]); hold off;
    title(cf.title_name, ' close look 3'); hold off;
    pause(cf.pause_time);
end
