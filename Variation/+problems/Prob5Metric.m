function M = Prob5Metric(x1, x2, cf)
    if nargin < 3 || isempty(cf)
        cf = struct('problemId', 5, 'Nx1', size(x1, 2), 'Nx2', size(x1, 1));
        cf = meshgen.defaults_dual(cf);
    else
        if isfield(cf, 'problem') && ~isfield(cf, 'problemId')
            cf.problemId = cf.problem;
        end
        if ~isfield(cf, 'Nx1') || isempty(cf.Nx1)
            cf.Nx1 = size(x1, 2);
        end
        if ~isfield(cf, 'Nx2') || isempty(cf.Nx2)
            cf.Nx2 = size(x1, 1);
        end
        if isfield(cf, 'problem')
            cf = meshgen.defaults_el(cf);
        else
            cf = meshgen.defaults_dual(cf);
        end
    end

    persistent cached_cf M11_fun M12_fun M22_fun
    if isempty(cached_cf) || ~isequaln(cached_cf, cf)
        [~, ~, M11_fun, M12_fun, M22_fun] = problems.InitProb5(cf);
        cached_cf = cf;
    end

    M.M11 = M11_fun(x1, x2);
    M.M12 = M12_fun(x1, x2);
    M.M22 = M22_fun(x1, x2);

    [M.dM11dx1, M.dM11dx2] = metric_grad(M.M11, x1, x2);
    [M.dM12dx1, M.dM12dx2] = metric_grad(M.M12, x1, x2);
    [M.dM22dx1, M.dM22dx2] = metric_grad(M.M22, x1, x2);
end
