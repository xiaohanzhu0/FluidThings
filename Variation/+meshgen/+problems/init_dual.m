function [x1, x2, M11_fun, M12_fun, M22_fun, boundary_points] = init_dual(params)
%INIT_DUAL Initialize physical grid and metric components for dual formulation.
    [x1, x2, M11_fun, M12_fun, M22_fun] = problems.Initialization(...
        params.problemId, params.Nx1, params.Nx2, params);
    boundary_points = extract_boundary_points(x1, x2);
end
