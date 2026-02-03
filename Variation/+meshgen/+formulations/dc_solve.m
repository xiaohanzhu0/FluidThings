function [x1, x2, info] = dc_solve(params)
%DUAL_SOLVE Dual formulation: harmonic map then metric-conforming inverse map.
    if nargin < 1
        params = struct();
    end

    [x1, x2, M11_fun, M12_fun, M22_fun, boundary_points] = meshgen.problems.init_dual(params);


    if params.with_harmonic_init
        [x1, x2, harmonic_info] = meshgen.formulations.dual_harmonic_map(x1, x2, boundary_points, params);
    else 
        harmonic_info = 0;
    end

    metric_data = meshgen.formulations.metric_transform_harmonic( ...
        x1, x2, M11_fun, M12_fun, M22_fun, params);

    [t1, t2] = meshgen.formulations.dc_map(metric_data);
    [t1_grid, t2_grid] = meshgrid(t1,t2);

    x1 = metric_data.x1_in_t(t1_grid, t2_grid);
    x2 = metric_data.x2_in_t(t1_grid, t2_grid);

    if nargout > 2
        info.harmonic = harmonic_info;
    end
end
