function [x1_samp, x2_samp, hist, info] = dual_solve(params)
%DUAL_SOLVE Dual formulation: harmonic map then metric-conforming inverse map.
    if nargin < 1
        params = struct();
    end

    [x1, x2, M11_fun, M12_fun, M22_fun, boundary_points] = meshgen.problems.init_dual(params);

    [x1, x2, harmonic_info] = meshgen.formulations.dual_harmonic_map( ...
        x1, x2, boundary_points, params);

    metric_data = meshgen.formulations.metric_transform_harmonic( ...
        x1, x2, M11_fun, M12_fun, M22_fun, params);
    
    metric_data = meshgen.formulations.metric_transform_inverse(metric_data, params);

    [x1_samp, x2_samp, hist] = meshgen.formulations.dual_metric_inverse_map( ...
        metric_data, M11_fun, M12_fun, M22_fun, params);

    if nargout > 3
        info.harmonic = harmonic_info;
    end
end
