function [x1_samp, x2_samp, hist, info] = run_dual(params)
%RUN_DUAL Entry point for dual/harmonic-map formulation.
    if nargin < 1
        params = struct();
    end
    params = meshgen.defaults_dual(params);
    [x1_samp, x2_samp, hist, info] = meshgen.formulations.dual_solve(params);
end
