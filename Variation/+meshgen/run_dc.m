function [x1_samp, x2_samp, info] = run_dc(params)
%RUN_DUAL Entry point for dual/harmonic-map formulation.
    if nargin < 1
        params = struct();
    end
    params = meshgen.defaults_dual(params);
    [x1_samp, x2_samp, info] = meshgen.formulations.dc_solve(params);
end
