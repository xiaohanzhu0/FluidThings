function [x1_samp, x2_samp] = nonlinear_dual(params)
%NONLINEAR_DUAL Legacy entry point for dual/harmonic-map formulation.
    if nargin < 1
        params = struct();
    end
    params = meshgen.defaults_dual(params);
    [x1_samp, x2_samp] = meshgen.formulations.dual_solve(params);
end
