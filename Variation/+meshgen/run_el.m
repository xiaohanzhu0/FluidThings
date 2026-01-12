function [x1, x2, res_list, cost_list, info] = run_el(cf)
%RUN_EL Entry point for EL formulation in physical coordinates.
    if nargin < 1
        cf = struct();
    end
    cf = meshgen.defaults_el(cf);
    [x1, x2, res_list, cost_list, info] = meshgen.formulations.el_solve(cf);
end
