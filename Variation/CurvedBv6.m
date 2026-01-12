function [res_list, cost_list] = CurvedBv6(cf)
%CURVEDBV6 Legacy entry point for EL formulation (physical coordinates).
    [~, ~, res_list, cost_list] = meshgen.formulations.el_solve(cf);
end
