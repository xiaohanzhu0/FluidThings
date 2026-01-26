function [x1, x2, info] = dual_harmonic_map(x1, x2, boundary_points, params)
%DUAL_HARMONIC_MAP Compute harmonic map for the dual formulation.
    if params.plotInitial
        figure
        plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k');
    end

    [x1, x2, info] = solve_harmonic(x1, x2, ...
        Omega=params.omega1, ShowPlot=params.showHarmonicPlot, ...
        ShowResPlot=params.showHarmonicRes, BoundaryPoints=boundary_points, ...
        PauseTime=params.pauseTime);
end
