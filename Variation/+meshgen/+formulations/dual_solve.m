function [x1_samp, x2_samp, hist, info] = dual_solve(params)
%DUAL_SOLVE Dual formulation: harmonic map then metric-conforming inverse map.
    if nargin < 1
        params = struct();
    end

    [x1, x2, M11_fun, M12_fun, M22_fun, boundary_points] = meshgen.problems.init_dual(params);

    if params.plotInitial
        figure
        plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k');
    end

    [x1, x2, harmonic_info] = solve_harmonic(x1, x2, ...
        Omega=params.omega1, ShowPlot=params.showHarmonicPlot, ...
        ShowResPlot=params.showHarmonicRes, BoundaryPoints=boundary_points, ...
        PauseTime=params.pauseTime);

    [params.Nx2, params.Nx1] = size(x1);
    [T1, T2] = ndgrid(linspace(0,1,params.Nx1), linspace(0,1,params.Nx2));
    x1_in_t = griddedInterpolant(T1, T2, x1', params.interpMethod);
    x2_in_t = griddedInterpolant(T1, T2, x2', params.interpMethod);

    [dx1dt1_samp, dx2dt2_samp] = DCentral(x1, x2, 1/(params.Nx1-1), 1/(params.Nx2-1));
    [dx2dt1_samp, dx1dt2_samp] = DCentral(x2, x1, 1/(params.Nx1-1), 1/(params.Nx2-1));

    [T1, T2] = meshgrid(linspace(0,1,params.Nx1), linspace(0,1,params.Nx2));
    t1 = linspace(0,1,50);
    t2 = linspace(0,1,50);
    [t1, t2] = meshgrid(t1, t2);
    dx1dt1 = interp2(T1, T2, dx1dt1_samp, t1, t2, params.interpMethod);
    dx1dt2 = interp2(T1, T2, dx1dt2_samp, t1, t2, params.interpMethod);
    dx2dt1 = interp2(T1, T2, dx2dt1_samp, t1, t2, params.interpMethod);
    dx2dt2 = interp2(T1, T2, dx2dt2_samp, t1, t2, params.interpMethod);

    x1_temp = x1_in_t(t1, t2);
    x2_temp = x2_in_t(t1, t2);
    M11 = M11_fun(x1_temp, x2_temp);
    M12 = M12_fun(x1_temp, x2_temp);
    M22 = M22_fun(x1_temp, x2_temp);

    Mp11 = (dx1dt1.*M11+dx2dt1.*M12).*dx1dt1 + (dx1dt1.*M12+dx2dt1.*M22).*dx2dt1;
    Mp22 = (dx1dt2.*M11+dx2dt2.*M12).*dx1dt2 + (dx1dt2.*M12+dx2dt2.*M22).*dx2dt2;
    Mp12 = (dx1dt1.*M11+dx2dt1.*M12).*dx1dt2 + (dx1dt1.*M12+dx2dt1.*M22).*dx2dt2;

    Mp11_fun = griddedInterpolant(t1', t2', Mp11', params.interpMethod);
    Mp12_fun = griddedInterpolant(t1', t2', Mp12', params.interpMethod);
    Mp22_fun = griddedInterpolant(t1', t2', Mp22', params.interpMethod);

    detJ = Mp11.*Mp22 - Mp12.^2;
    Mp_inv11 = Mp22 ./ detJ;
    Mp_inv22 = Mp11 ./ detJ;
    Mp_inv12 = -Mp12 ./ detJ;

    if params.diagonal_reg == 1
        Mp_inv11 = Mp_inv11 + 0.01*max(Mp_inv11,[],'all');
        Mp_inv22 = Mp_inv22 + 0.01*max(Mp_inv22,[],'all');
    end

    s1 = linspace(0,1,params.Nt1);
    s2 = linspace(0,1,params.Nt2);
    [s1, s2] = meshgrid(s1, s2);
    t1_grid = linspace(0,1,params.Nt1);
    t2_grid = linspace(0,1,params.Nt2);
    [t1_grid, t2_grid] = meshgrid(t1_grid, t2_grid);

    hist = struct('Lt',zeros(1,params.max_iter), ...
                  'Lx',zeros(1,params.max_iter), ...
                  'Theta1',zeros(1,params.max_iter), ...
                  'Theta2',zeros(1,params.max_iter), ...
                  'ThetaInf',zeros(1,params.max_iter));

    warning('off', 'MATLAB:griddedInterpolant:MeshgridEval2DWarnId');
    for i = 1:params.max_iter
        [ds1dt1, ds2dt2] = DCentral(s1, s2, 1/(params.Nt1-1), 1/(params.Nt2-1));
        [ds2dt1, ds1dt2] = DCentral(s2, s1, 1/(params.Nt1-1), 1/(params.Nt2-1));
        J = abs(ds1dt1.*ds2dt2 - ds2dt1.*ds1dt2);

        Mp_inv11_fun = griddedInterpolant(t1_grid', t2_grid', Mp_inv11'.*J', params.interpMethod);
        Mp_inv22_fun = griddedInterpolant(t1_grid', t2_grid', Mp_inv22'.*J', params.interpMethod);
        Mp_inv12_fun = griddedInterpolant(t1_grid', t2_grid', Mp_inv12'.*J', params.interpMethod);

        [s1_new, s2_new] = AssembleLinearSystemDual(t1_grid, t2_grid, ...
            Mp_inv11_fun, Mp_inv12_fun, Mp_inv22_fun);

        s1 = s1 + params.omega2*(s1_new-s1);
        s2 = s2 + params.omega2*(s2_new-s2);

        s1_samp = linspace(0,1,params.sampleNs);
        s2_samp = linspace(0,1,params.sampleNs);
        [s1_samp, s2_samp] = meshgrid(s1_samp, s2_samp);
        t1_samp = griddata(s1, s2, t1_grid, s1_samp, s2_samp);
        t2_samp = griddata(s1, s2, t2_grid, s1_samp, s2_samp);
        x1_samp = x1_in_t(t1_samp, t2_samp);
        x2_samp = x2_in_t(t1_samp, t2_samp);

        if params.doPlot && mod(i, params.plotEvery)==0
            plotIterationMeshes(s1, s2, t1_samp, t2_samp, x1_samp, x2_samp);
            pause(params.pauseTime);
        end

        M11_samp = Mp11_fun(t1_samp, t2_samp);
        M12_samp = Mp12_fun(t1_samp, t2_samp);
        M22_samp = Mp22_fun(t1_samp, t2_samp);
        [Lt1, Lt2] = Cost(t1_samp, t2_samp, M11_samp, M12_samp, M22_samp);
        hist.Lt(i) = mean(Lt1+Lt2,'all');

        M11_samp = M11_fun(x1_samp, x2_samp);
        M12_samp = M12_fun(x1_samp, x2_samp);
        M22_samp = M22_fun(x1_samp, x2_samp);
        [Lx1, Lx2] = Cost(x1_samp, x2_samp, M11_samp, M12_samp, M22_samp);
        hist.Lx(i) = mean(Lx1+Lx2,'all');

        [Theta, Theta_1, Theta_2, Theta_inf] = analysis.skewness(x1_samp, x2_samp);
        hist.Theta1(i) = Theta_1;
        hist.Theta2(i) = Theta_2;
        hist.ThetaInf(i) = Theta_inf;
    end
    warning('on', 'MATLAB:griddedInterpolant:MeshgridEval2DWarnId');

    if nargout > 3
        info.harmonic = harmonic_info;
    end
end

function plotIterationMeshes(s1, s2, t1_samp, t2_samp, x1_samp, x2_samp)
    figure(8)
    plot(s1, s2, 'k'); hold on; plot(s1', s2', 'k'); hold off
    figure(9)
    plot(t1_samp, t2_samp, 'k'); hold on; plot(t1_samp', t2_samp', 'k'); hold off
    figure(10)
    plot(x1_samp, x2_samp, 'k'); hold on; plot(x1_samp', x2_samp', 'k'); hold off
    axis equal
end
