function [x1_samp, x2_samp, hist] = dual_metric_inverse_map(metric_data, M11_fun, M12_fun, M22_fun, params)
%DUAL_METRIC_INVERSE_MAP Metric conforming inverse map after harmonic map.
    x1_in_t = metric_data.x1_in_t;
    x2_in_t = metric_data.x2_in_t;
    Mp11_fun = metric_data.Mp11_fun;
    Mp12_fun = metric_data.Mp12_fun;
    Mp22_fun = metric_data.Mp22_fun;
    Mp_inv11 = metric_data.Mp_inv11;
    Mp_inv12 = metric_data.Mp_inv12;
    Mp_inv22 = metric_data.Mp_inv22;

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
