function metric_data = metric_transform_harmonic(x1, x2, M11_fun, M12_fun, M22_fun, params)
%METRIC_TRANSFORM_HARMONIC Compute the harmonic-transformed metric Mp.
    [Nx2, Nx1] = size(x1);
    [T1, T2] = ndgrid(linspace(0,1,Nx1), linspace(0,1,Nx2));
    metric_data.x1_in_t = griddedInterpolant(T1, T2, x1', params.interpMethod);
    metric_data.x2_in_t = griddedInterpolant(T1, T2, x2', params.interpMethod);

    [dx1dt1_samp, dx2dt2_samp] = DCentral(x1, x2, 1/(Nx1-1), 1/(Nx2-1));
    [dx2dt1_samp, dx1dt2_samp] = DCentral(x2, x1, 1/(Nx1-1), 1/(Nx2-1));

    [T1, T2] = meshgrid(linspace(0,1,Nx1), linspace(0,1,Nx2));
    t1 = linspace(0,1,size(x1,2));
    t2 = linspace(0,1,size(x1,1));
    [t1, t2] = meshgrid(t1, t2);
    dx1dt1 = interp2(T1, T2, dx1dt1_samp, t1, t2, params.interpMethod);
    dx1dt2 = interp2(T1, T2, dx1dt2_samp, t1, t2, params.interpMethod);
    dx2dt1 = interp2(T1, T2, dx2dt1_samp, t1, t2, params.interpMethod);
    dx2dt2 = interp2(T1, T2, dx2dt2_samp, t1, t2, params.interpMethod);

    x1_temp = metric_data.x1_in_t(t1, t2);
    x2_temp = metric_data.x2_in_t(t1, t2);
    M11 = M11_fun(x1_temp, x2_temp);
    M12 = M12_fun(x1_temp, x2_temp);
    M22 = M22_fun(x1_temp, x2_temp);

    metric_data.M11 = M11;
    metric_data.M12 = M12;
    metric_data.M22 = M22;

    metric_data.Mp11 = (dx1dt1.*M11+dx2dt1.*M12).*dx1dt1 + (dx1dt1.*M12+dx2dt1.*M22).*dx2dt1;
    metric_data.Mp22 = (dx1dt2.*M11+dx2dt2.*M12).*dx1dt2 + (dx1dt2.*M12+dx2dt2.*M22).*dx2dt2;
    metric_data.Mp12 = (dx1dt1.*M11+dx2dt1.*M12).*dx1dt2 + (dx1dt1.*M12+dx2dt1.*M22).*dx2dt2;

    metric_data.Mp11_fun = griddedInterpolant(t1', t2', metric_data.Mp11', params.interpMethod);
    metric_data.Mp12_fun = griddedInterpolant(t1', t2', metric_data.Mp12', params.interpMethod);
    metric_data.Mp22_fun = griddedInterpolant(t1', t2', metric_data.Mp22', params.interpMethod);
end
