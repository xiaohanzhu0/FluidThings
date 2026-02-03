function [x1_s, x2_s] = dc_map(metric_data)
%DC_MAP Decoupled metric-conforming map using averaged diagonal metric terms.
%   The metric is averaged to M11(x1) and M22(x2), then the inverse maps
%   s1(x1), s2(x2) are built and inverted to produce x1(s1), x2(s2).

    %#ok<*NASGU> % M11_fun, M22_fun, params kept for interface consistency.

    M11 = metric_data.Mp11;
    M22 = metric_data.Mp22;

    [Nx2, Nx1] = size(M11);
    %t1 = linspace(0, 1, Nx1);
    %t2 = linspace(0, 1, Nx2);
    %[t1_grid, t2_grid] = meshgrid(t1, t2);
    %x1 = metric_data.x1_in_t(t1_grid, t2_grid);
    %x2 = metric_data.x2_in_t(t1_grid, t2_grid);
    x1 = linspace(0, 1, Nx1);
    x2 = linspace(0, 1, Nx2);
    dx1 = x1(2) - x1(1);
    dx2 = x2(2) - x2(1);

    % Average M11 over x2 and M22 over x1.
    M11_bar = mean(M11, 1, 'omitnan');
    M22_bar = mean(M22, 2, 'omitnan');

    sqrtM11 = sqrt(max(M11_bar, eps));
    sqrtM22 = sqrt(max(M22_bar, eps));

    I1 = trapz(dx1, sqrtM11);
    I2 = trapz(dx2, sqrtM22);
    if I1 <= 0
        I1 = eps;
    end
    if I2 <= 0
        I2 = eps;
    end
    Delta1 = 1 / I1;
    Delta2 = 1 / I2;

    s1_of_x1 = Delta1 * cumtrapz(dx1, sqrtM11);
    s2_of_x2 = Delta2 * cumtrapz(dx2, sqrtM22);

    s1_of_x1 = s1_of_x1 - s1_of_x1(1);
    s2_of_x2 = s2_of_x2 - s2_of_x2(1);
    if s1_of_x1(end) ~= 0
        s1_of_x1 = s1_of_x1 / s1_of_x1(end);
    end
    if s2_of_x2(end) ~= 0
        s2_of_x2 = s2_of_x2 / s2_of_x2(end);
    end

    N1 = max(round(1 / Delta1), 2);
    N2 = max(round(1 / Delta2), 2);

    s1 = linspace(0, 1, N1);
    s2 = linspace(0, 1, N2);
    x1_s = interp1(s1_of_x1, x1, s1, 'linear');
    x2_s = interp1(s2_of_x2, x2, s2, 'linear');

    x1_s = x1_s(:).';
    x2_s = x2_s(:).';
end
