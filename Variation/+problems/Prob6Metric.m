function M = Prob6Metric(x1, x2)
    theta = pi/6;
    s1 = x1*cos(theta) + x2*sin(theta);
    s2 = -x1*sin(theta) + x2*cos(theta);

    M11_temp = 40000 * (1 + C*15 * s1).^(-2);
    M22_temp = 40000 * (1 + C*15 * s2).^(-2);
    M12_temp = zeros(size(s1));

    M.M11 = cos(theta)*(M11_temp*cos(theta)-M12_temp*sin(theta)) - sin(theta)*(M12_temp*cos(theta)-M22_temp*sin(theta));
    M.M22 = sin(theta)*(M11_temp*sin(theta)+M12_temp*cos(theta)) + cos(theta)*(M12_temp*sin(theta)+M22_temp*cos(theta));
    M.M12 = cos(theta)*(M11_temp*sin(theta)+M12_temp*cos(theta)) - sin(theta)*(M12_temp*sin(theta)+M22_temp*cos(theta));

    [M.dM11dx1, M.dM11dx2] = metric_grad(M.M11, x1, x2);
    [M.dM22dx1, M.dM22dx2] = metric_grad(M.M22, x1, x2);
    [M.dM12dx1, M.dM12dx2] = metric_grad(M.M12, x1, x2);
end