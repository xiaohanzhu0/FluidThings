function M = Prob2Metric(x1, x2)
    M.M11 = 1000 + 600*sin(2*pi*x1).*sin(2*pi*x2);
    M.M12 = zeros(size(x1));
    M.M22 = 1000 - 600*sin(2*pi*x1).*sin(2*pi*x2);

    M.dM11dx1 = 1200*pi*cos(2*pi*x1).*sin(2*pi*x2);
    M.dM11dx2 = 1200*pi*sin(2*pi*x1).*cos(2*pi*x2);
    M.dM12dx1 = zeros(size(x1));
    M.dM12dx2 = zeros(size(x1));
    M.dM22dx1 = -1200*pi*cos(2*pi*x1).*sin(2*pi*x2);
    M.dM22dx2 = -1200*pi*sin(2*pi*x1).*cos(2*pi*x2);
end