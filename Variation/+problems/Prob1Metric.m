function M = Prob1Metric(x1, x2)
    M.M11 = 40000 * (1 + 15 * x1).^(-2);
    M.M22 = 40000 * (1 + 15 * x2).^(-2);
    M.M12 = zeros(size(x1));

    M.dM11dx1 = -1200000 * (1 + 15 * x1).^(-3);
    M.dM11dx2 = zeros(size(x1));
    M.dM22dx1 = zeros(size(x1));
    M.dM22dx2 = -1200000 * (1 + 15 * x2).^(-3);
    M.dM12dx1 = zeros(size(x1));
    M.dM12dx2 = zeros(size(x1));
end