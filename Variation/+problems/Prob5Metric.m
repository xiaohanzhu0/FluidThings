function M = Prob5Metric(x1, x2)
    [Nx2, Nx1] = size(x1);
    M.M11 = 400*ones(Nx2, Nx1);
    M.M22 = 400*ones(Nx2, Nx1);
    M.M12 = zeros(Nx2, Nx1);

    M.dM11dx1 = zeros(Nx2, Nx1);
    M.dM11dx2 = zeros(Nx2, Nx1);
    M.dM22dx1 = zeros(Nx2, Nx1);
    M.dM22dx2 = zeros(Nx2, Nx1);
    M.dM12dx1 = zeros(Nx2, Nx1);
    M.dM12dx2 = zeros(Nx2, Nx1);
end