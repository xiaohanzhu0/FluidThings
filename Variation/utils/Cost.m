function L = Cost(x1, x2, sigma1, sigma2, problem)
    [Nx1, Nx2] = size(x1);

    [M11, M22] = M(x1, x2, problem);
    [dx1ds1, dx2ds2] = DCentral(x1, x2, sigma1, sigma2);
    [dx2ds1, dx1ds2] = DCentral(x2, x1, sigma1, sigma2);
    L1 = sigma1^2 * (M11.*dx1ds1.^2+M22.*dx2ds1.^2);
    L2 = sigma2^2 * (M22.*dx2ds2.^2+M11.*dx1ds2.^2);
    L = norm(L1(:)+L2(:), 1) / (Nx1*Nx2);
end