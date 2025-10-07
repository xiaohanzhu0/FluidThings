function [L, gtilde2, sigma1, sigma2] = CostExact(x1, x2, Mfun)
    [Nx1, Nx2] = size(x1);

    M = Mfun(x1,x2);
    [dx1ds1, dx2ds2] = DCentral(x1, x2, 1/Nx1, 1/Nx2);
    [dx2ds1, dx1ds2] = DCentral(x2, x1, 1/Nx1, 1/Nx2);

    m1 = (1/Nx1)^2 * (M.M11.*dx1ds1.^2+M.M22.*dx2ds1.^2+2*M.M12.*dx1ds1.*dx2ds1) - 1;
    m2 = (1/Nx2)^2 * (M.M22.*dx2ds2.^2+M.M11.*dx1ds2.^2+2*M.M12.*dx1ds2.*dx2ds2) - 1;

    L = m1.^2 + m2.^2;
    L = norm(L(:),1) / (Nx1*Nx2);

    g12 = dx1ds1.*dx1ds2 + dx2ds1.*dx2ds2;
    g11 = dx1ds1.*dx1ds1 + dx2ds1.*dx2ds1;
    g22 = dx1ds2.*dx1ds2 + dx2ds2.*dx2ds2;
    gtilde = g12 ./ sqrt(g11 .* g22);
    gtilde2 = norm(gtilde(:));

    %% Explore with ideal number of mesh points
    p1 = M.M11.*dx1ds1.^2 + M.M22.*dx2ds1.^2 + 2*M.M12.*dx1ds1.*dx2ds1;
    p2 = M.M11.*dx1ds2.^2 + M.M22.*dx2ds2.^2 + 2*M.M12.*dx1ds2.*dx2ds2;
    sigma1 = sqrt( sum(abs(p1),'all') / sum(p1.^2,'all') );
    sigma2 = sqrt( sum(abs(p2),'all') / sum(p2.^2,'all') );
end