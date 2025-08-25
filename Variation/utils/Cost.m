function [L1, L2, Linf, gtilde2, sigma1, sigma2] = Cost(x1, x2, M_type, C)
    [Nx1, Nx2] = size(x1);

    M = GetM(x1, x2, M_type, C);
    [dx1ds1, dx2ds2] = DCentral(x1, x2, 1/Nx1, 1/Nx2);
    [dx2ds1, dx1ds2] = DCentral(x2, x1, 1/Nx1, 1/Nx2);

    Lx1 = (1/Nx1)^2 * (M.M11.*dx1ds1.^2+M.M22.*dx2ds1.^2+2*M.M12.*dx1ds1.*dx2ds1);
    Lx2 = (1/Nx2)^2 * (M.M22.*dx2ds2.^2+M.M11.*dx1ds2.^2+2*M.M12.*dx1ds2.*dx2ds2);

    Lx = abs(Lx1) + abs(Lx2);
    L1 = norm(Lx(:), 1) / (Nx1*Nx2);

    Lx = sqrt(Lx1.^2 + Lx2.^2);
    L2 = norm(Lx(:), 2) / sqrt(Nx1*Nx2);

    Lx = max(Lx1, Lx2);
    Linf = norm(Lx(:), 'inf');

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