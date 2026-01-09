function [Lx1, Lx2] = Cost_dual(x1, x2, M)
    [Nx1, Nx2] = size(x1);
    [dx1ds1, dx2ds2] = DCentral(x1, x2, 1/Nx1, 1/Nx2);
    [dx2ds1, dx1ds2] = DCentral(x2, x1, 1/Nx1, 1/Nx2);

    detJ = dx1ds1.*dx2ds2 - dx2ds1.*dx1ds2;
    ds1dx1 = dx2ds2 ./ detJ;
    ds2dx2 = dx1ds1 ./ detJ;
    ds2dx1 = -dx2ds1 ./ detJ;
    ds1dx2 = -dx1ds2 ./ detJ;

    M11 = M.M11(x1, x2);
    M12 = M.M12(x1, x2);
    M22 = M.M22(x1, x2);

    Lx1 = (1/Nx1)^2 * (M11.*ds1dx1.^2+M22.*ds1dx2.^2+2*M12.*ds1dx1.*ds1dx2);
    Lx2 = (1/Nx2)^2 * (M11.*ds2dx1.^2+M22.*ds2dx2.^2+2*M12.*ds2dx1.*ds2dx2);
end