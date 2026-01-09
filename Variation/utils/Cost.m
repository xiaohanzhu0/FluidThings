function [Lx1, Lx2] = Cost(x1, x2, M11_samp, M12_samp, M22_samp)
    [Nx1, Nx2] = size(x1);
    [dx1ds1, dx2ds2] = DCentral(x1, x2, 1/Nx1, 1/Nx2);
    [dx2ds1, dx1ds2] = DCentral(x2, x1, 1/Nx1, 1/Nx2);

    M11 = M11_samp;
    M12 = M12_samp;
    M22 = M22_samp;

    Lx1 = (1/(Nx1+1))^2 * (M11.*dx1ds1.^2+M22.*dx2ds1.^2+2*M12.*dx1ds1.*dx2ds1);
    Lx2 = (1/(Nx2+1))^2 * (M22.*dx2ds2.^2+M11.*dx1ds2.^2+2*M12.*dx1ds2.*dx2ds2);
end