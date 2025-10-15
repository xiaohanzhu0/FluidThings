Nx1 = 40;
Nx2 = 40;

[x1, x2] = InitProb3(Nx1, Nx2, 0.1);
M.M11 = @(x1,x2) ones(Nx1, Nx2);
M.M12 = @(x1,x2) zeros(Nx1, Nx2);
M.M22 = @(x1,x2) ones(Nx1, Nx2);

[Lx1, Lx2] = compute_conformality_field(x1, x2, M);

pcolor(x1, x2, abs(Lx1)+abs(Lx2)); colorbar
figure
pcolor(x1, x2, sqrt(Lx1.^2+Lx2.^2)); colorbar
figure
pcolor(x1, x2, max(Lx1, Lx2)); colorbar
title("Conformality")

function [Lx1, Lx2] = compute_conformality_field(x1, x2, M)
    [Nx1, Nx2] = size(x1);
    [dx1ds1, dx2ds2] = DCentral(x1, x2, 1/Nx1, 1/Nx2);
    [dx2ds1, dx1ds2] = DCentral(x2, x1, 1/Nx1, 1/Nx2);

    M11 = M.M11(x1, x2);
    M12 = M.M12(x1, x2);
    M22 = M.M22(x1, x2);

    Lx1 = (1/Nx1)^2 * (M11.*dx1ds1.^2+M22.*dx2ds1.^2+2*M12.*dx1ds1.*dx2ds1);
    Lx2 = (1/Nx2)^2 * (M22.*dx2ds2.^2+M11.*dx1ds2.^2+2*M12.*dx1ds2.*dx2ds2);
end