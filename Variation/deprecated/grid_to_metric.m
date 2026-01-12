[x1, x2] = InitProb8(80, 40);

[Nx1, Nx2] = size(x1);
ds1 = 1 / Nx1;
ds2 = 1 / Nx2;


function [g11, g12, g22] = grid_to_metric(x1,x2)
    [Nx1, Nx2] = size(x1);
    ds1 = 1 / Nx1;
    ds2 = 1 / Nx2;
    [dx1ds1, dx2ds2] = DCentral(x1, x2, ds1, ds2);
    [dx2ds1, dx1ds2] = DCentral(x2, x1, ds1, ds2);
    
    g11 = dx1ds1.*dx1ds1 + dx2ds1.*dx2ds1;
    g12 = dx1ds1.*dx1ds2 + dx2ds1.*dx2ds2;
    g22 = dx1ds2.*dx1ds2 + dx2ds2.*dx2ds2;
end
