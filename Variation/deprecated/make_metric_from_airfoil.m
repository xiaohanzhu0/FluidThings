clear
airfoil_datapath = '~/Files/data/Mesh_Generation/Airfoil/foil2/airfoil_18M_coarseIJK.grid';

[x_cv, ~, ~, ~, ~, ~, ~, ~, ~, ~] = readTurtleGrid(airfoil_datapath);

block_idx = 9;

x_metric = x_cv{1,block_idx}(:,:,1,1)';
y_metric = x_cv{1,block_idx}(:,:,1,2)';

[g11, g12, g22] = grid_to_metric(x_metric,y_metric);
g_det = g11.*g22 - g12.^2;
M11 = g22 ./ g_det;
M22 = g11 ./ g_det;
M12 = -g12 ./ g_det;

M.M11 = g11;
M.M12 = g12;
M.M22 = g22;

h = pcolor(x_metric,y_metric,M22);
h.EdgeColor = 'none';

%M.M11 = M11; M.M12 = M12; M.M22 = M22; 

analysis.plot_metric(x_metric, y_metric, M);


function [g11, g12, g22] = grid_to_metric(x1,x2)
    [Nx2, Nx1] = size(x1);
    ds1 = 1 / Nx1;
    ds2 = 1 / Nx2;
    [dx1ds1, dx2ds2] = DCentral(x1, x2, ds1, ds2);
    [dx2ds1, dx1ds2] = DCentral(x2, x1, ds1, ds2);
    
    g11 = dx1ds1.*dx1ds1 + dx2ds1.*dx2ds1;
    g12 = dx1ds1.*dx1ds2 + dx2ds1.*dx2ds2;
    g22 = dx1ds2.*dx1ds2 + dx2ds2.*dx2ds2;
end