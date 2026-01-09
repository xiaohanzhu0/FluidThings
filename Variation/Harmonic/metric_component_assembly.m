function g = metric_component_assembly(g11, g12, g22)
    [Nx, Ny] = size(g11);
    g = zeros(2, 2, Nx, Ny);
    g(1,1,:,:) = reshape(g11, [1 1 Nx Ny]);
    g(1,2,:,:) = reshape(g12, [1 1 Nx Ny]);
    g(2,1,:,:) = reshape(g12, [1 1 Nx Ny]);
    g(2,2,:,:) = reshape(g22, [1 1 Nx Ny]);
end
