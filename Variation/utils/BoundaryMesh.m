% Used for symbolically parametrized curve (x1(s), x2(s)) with 0<=s<=1,
% or for numerically parametrized curve as a vector
% F is weight*norm(dx1/ds, dx2/ds),
% Nx is number of desired boundary points.
function phi_solution = BoundaryMesh(F, Nx)
    phi_grid = linspace(0, 1, length(F));
    if isa(F,'double')
        c = trapz(1/Nx,F);
        xi_grid = cumtrapz(1/Nx, F) / c;
    else
        c = integral(F, 0, 1);
        F_values = F(phi_grid);
        xi_grid = cumtrapz(phi_grid, F_values) / c;
    end
    xi_desired = linspace(0, 1, Nx); % Desired xi points
    phi_solution = interp1(xi_grid, phi_grid, xi_desired, 'pchip');

end