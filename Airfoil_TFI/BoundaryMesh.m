clear
BCtype = 1; % 1 for linear, 2 for curved
BCspace = 1; % 1 for equi-spaced, 2 for non equi-spaced
xy = geometry();
bottom = xy{1};
top = xy{2};
left = xy{3};
right = xy{4};

%% Make equispaced data points in i boundary
Nx1 = 50;
Nx2 = size(left, 1);
[bottom, top] = make_boundary_data(left, right, Nx1, BCspace);

if BCtype == 1
    bottom(2:end-1,1) = bottom(2:end-1,1);
elseif BCtype == 2
    bottom(2:end-1,1) = bottom(2:end-1,1) - 0.3*sin(bottom(2:end-1,2)*pi/1.5);
end

figure(1)
plot(bottom(:,1),bottom(:,2),'-o'); hold on
plot(top(:,1),top(:,2),'-o'); hold on
plot(left(:,1),left(:,2)); hold on
plot(right(:,1),right(:,2)); hold on
legend('bottom', 'top', 'left', 'right', Location='best');

%%
function phi_solution = NumBoundaryMesh(F, Nx)
    %dx1dphi = gradient(right(:,1), 1/Nx);
    %dx2dphi = gradient(right(:,2), 1/Nx);
    %F = sqrt(dx1dphi.^2 + dx2dphi.^2);

    c = trapz(1/Nx,F);
    phi_grid = linspace(0, 1, Nx2);
    
    xi_grid = cumtrapz(1/Nx2, F) / c;
    
    xi_desired = linspace(0, 1, 200);
    phi_solution = interp1(xi_grid, phi_grid, xi_desired, 'pchip');
    
    x1_solution = interp1(phi_grid, right(:,1), phi_solution, 'pchip');
    x2_solution = interp1(phi_grid, right(:,2), phi_solution, 'pchip');
    
    plot(x1_solution,x2_solution,'-o');
end
%%
left_dx1dphi = gradient(left(:,1), 1/Nx2);
left_dx2dphi = gradient(left(:,2), 1/Nx2);
left_F = sqrt(left_dx1dphi.^2 + left_dx2dphi.^2);
left_c = trapz(1/Nx2,left_F);

left_phi_grid = linspace(0, 1, 1000);
%%
function [bottom, top] = make_boundary_data(left, right, Nx1, BCspace)
    if BCspace == 1
        Ix_neg = linspace(left(1,1), right(1,1), Nx1);
        Iy_neg = linspace(left(1,2), right(1,2), Nx1);
        Ix_pos = linspace(left(end,1), right(end,1), Nx1);
        Iy_pos = linspace(left(end,2), right(end,2), Nx1);
    elseif BCspace == 2
        cheb = -cos(linspace(0,pi,Nx1)); % Use Chebyshev points for example
        Ix_neg = (left(1,1)+ right(1,1))/2 + (right(1,1)-left(1,1))*cheb/2;
        Iy_neg = (left(1,2)+ right(1,2))/2 + (right(1,2)-left(1,2))*cheb/2;
        Ix_pos = (left(end,1)+ right(end,1))/2 + (right(end,1)-left(end,1))*cheb/2;
        Iy_pos = (left(end,2)+ right(end,2))/2 + (right(end,2)-left(end,2))*cheb/2;
    end
    
    bottom = [Ix_neg', Iy_neg'];
    top = [Ix_pos', Iy_pos'];
end