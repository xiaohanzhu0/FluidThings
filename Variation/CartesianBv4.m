% Second version:  diagonal metric tensor, m11 only depends on x1,
%                  m22 only depends on x2, unit square boundary,
%                  grid points x1 and x2 coordinates are solved in ONE
%                  SYSTEM
clear

addpath('./','./utils')
animation = 1;
pause_time = 0.2;
make_gif = 0;
problem = 1;
method = 1;
gif_name = 'example1.gif';
title_name = 'Cartesian Boundary';

max_iter = 150;
tolerance = 1e-6;
epsilon = 0.2; % Under-relaxation factor

% Grid size
Nx1 = 30;
Nx2 = 30;
N = Nx1*Nx2;

% Equispaced computational coordinates
s1 = linspace(0, 1, Nx1);
s2 = linspace(0, 1, Nx2);
[s1, s2] = meshgrid(s1, s2);

% Initial physical coordinates (equispaced)
x1 = s1;
x2 = s2;

% Plot initial grid
plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k');
title(title_name); xlim([0,1]); ylim([0,1]);
axis equal; hold off
if make_gif; exportgraphics(gcf, gif_name); end

[M11, M22] = M(x1, x2, problem);

min_cost = Inf;

for iter = 1:max_iter
    [A, b] = AssembleLinearSystem(x1, x2, s1, s2, problem, method);
    x_star = A \ b;

    x1_star = x_star(1:N);
    x2_star = x_star(N+1:end);

    dx1 = reshape(x1_star, Nx2, Nx1) - x1;
    dx2 = reshape(x2_star, Nx2, Nx1) - x2;

    % Check convergence
    cost = Cost(x1, x2, s1, s2, problem);
    %if cost > min_cost
    %    break
    %end

    min_cost = cost;
    x1 = x1 + epsilon * dx1;
    x2 = x2 + epsilon * dx2;

    plot(x1, x2, 'k'); hold on; plot(x1', x2', 'k');
    title(title_name); xlim([0,1]); ylim([0,1]);
    axis equal; hold off;
    pause(pause_time);
    if make_gif; exportgraphics(gcf, gif_name, Append=true); end
end


%%
function [A, b] = AssembleLinearSystem(x1, x2, s1, s2, problem, method)
    [Nx1, Nx2] = size(x1);
    N = Nx1*Nx2;
    N_all = 2*N;

    sigma1 = 1 / Nx1;
    sigma2 = 1 / Nx2;
    [M11, M22] = M(x1, x2, problem);
    [dM11dx1, dM11dx2, dM22dx1, dM22dx2] = dMdx(x1, x2, problem);
    [dx1ds1, dx2ds2] = DCentral(x1, x2, s1, s2);
    [dx2ds1, dx1ds2] = DCentral(x2, x1, s1, s2);

    A = sparse(N_all, N_all);
    id = GetIndex(Nx1, Nx2);
    
    if method == 2 % For approximate cost function
        coef1 = -2*sigma1^2*M11.*dx1ds1.*M11.*dx1ds1;
        for i = id.inner
            A(i,i) = -2;
            A(i,i-Nx2) = 1;
            A(i,i+Nx2) = 1;
        end
    
        coef2 = -2*sigma2^2*M22.*dx2ds2.*M22.*dx2ds2;
        for i = N+id.inner
            A(i,i) = -2;
            A(i,i-1) = 1;
            A(i,i+1) = 1;
        end

        for i = N+[id.l, id.r]  % s2 component of left and right boundary: replace Laciancian by dy
            A(i,i) = -2;
            A(i,i-1) = 1;
            A(i,i+1) = 1;
        end
    
        for i = [id.b, id.t]  % s1 component of top and bottom boundary: replace Laciancian by dx
            A(i,i) = -2;
            A(i,i-Nx2) = 1;
            A(i,i+Nx2) = 1;
        end
    
        coef = [coef1(:); coef2(:)];
        A = A.*coef;
        % For Dirichlet conditions
        for i =   [id.l, id.r, id.corner]; A(i, i) = 1; end
        for i = N+[id.b, id.t, id.corner]; A(i, i) = 1; end
    
        b1 = sigma1^4 * M11 .* dx1ds1 .* (dM11dx1 .* dx1ds1 .* dx1ds1 .* dx1ds1);
        b2 = sigma2^4 * M22 .* dx2ds2 .* (dM22dx2 .* dx2ds2 .* dx2ds2 .* dx2ds2);

    elseif method == 1 % For alternative cost function
        for i = [id.inner, N+id.inner]  % N+id.inner for the second component
            A(i,i) = -4;
            A(i,i-1) = 1;
            A(i,i+1) = 1;
            A(i,i-Nx2) = 1;
            A(i,i+Nx2) = 1;
        end
    
        for i = N+[id.l, id.r]  % s2 component of left and right boundary: replace Laciancian by dy
            A(i,i) = -2;
            A(i,i-1) = 1;
            A(i,i+1) = 1;
        end
    
        for i = [id.b, id.t]  % s1 component of top and bottom boundary: replace Laciancian by dx
            A(i,i) = -2;
            A(i,i-Nx2) = 1;
            A(i,i+Nx2) = 1;
        end
        Mii = [M11(:); M22(:)];
        A = -2*A.*Mii;
    
        % For Dirichlet conditions
        for i =   [id.l, id.r, id.corner]; A(i, i) = 1; end
        for i = N+[id.b, id.t, id.corner]; A(i, i) = 1; end
    
        % Assemble the vector b
        b1 = dM11dx1.* ((dx1ds1.^2)*sigma1^2 + (dx1ds2.^2)*sigma2^2);
        b1 = b1 + 2*dM11dx2.*(dx2ds1.*dx1ds1*sigma1^2 + dx2ds2.*dx1ds2*sigma2^2);
        b1 = b1 - dM22dx1.*(dx2ds1.^2*sigma1^2 + dx2ds2.^2*sigma2^2);
    
        b2 = dM22dx2.* ((dx2ds1.^2)*sigma1^2 + (dx2ds2.^2)*sigma2^2);
        b2 = b2 + 2*dM22dx1.*(dx1ds1.*dx2ds1*sigma1^2 + dx1ds2.*dx2ds2*sigma2^2);
        b2 = b2 - dM11dx2.*(dx1ds1.^2*sigma1^2 + dx1ds2.^2*sigma2^2);
    end
    % Apply x1=0 and x1=1 to the left and right boundary
    b1(:,1) = 0;
    b1(:,end) = 1;
    % Apply x2=0 and x2=1 to the top and bottom boundary
    b2(1,:) = 0;
    b2(end,:) = 1;
    b = [b1(:); b2(:)];
end
