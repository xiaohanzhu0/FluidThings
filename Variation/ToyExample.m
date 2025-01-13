clear
% MATLAB code to solve a 2D nonlinear elliptic PDE using fixed-point iteration and finite differences

% Problem: Solve -div(a(u) grad(u)) = f in a unit square domain with Dirichlet boundary conditions

% Parameters
N = 50;                   % Number of grid points in each direction
h = 1/(N+1);              % Grid spacing
x = linspace(0, 1, N+2);  % Grid points including boundaries
y = linspace(0, 1, N+2);
omega = 10;               % Ralaxization parameter

% Source term f(x, y) and nonlinear coefficient a(u)
f = @(x, y) sin(pi*x).*sin(pi*y);    % Example source function
a = @(u) 1 + u.^2;                   % Nonlinear coefficient

% Initialization
u = zeros(N+2, N+2);      % Solution initialized to zero
u_new = zeros(N+2, N+2);  % Placeholder for updates

% Boundary conditions
u(1, :) = 0;              % u = 0 on bottom boundary
u(:, 1) = 0;              % u = 0 on left boundary
u(end, :) = 0;            % u = 0 on top boundary
u(:, end) = 0;            % u = 0 on right boundary

% Iteration parameters
tol = 1e-6;               % Convergence tolerance
max_iter = 1000;          % Maximum iterations
iter = 0;
err = 100;
min_err = Inf;
err_list = [];

% Fixed-point iteration
while err > tol && iter < max_iter && err < min_err
    min_err = err;
    for i = 2:N+1
        for j = 2:N+1
            % Compute coefficients based on current values of u
            aE = a(0.5 * (u(i+1, j) + u(i, j)));
            aW = a(0.5 * (u(i-1, j) + u(i, j)));
            aN = a(0.5 * (u(i, j+1) + u(i, j)));
            aS = a(0.5 * (u(i, j-1) + u(i, j)));

            % Fixed-point update based on finite differences
            u_new(i, j) = (1/(aE + aW + aN + aS)) * ( ...
                aE * u(i+1, j) + aW * u(i-1, j) + ...
                aN * u(i, j+1) + aS * u(i, j-1) - h^2 * f(x(i), y(j)));
        end
    end
    
    % Compute error and update
    err = max(max(abs(u_new - u)));
    err_list = [err_list, err];
    %u = u_new;
    u = u + omega*(u_new - u);
    iter = iter + 1;
end

% Display results
fprintf('Converged in %d iterations with error %.2e\n', iter, err);
[X, Y] = meshgrid(x, y);
figure()
surf(X, Y, u);
xlabel('x'); ylabel('y'); zlabel('u(x,y)'); title('2D Nonlinear Elliptic PDE Solution');

%%
figure();
plot(1:iter, err_list);
xlabel('Iteration'); ylabel('Residual');
